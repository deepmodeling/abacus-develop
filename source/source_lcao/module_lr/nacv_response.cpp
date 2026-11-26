/** @file nacv_response.cpp
 * @brief Build and reuse the real-Gamma NAC orbital response for one geometry.
 */
#include "nacv_response.h"
#include "gradient_inputs.h"
#include "hamilt_zeq_l.h"
#include "hamilt_zeq_r.h"
#include "nacv_solver.h"
#include "nacv_tools.h"
#include "zeq_solver.h"
#include <algorithm>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace LR
{
struct NacvResponse::Impl
{
    const GradientInputs<double> &inputs;
    const int channel;
    const bool open_shell;
    const std::string solver;
    const NacvResponseControls controls;
    const std::string input_directory;
    std::unique_ptr<Z_vector_L<double>> closed_left;
    std::unique_ptr<Z_vector_UL<double>> open_left;
    std::unique_ptr<Z_vector_R<double>> closed_right;
    std::unique_ptr<Z_vector_UR<double>> open_right;
    Parallel_Blkcyclic matrix_layout;
    std::unique_ptr<NacvFactorization> factorization;

    Impl(
        const GradientInputs<double> &data,
        int selected_channel,
        bool unrestricted,
        const std::string &selected_solver,
        const NacvResponseControls &settings)
        : inputs(data), channel(selected_channel), open_shell(unrestricted), solver(selected_solver), controls(settings)
    {
    }

    void prepare_left()
    {
        const NacvTimerScope timer("NacvResponse", "prepare_left");
        const GradientInputs<double> &g = inputs;
        if (open_shell)
        {
            if (open_left)
            {
                return;
            }
            const int dimension = g.nk * (g.nocc[0] * g.nvirt[0] + g.nocc[1] * g.nvirt[1]);
            open_left.reset(new Z_vector_UL<double>(
                g.xc_kernel,
                g.nspin,
                g.nbasis,
                g.nocc,
                g.nvirt,
                g.ucell,
                g.orb_cutoff,
                g.gd,
                g.psi_ks,
                g.eig_ks,
#ifdef __EXX
                g.exx_lri,
                g.hybrid_alpha,
#endif
                g.pot_hxc_gs,
                g.kv,
                g.px,
                g.pc,
                g.pmat,
                g.dft_functional,
                g.nk,
                g.nloc,
                dimension));
        }
        else
        {
            if (closed_left)
            {
                return;
            }
            const std::string &spin = g.spin_types[channel];
            closed_left.reset(new Z_vector_L<double>(
                g.xc_kernel,
                g.nspin,
                g.nbasis,
                g.nocc,
                g.nvirt,
                g.ucell,
                g.orb_cutoff,
                g.gd,
                g.psi_ks,
                g.eig_ks,
#ifdef __EXX
                g.exx_lri,
                g.hybrid_alpha,
#endif
                g.pot_hxc_gs,
                g.kv,
                g.px,
                g.pc,
                g.pmat,
                spin,
                input_directory,
                g.out_dir,
                g.dft_functional));
        }
    }

    void prepare_right()
    {
        const NacvTimerScope timer("NacvResponse", "prepare_right");
        const GradientInputs<double> &g = inputs;
        const std::string &spin = g.spin_types[channel];
        const std::weak_ptr<PotHxcLR> potential = g.pot[channel];
        if (open_shell)
        {
            if (open_right)
            {
                return;
            }
            const int dimension = g.nk * (g.nocc[0] * g.nvirt[0] + g.nocc[1] * g.nvirt[1]);
            open_right.reset(new Z_vector_UR<double>(
                g.xc_kernel,
                g.nspin,
                g.nbasis,
                g.nocc,
                g.nvirt,
                g.ucell,
                g.orb_cutoff,
                g.gd,
                g.psi_ks,
                g.eig_ks,
#ifdef __EXX
                g.exx_lri,
                g.hybrid_alpha,
#endif
                potential,
                g.pot_hxc_gs,
                g.kv,
                g.px,
                g.pc,
                g.pmat,
                g.ks_solver,
                g.dft_functional,
                g.nk,
                g.nloc,
                dimension));
        }
        else
        {
            if (closed_right)
            {
                return;
            }
            closed_right.reset(new Z_vector_R<double>(
                g.xc_kernel,
                g.nspin,
                g.nbasis,
                g.nocc,
                g.nvirt,
                g.ucell,
                g.orb_cutoff,
                g.gd,
                g.psi_ks,
                g.eig_ks,
#ifdef __EXX
                g.exx_lri,
                g.hybrid_alpha,
#endif
                potential,
                g.pot_hxc_gs,
                g.kv,
                g.px,
                g.pc,
                g.pmat,
                input_directory,
                g.out_dir,
                g.dft_functional,
                spin));
        }
    }

    template <typename Operator> void prepare_matrix(Operator &left, int spin_count)
    {
        const NacvTimerScope timer("NacvResponse", "prepare_matrix");
        const int dimension = zvec_global_dim(left, spin_count);
        if (solver == "lapack")
        {
            std::vector<double> matrix = left.matrix();
            factorization.reset(new NacvFactorization(std::move(matrix), dimension));
        }
#ifdef __MPI
        else
        {
            const Parallel_Blkcyclic &px = inputs.px[0];
            const int process_rows = px.nrproc();
            const int process_columns = px.ncproc();
            const int process_dimension = std::max(process_rows, process_columns);
            const int local_dimension = dimension / process_dimension;
            const int bounded_dimension = std::min(32, local_dimension);
            const int block = std::max(1, bounded_dimension);
            LR_Util::setup_2d_division(
                matrix_layout,
                block,
                dimension,
                dimension,
                px.ictxt());
            std::vector<double> matrix = zvec_hessian_2d<double>(
                left,
                spin_count,
                inputs.nloc,
                matrix_layout);
            const bool cholesky = solver == "scalapack_chol";
            factorization.reset(new NacvFactorization(std::move(matrix), matrix_layout, cholesky));
        }
#else
        else
        {
            throw std::logic_error("ScaLAPACK NAC response requires an MPI build");
        }
#endif
    }

    template <typename Operator>
    void solve_with(
        Operator &left,
        double *rhs,
        double *solution,
        int count,
        int spin_count)
    {
        const int ld = inputs.nloc;
        if (solver == "cg")
        {
            std::function<void(const double *, double *)> apply = [&left, ld, count](const double *in, double *out)
            {
                left.hPsi(
                    in,
                    out,
                    ld,
                    count);
            };
            solve_nac_cg(
                solution,
                rhs,
                ld,
                count,
                std::move(apply),
                inputs.test_force,
                controls.tolerance,
                controls.max_iterations);
            return;
        }
        if (!factorization)
        {
            prepare_matrix(left, spin_count);
        }
        const int dimension = zvec_global_dim(left, spin_count);
        const std::size_t full_size = static_cast<std::size_t>(dimension) * count;
        std::vector<double> full(full_size, 0.0);
        for (int state = 0; state < count; ++state)
        {
            double *target = full.data() + static_cast<std::size_t>(state) * dimension;
            const double *source = rhs + state * ld;
#ifdef __MPI
            zvec_local_to_full(
                left,
                spin_count,
                source,
                target);
#else
            std::copy_n(source, dimension, target);
#endif
        }
        if (solver == "lapack")
        {
            factorization->solve(full, count);
        }
#ifdef __MPI
        else
        {
            Parallel_Blkcyclic rhs_layout;
            LR_Util::setup_2d_division(
                rhs_layout,
                matrix_layout.nb(),
                dimension,
                count,
                matrix_layout.ictxt());
            const std::size_t storage = std::max<std::size_t>(1, rhs_layout.size());
            std::vector<double> local(storage, 0.0);
            for (int col = 0; col < rhs_layout.ncloc(); ++col)
            {
                for (int row = 0; row < rhs_layout.nrloc(); ++row)
                {
                    local[col * rhs_layout.nrloc() + row] =
                        full[rhs_layout.icglb(col) * dimension + rhs_layout.irglb(row)];
                }
            }
            factorization->solve(local, rhs_layout);
            std::fill(full.begin(), full.end(), 0.0);
            LR_Util::gather_2d_to_full(
                rhs_layout,
                local.data(),
                full.data(),
                false,
                dimension,
                count);
        }
#endif
        for (int state = 0; state < count; ++state)
        {
            const double *source = full.data() + static_cast<std::size_t>(state) * dimension;
            double *target = solution + state * ld;
#ifdef __MPI
            zvec_full_to_local(
                left,
                spin_count,
                source,
                target);
#else
            std::copy_n(source, dimension, target);
#endif
        }
        if (inputs.test_force)
        {
            LR_Util::print_value(solution, count, ld);
        }
    }
};

NacvResponse::NacvResponse(
    const GradientInputs<double> &inputs,
    int channel,
    bool open_shell,
    const std::string &solver,
    const NacvResponseControls &controls)
{
    impl_.reset(new Impl(
        inputs,
        channel,
        open_shell,
        solver,
        controls));
}
NacvResponse::~NacvResponse() = default;

void NacvResponse::solve(double *rhs, double *solution, int count)
{
    const NacvTimerScope timer("NacvResponse", "solve");
    if (!impl_->closed_left && !impl_->open_left)
    {
        impl_->prepare_left();
    }
    if (impl_->open_shell)
    {
        impl_->solve_with(
            *impl_->open_left,
            rhs,
            solution,
            count,
            2);
    }
    else
    {
        impl_->solve_with(
            *impl_->closed_left,
            rhs,
            solution,
            count,
            1);
    }
}
void NacvResponse::solve_polarized(const double *amplitudes, double *solution, int count)
{
    const NacvTimerScope timer("NacvResponse", "solve_polarized");
    if (!impl_->closed_right && !impl_->open_right)
    {
        impl_->prepare_right();
    }
    const int ld = impl_->inputs.nloc;
    ct::Tensor rhs = LR_Util::newTensor<double>({count, ld});
    rhs.zero();
    if (impl_->open_shell)
    {
        impl_->open_right->hPsi(
            amplitudes,
            rhs.data<double>(),
            ld,
            count);
    }
    else
    {
        impl_->closed_right->hPsi(
            amplitudes,
            rhs.data<double>(),
            ld,
            count);
    }
    solve(rhs.data<double>(), solution, count);
}
} // namespace LR
