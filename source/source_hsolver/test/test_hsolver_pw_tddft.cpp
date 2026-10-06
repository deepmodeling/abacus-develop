#include "source_base/parallel_global.h"
#include "source_hsolver/hsolver_pw_tddft.h"
#include "source_hsolver/test/linear_test_utils.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <gtest/gtest.h>
#include <sstream>
#include <type_traits>

namespace
{
using Device = base_device::DEVICE_CPU;

template <typename T>
class DenseHamiltonian final : public hsolver::HSOperator<T>
{
  public:
    std::vector<std::vector<T>> matrices;
    std::vector<int> dimensions;

    explicit DenseHamiltonian(const std::vector<int>& sizes) : dimensions(sizes)
    {
        for (const int n: sizes)
        {
            matrices.push_back(linear_test::hermitian<T>(n));
        }
    }
    void update_k(const int ik) override
    {
        current_k_ = ik;
    }
    void hpsi(const T* x, T* y, const int ld, const int nvec) const override
    {
        const int n = dimensions.at(current_k_);
        const std::vector<T>& h = matrices.at(current_k_);
        for (int band = 0; band < nvec; ++band)
        {
            for (int i = 0; i < n; ++i)
            {
                y[band * ld + i] = T(0);
                for (int j = 0; j < n; ++j)
                {
                    y[band * ld + i] += h[j * n + i] * x[band * ld + j];
                }
            }
        }
    }
    void spsi(const T* x, T* y, const int ld, const int nvec) const override
    {
        FAIL() << "PW CN propagation must not request an overlap operator";
    }

  private:
    int current_k_ = -1;
};

void initialize_basis(const std::vector<int>& sizes, ModulePW::PW_Basis_K* basis)
{
    basis->nks = sizes.size();
    basis->npwk_max = *std::max_element(sizes.begin(), sizes.end());
    basis->tpiba = 0.8;
    basis->tpiba2 = 0.64;
    // These arrays are owned and released by PW_Basis/PW_Basis_K.
    basis->npwk = new int[basis->nks];
    basis->kvec_c = new ModuleBase::Vector3<double>[basis->nks];
    const int size = basis->nks * basis->npwk_max;
    basis->gcar = new ModuleBase::Vector3<double>[size];
    basis->gk2 = new double[size]();
    for (int ik = 0; ik < basis->nks; ++ik)
    {
        basis->npwk[ik] = sizes[ik];
        basis->kvec_c[ik] = ModuleBase::Vector3<double>(0.1 * ik, 0.07, -0.03);
        for (int i = 0; i < sizes[ik]; ++i)
        {
            const int index = ik * basis->npwk_max + i;
            basis->gcar[index] = ModuleBase::Vector3<double>(0.2 * i, 0.1 * (i % 3), -0.15);
            const ModuleBase::Vector3<double> q = basis->gcar[index] + basis->kvec_c[ik];
            basis->gk2[index] = q.norm2();
        }
    }
}

template <typename T>
std::vector<linear_test::Complex> dense_cn(const std::vector<T>& h, const std::vector<T>& previous, int n, int bands, double dt)
{
    std::vector<T> lhs(n * n);
    std::vector<T> rhs_matrix(n * n);
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            const T identity = i == j ? T(1) : T(0);
            lhs[j * n + i] = identity + T(0, dt / 4.0) * h[j * n + i];
            rhs_matrix[j * n + i] = identity - T(0, dt / 4.0) * h[j * n + i];
        }
    }
    const std::vector<T> rhs = linear_test::multiply(rhs_matrix, previous, n, bands);
    return linear_test::lapack_solve(lhs, rhs, n, bands);
}

template <typename T>
void check_dense_step(const DenseHamiltonian<T>& op, const psi::Psi<T>& previous, const psi::Psi<T>& current, const double dt)
{
    const double tolerance = std::is_same<T, std::complex<float>>::value ? 1e-4 : 1e-10;
    const int bands = current.get_nbands();
    for (int ik = 0; ik < current.get_nk(); ++ik)
    {
        previous.fix_k(ik);
        current.fix_k(ik);
        const int n = op.dimensions[ik];
        std::vector<T> packed(n * bands);
        for (int band = 0; band < bands; ++band)
        {
            for (int i = 0; i < n; ++i)
            {
                packed[band * n + i] = previous(band, i);
            }
        }
        const std::vector<linear_test::Complex> reference = dense_cn(op.matrices[ik], packed, n, bands, dt);
        for (int band = 0; band < bands; ++band)
        {
            for (int i = 0; i < n; ++i)
            {
                EXPECT_LT(std::abs(linear_test::Complex(current(band, i)) - reference[band * n + i]), tolerance);
            }
            for (int other = 0; other < bands; ++other)
            {
                T overlap(0);
                for (int i = 0; i < n; ++i)
                {
                    overlap += std::conj(current(other, i)) * current(band, i);
                }
                const T expected = band == other ? T(1) : T(0);
                EXPECT_LT(std::abs(overlap - expected), tolerance);
            }
        }
    }
}

template <typename Real>
class PWTDDFTTest : public testing::Test
{
};
using Precisions = testing::Types<double, float>;
TYPED_TEST_SUITE(PWTDDFTTest, Precisions);

TYPED_TEST(PWTDDFTTest, DenseCNCorrectorsAndConservation)
{
    using T = std::complex<TypeParam>;
    const bool single_precision = std::is_same<TypeParam, float>::value;
    const double energy_tolerance = single_precision ? 1e-4 : 1e-12;
    const std::vector<int> sizes{9, 6};
    const int bands = 3;
    const int ld = 12;
    const double dt = 0.27;
    const ModuleBase::Vector3<double> shift(0.13, -0.07, 0.02);
#ifdef __MPI
    const hsolver::diag_comm_info comm(MPI_COMM_SELF, 0, 1);
#else
    const hsolver::diag_comm_info comm(0, 1);
#endif
    for (const std::string method: {"bicgstab", "cgs", "gmres", "bicgstab_cn", "cgs_cn"})
    {
        for (const std::string precond: {"none", "kinetic", "kinetic_recycle", "kinetic_subspace"})
        {
            const bool recycle = precond == "kinetic_recycle";
            const bool cn_variant = method == "bicgstab_cn" || method == "cgs_cn";
            if (cn_variant && precond == "none")
            {
                continue;
            }
            if (!cn_variant && method != "gmres" && (recycle || precond == "kinetic_subspace"))
            {
                continue;
            }
            SCOPED_TRACE(method + "/" + precond);
            ModulePW::PW_Basis_K basis;
            initialize_basis(sizes, &basis);
            DenseHamiltonian<T> op(sizes);
            std::ostringstream log;
            hsolver::PWLinearOptions options;
            const std::string solver_name = cn_variant ? method.substr(0, method.size() - 3) : method;
            options.linear.method = hsolver::parse_linear_method(solver_name);
            options.linear.tolerance = single_precision ? 2e-6 : 1e-13;
            options.linear.max_iterations = 100;
            options.preconditioner = hsolver::parse_pw_precond(precond);
            options.kinetic_enabled = true;
            options.linear.restart = 2;
            options.linear.reconstruct = method == "gmres" && recycle;
            options.cn_init = method == "gmres" || cn_variant;
            hsolver::HSolverPWTDDFT<T, Device> solver(basis, options, comm, log);
            psi::Psi<T> previous(2, bands, ld, sizes, true);
            for (int ik = 0; ik < 2; ++ik)
            {
                previous.fix_k(ik);
                for (int band = 0; band < bands; ++band)
                {
                    for (int i = 0; i < ld; ++i)
                    {
                        const double phase = 2.0 * std::acos(-1.0) * i * band / sizes[ik];
                        previous(band, i) = i < sizes[ik] ? T(std::polar(1.0 / std::sqrt(sizes[ik]), phase)) : T(0);
                    }
                }
            }
            psi::Psi<T> current(previous);
            for (int step = 1; step <= 12; ++step)
            {
                // A time-step change must invalidate both k-point histories before reuse.
                const double step_dt = step < 4 ? dt : dt * 0.5;
                log.str("");
                log.clear();
                solver.solve(op, previous, &current, step_dt, shift, step, 1, true, log);
                check_dense_step(op, previous, current, step_dt);
                if (recycle || precond == "kinetic_subspace")
                {
                    std::istringstream records(log.str());
                    std::string line;
                    int count = 0;
                    while (std::getline(records, line))
                    {
                        const std::string key = "coarse_rank=";
                        const std::size_t position = line.find(key);
                        if (position == std::string::npos)
                        {
                            continue;
                        }
                        const int rank = std::stoi(line.substr(position + key.size()));
                        if (recycle && (step == 1 || step == 4))
                        {
                            EXPECT_EQ(rank, 0);
                        }
                        else
                        {
                            EXPECT_GT(rank, 0);
                        }
                        EXPECT_NE(line.find("cn_initial=1"), std::string::npos);
                        EXPECT_NE(line.find("kinetic_retry=0"), std::string::npos);
                        ++count;
                    }
                    EXPECT_EQ(count, sizes.size());
                }
                if (step == 1)
                {
                    // A corrector rebuilds the RHS from the fixed previous state, even with a changed H and guess.
                    for (int ik = 0; ik < 2; ++ik)
                    {
                        const int n = sizes[ik];
                        current.fix_k(ik);
                        for (int i = 0; i < n; ++i)
                        {
                            op.matrices[ik][i * n + i] += 0.15 * (i + 1);
                        }
                        for (int band = 0; band < bands; ++band)
                        {
                            for (int i = 0; i < n; ++i)
                            {
                                current(band, i) *= T(0.8, 0.1);
                            }
                        }
                    }
                    log.str("");
                    log.clear();
                    solver.solve(op, previous, &current, step_dt, shift, step, 2, true, log);
                    check_dense_step(op, previous, current, step_dt);
                    EXPECT_EQ(log.str().find("cn_initial=1"), std::string::npos);
                }
                previous = current;
            }
            // Evaluate expectations with a different endpoint H, not the last propagation H.
            for (int ik = 0; ik < 2; ++ik)
            {
                for (int i = 0; i < sizes[ik]; ++i)
                {
                    op.matrices[ik][i * sizes[ik] + i] += 0.3;
                }
            }
            ModuleBase::matrix energies(2, bands);
            solver.cal_band_energy(op, current, &energies);
            for (int ik = 0; ik < 2; ++ik)
            {
                current.fix_k(ik);
                const int n = sizes[ik];
                for (int band = 0; band < bands; ++band)
                {
                    T expected(0);
                    for (int j = 0; j < n; ++j)
                    {
                        for (int i = 0; i < n; ++i)
                        {
                            expected += std::conj(current(band, i)) * op.matrices[ik][j * n + i] * current(band, j);
                        }
                    }
                    EXPECT_NEAR(energies(ik, band), expected.real(), energy_tolerance);
                }
            }
        }
    }
}
TYPED_TEST(PWTDDFTTest, UnequalOccupationsAndDiagnosticsPreservePropagation)
{
    using T = std::complex<TypeParam>;
    using Wide = linear_test::Complex;
    const bool single = std::is_same<TypeParam, float>::value;
    const double tolerance = single ? 1e-4 : 1e-10;
    const int n = 6;
    const int bands = 3;
    const int ld = 8;
    const double dt = 0.27;
    const std::vector<int> sizes{n};
    const ModuleBase::Vector3<double> shift(0, 0, 0);
    ModulePW::PW_Basis_K basis;
    initialize_basis(sizes, &basis);
    DenseHamiltonian<T> op(sizes);
    const std::vector<Wide> h(op.matrices[0].begin(), op.matrices[0].end());
#ifdef __MPI
    const hsolver::diag_comm_info comm(MPI_COMM_SELF, 0, 1);
#else
    const hsolver::diag_comm_info comm(0, 1);
#endif
    ModuleBase::matrix occupations(1, bands);
    occupations(0, 0) = 1.0;
    occupations(0, 1) = 0.4;
    occupations(0, 2) = 0.0;
    std::vector<T> without_stat;
    for (const bool out_stat: {false, true})
    {
        std::ostringstream log;
        hsolver::PWLinearOptions options;
        options.linear.method = hsolver::LinearMethod::gmres;
        options.linear.tolerance = single ? 2e-6 : 1e-13;
        options.linear.max_iterations = 100;
        options.linear.restart = n;
        options.preconditioner = hsolver::PWPreconditioner::none;
        options.out_stat = out_stat;
        hsolver::HSolverPWTDDFT<T, Device> solver(basis, options, comm, log);
        psi::Psi<T> previous(1, bands, ld, sizes, true);
        std::fill_n(previous.get_pointer(), ld * bands, T(0));
        std::vector<Wide> reference(n * bands, Wide(0));
        for (int band = 0; band < bands; ++band)
        {
            previous(band, band) = T(1);
            reference[band * n + band] = Wide(1);
        }
        psi::Psi<T> initial(previous);
        initial(0, 0) = T(1.01);
        const psi::Psi<T> original(initial);
        solver.check_initial(initial, 1);
        EXPECT_EQ(std::memcmp(initial.get_pointer(), original.get_pointer(), ld * bands * sizeof(T)), 0);
        // The finite initial defect is diagnostic; propagate the orthonormal reference state.
        psi::Psi<T> current(previous);
        for (int step = 1; step <= 8; ++step)
        {
            solver.reset_orth_stats();
            solver.solve(op, previous, &current, dt, shift, step, 1, false, log);
            reference = dense_cn(h, reference, n, bands, dt);
            for (int band = 0; band < bands; ++band)
            {
                for (int i = 0; i < n; ++i)
                {
                    EXPECT_LT(std::abs(Wide(current(band, i)) - reference[band * n + i]), tolerance);
                }
            }
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i < n; ++i)
                {
                    Wide actual(0);
                    Wide expected(0);
                    for (int band = 0; band < bands; ++band)
                    {
                        actual += occupations(0, band) * Wide(current(band, i)) * std::conj(Wide(current(band, j)));
                        expected += occupations(0, band) * reference[band * n + i] * std::conj(reference[band * n + j]);
                    }
                    // Reference columns have unit norm, hence every coefficient has magnitude <= 1.
                    const double density_bound = 1.4 * (2.0 * tolerance + tolerance * tolerance);
                    EXPECT_LT(std::abs(actual - expected), density_bound);
                }
            }
            const double expected_electrons = out_stat ? 1.4 : 0.0;
            EXPECT_NEAR(solver.wave_electrons(occupations), expected_electrons, tolerance);
            const T* data = current.get_pointer();
            if (!out_stat)
            {
                without_stat.insert(without_stat.end(), data, data + ld * bands);
            }
            else
            {
                const T* baseline = without_stat.data() + (step - 1) * ld * bands;
                EXPECT_EQ(std::memcmp(data, baseline, ld * bands * sizeof(T)), 0);
            }
            previous = current;
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    int nproc = 1;
    int nthread = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, nproc, nthread, rank);
#ifdef __MPI
    // This numerical test creates no application pools; the shared finalizer owns these handles.
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
