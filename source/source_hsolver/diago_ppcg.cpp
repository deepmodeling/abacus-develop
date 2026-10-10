#include "diago_ppcg.h"
#include "source_base/parallel_reduce.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <ATen/kernels/lapack.h>
#include <stdexcept>
#include "source_base/kernels/math_kernel_op.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>

namespace hsolver {
namespace {

const int ppcg_openmp_work_threshold = 4096;
const int ppcg_openmp_column_threshold = 16;
const double ppcg_minimum_diagonalization_threshold = 1.0e-14;
const double ppcg_preconditioner_threshold = 1.0e-12;
const double ppcg_numerical_threshold = 1.0e-30;
const double ppcg_scaling_threshold = 1.0e-15;
const int ppcg_stagnation_restart_interval = 15;
// Relative cutoff for overlap eigenvalues in a local [X, W, P] subspace.
// Directions below this level are numerically dependent; retaining them can
// make the reduced generalized solve stagnate instead of improving a Ritz pair.
const double ppcg_subspace_rank_threshold = 1.0e-8;

// Diagonal shifts used by the small projected eigensolve fallback.
const double ppcg_subspace_shifts[] = {0.0, 1.0e-10, 1.0e-8, 1.0e-6};

} // namespace
} // namespace hsolver



namespace hsolver {
namespace {

template <typename Value>
void reduce_pool_if_mpi_ready(Value& value)
{
#ifdef __MPI
    if (Parallel_Reduce::mpi_ready())
    {
        Parallel_Reduce::reduce_pool(value);
    }
#endif
}

template <typename Value>
void reduce_pool_if_mpi_ready(Value* value, const int n)
{
#ifdef __MPI
    if (Parallel_Reduce::mpi_ready())
    {
        Parallel_Reduce::reduce_pool(value, n);
    }
#endif
}

bool all_pool_operations_succeeded(const bool local_success)
{
    int failure_count = local_success ? 0 : 1;
    reduce_pool_if_mpi_ready(&failure_count, 1);
    return failure_count == 0;
}

template <typename T, typename Real>
Real max_generalized_residual(
    const T* hpsi,
    const T* spsi,
    const Real* eigenvalue,
    int ld,
    int n_dim,
    int ncol)
{
    Real max_res = 0;
    std::vector<double> nrm2_all(ncol, 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim * ncol > ppcg_openmp_work_threshold)
#endif
    for (int j = 0; j < ncol; ++j)
    {
        double nrm2 = 0.0;
        for (int ig = 0; ig < n_dim; ++ig)
        {
            const T r = hpsi[ig + j * ld] - T(eigenvalue[j]) * spsi[ig + j * ld];
            nrm2 += double(std::norm(r));
        }
        nrm2_all[j] = nrm2;
    }
    reduce_pool_if_mpi_ready(nrm2_all.data(), ncol);
    for (int j = 0; j < ncol; ++j)
    {
        const Real residual_norm = std::sqrt(Real(nrm2_all[j]));
        if (!std::isfinite(residual_norm))
        {
            return std::numeric_limits<Real>::infinity();
        }
        max_res = std::max(max_res, residual_norm);
    }
    return max_res;
}

template <typename T>
inline void set_zero(std::vector<T>& x)
{
    const int n = int(x.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n > ppcg_openmp_work_threshold)
#endif
    for (int i = 0; i < n; ++i)
    {
        x[i] = T(0);
    }
}

} // anonymous namespace
} // namespace hsolver



namespace hsolver {
namespace {

template <typename Scalar>
struct HermitianLapack
{
    using Real = typename container::GetTypeReal<Scalar>::type;
    using Device = container::DEVICE_CPU;

    static void sygvd(int n, Scalar* a, Scalar* b, Real* w)
    {
        std::vector<Scalar> eigenvectors(n * n);
        container::kernels::lapack_hegvd<Scalar, Device>()(
            n, n, a, b, w, eigenvectors.data());
        std::copy(eigenvectors.begin(), eigenvectors.end(), a);
    }
};

} // anonymous namespace
} // namespace hsolver



namespace hsolver {
namespace {

inline bool ppcg_contiguous_cols(const std::vector<int>& cols, int& first)
{
    if (cols.empty())
    {
        return false;
    }

    first = cols.front();
    for (int j = 0; j < int(cols.size()); ++j)
    {
        if (cols[j] != first + j)
        {
            return false;
        }
    }
    return true;
}

} // anonymous namespace

// =============================================================================
// Constructor
// =============================================================================
template <typename T, typename Device>
DiagoPPCG<T, Device>::DiagoPPCG(const Real& diag_thr,
                                 const int& diag_iter_max,
                                 const int& sbsize,
                                 const int& rr_step,
                                 const bool gamma_g0_real)
    : maxiter_(diag_iter_max),
      sbsize_(std::max(1, sbsize)),
      rr_step_(std::max(1, rr_step)),
      diag_thr_(std::max(diag_thr, Real(ppcg_minimum_diagonalization_threshold))),
      gamma_g0_real_(gamma_g0_real)
{
}

// =============================================================================
// Input validation
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::validate_input(
    const HPsiFunc& hpsi_func,
    const T* psi_in,
    const Real* eigenvalue_in,
    const std::vector<double>& ethr_band,
    const Real* prec) const
{
    if (!hpsi_func)
    {
        throw std::invalid_argument("PPCG: H operator is empty.");
    }
    if (psi_in == nullptr || eigenvalue_in == nullptr)
    {
        throw std::invalid_argument("PPCG: psi/eigenvalue pointer is null.");
    }
    if (prec == nullptr)
    {
        throw std::invalid_argument("PPCG: preconditioner pointer is null.");
    }
    if (ld_psi_ <= 0 || n_band_ <= 0 || n_dim_ <= 0)
    {
        throw std::invalid_argument("PPCG: invalid dimensions.");
    }
    if (n_dim_ > ld_psi_)
    {
        throw std::invalid_argument("PPCG: dim must not exceed ld_psi.");
    }
    if (ethr_band.size() < size_t(n_band_))
    {
        throw std::invalid_argument("PPCG: ethr_band size is smaller than nband.");
    }
    for (int i = 0; i < n_band_; ++i)
    {
        if (!std::isfinite(ethr_band[i]))
        {
            throw std::invalid_argument("PPCG: ethr_band contains non-finite value.");
        }
    }
    for (int i = 0; i < n_dim_; ++i)
    {
        if (!std::isfinite(prec[i]))
        {
            throw std::invalid_argument("PPCG: preconditioner contains non-finite value.");
        }
    }
}

// =============================================================================
// Gamma-point symmetry: enforce real-valued first element
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::force_g0_real(T* x, int ncol) const
{
    if (!gamma_g0_real_ || n_dim_ <= 0)
    {
        return;
    }
    for (int j = 0; j < ncol; ++j)
    {
        x[idx(0, j, ld_psi_)] = T(std::real(x[idx(0, j, ld_psi_)]), 0.0);
    }
}

// =============================================================================
// Operator application
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::apply_h(const HPsiFunc& hpsi_func,
                                    T* psi_in, T* hpsi_out,
                                    int ncol) const
{
    hpsi_func(psi_in, hpsi_out, ld_psi_, ncol);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::apply_s(const SPsiFunc& spsi_func,
                                    T* psi_in, T* spsi_out,
                                    int ncol) const
{
    if (spsi_func)
    {
        spsi_func(psi_in, spsi_out, ld_psi_, ncol);
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ld_psi_ * ncol > ppcg_openmp_work_threshold)
#endif
        for (int j = 0; j < ncol; ++j)
        {
            std::copy(psi_in + j * ld_psi_, psi_in + (j + 1) * ld_psi_,
                      spsi_out + j * ld_psi_);
        }
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::apply_s_current(T* psi_in, T* spsi_out,
                                            int ncol) const
{
    apply_s(spsi_func_, psi_in, spsi_out, ncol);
}

// =============================================================================
// Inner product <x|y> (real part only, for Hermitian operators)
// =============================================================================
template <typename T, typename Device>
typename DiagoPPCG<T, Device>::Real
DiagoPPCG<T, Device>::gamma_dot(const T* x, const T* y) const
{
    Real result = ModuleBase::dot_real_op<T, Device>()(n_dim_, x, y, false);
    reduce_pool_if_mpi_ready(result);
    return result;
}

// =============================================================================
// Gram matrix: out[i, j] = <a_i | b_j>
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::gram(const T* mat_a, const T* mat_b,
                                 int ncol_a, int ncol_b,
                                 std::vector<T>& out,
                                 int ld_out) const
{
    out.resize(ld_out * ncol_b);
    const T one = T(1);
    const T zero = T(0);
    ModuleBase::gemm_op<T, Device>()('C',
                                     'N',
                                     ncol_a,
                                     ncol_b,
                                     n_dim_,
                                     &one,
                                     mat_a,
                                     ld_psi_,
                                     mat_b,
                                     ld_psi_,
                                     &zero,
                                     out.data(),
                                     ld_out);
    reduce_pool_if_mpi_ready(out.data(), ld_out * ncol_b);
}

// =============================================================================
// Column gather: extract selected columns into contiguous storage
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::copy_cols(const T* src,
                                      const std::vector<int>& cols,
                                      std::vector<T>& dst) const
{
    const int ncols = int(cols.size());
    dst.resize(ld_psi_ * ncols);
    if (ncols == 0)
    {
        return;
    }

    int first = 0;
    if (ppcg_contiguous_cols(cols, first))
    {
        std::copy(src + first * ld_psi_,
                  src + (first + ncols) * ld_psi_,
                  dst.begin());
        return;
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ld_psi_ * ncols > ppcg_openmp_work_threshold)
#endif
    for (int j = 0; j < ncols; ++j)
    {
        const int c = cols[j];
        std::copy(src + c * ld_psi_, src + c * ld_psi_ + ld_psi_,
                  dst.begin() + j * ld_psi_);
    }
}

// =============================================================================
// Column scatter: write contiguous storage back into selected columns
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::scatter_cols(
    T* dst,
    const std::vector<int>& cols,
    const std::vector<T>& src) const
{
    const int ncols = int(cols.size());
    if (ncols == 0)
    {
        return;
    }

    int first = 0;
    if (ppcg_contiguous_cols(cols, first))
    {
        std::copy(src.begin(),
                  src.begin() + ld_psi_ * ncols,
                  dst + first * ld_psi_);
        return;
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ld_psi_ * ncols > ppcg_openmp_work_threshold)
#endif
    for (int j = 0; j < ncols; ++j)
    {
        const int c = cols[j];
        std::copy(src.begin() + j * ld_psi_,
                  src.begin() + (j + 1) * ld_psi_,
                  dst + c * ld_psi_);
    }
}

// =============================================================================
// Project x onto vectors orthogonal to S-orthonormal basis
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::project_against(
    const T* basis, const T* sbasis,
    const std::vector<int>& basis_cols,
    std::vector<T>& x, std::vector<T>& sx,
    const std::vector<int>& x_cols) const
{
    if (basis_cols.empty() || x_cols.empty())
    {
        return;
    }

    const int nbasis = int(basis_cols.size());
    const int nx = int(x_cols.size());

    int x_first = 0;
    const bool contiguous_x = ppcg_contiguous_cols(x_cols, x_first);

    std::vector<T> x_l;
    std::vector<T> sx_l;
    T* x_data = x.data() + x_first * ld_psi_;
    T* sx_data = sx.data() + x_first * ld_psi_;
    if (!contiguous_x)
    {
        x_l.reserve(ld_psi_ * nx);
        sx_l.reserve(ld_psi_ * nx);
        copy_cols(x.data(), x_cols, x_l);
        copy_cols(sx.data(), x_cols, sx_l);
        x_data = x_l.data();
        sx_data = sx_l.data();
    }

    int basis_first = 0;
    const bool contiguous_basis =
        ppcg_contiguous_cols(basis_cols, basis_first);

    std::vector<T> basis_l;
    std::vector<T> sbasis_l;
    const T* basis_data = basis + basis_first * ld_psi_;
    const T* sbasis_data = sbasis + basis_first * ld_psi_;
    if (!contiguous_basis)
    {
        basis_l.reserve(ld_psi_ * nbasis);
        sbasis_l.reserve(ld_psi_ * nbasis);
        copy_cols(basis, basis_cols, basis_l);
        copy_cols(sbasis, basis_cols, sbasis_l);
        basis_data = basis_l.data();
        sbasis_data = sbasis_l.data();
    }

    std::vector<T> coeff(nbasis * nx, T(0));
    gram(basis_data, sx_data, nbasis, nx, coeff, nbasis);

    const T minus_one = T(-1);
    const T one = T(1);
    ModuleBase::gemm_op<T, Device>()('N',
                                     'N',
                                     n_dim_,
                                     nx,
                                     nbasis,
                                     &minus_one,
                                     basis_data,
                                     ld_psi_,
                                     coeff.data(),
                                     nbasis,
                                     &one,
                                     x_data,
                                     ld_psi_);
    ModuleBase::gemm_op<T, Device>()('N',
                                     'N',
                                     n_dim_,
                                     nx,
                                     nbasis,
                                     &minus_one,
                                     sbasis_data,
                                     ld_psi_,
                                     coeff.data(),
                                     nbasis,
                                     &one,
                                     sx_data,
                                     ld_psi_);

    if (!contiguous_x)
    {
        scatter_cols(x.data(), x_cols, x_l);
        scatter_cols(sx.data(), x_cols, sx_l);
    }
}

// =============================================================================
// Preconditioner: x[c] /= max(prec, eps) for each active column c
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::divide_by_preconditioner(
    const std::vector<int>& active_cols,
    const Real* prec,
    std::vector<T>& x) const
{
    const int ncols = int(active_cols.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ * ncols > ppcg_openmp_work_threshold)
#endif
    for (int j = 0; j < ncols; ++j)
    {
        const int c = active_cols[j];
        for (int ig = 0; ig < n_dim_; ++ig)
        {
            x[idx(ig, c, ld_psi_)] /=
                std::max(prec[ig], Real(ppcg_preconditioner_threshold));
        }
    }
}

} // namespace hsolver


namespace hsolver {

//==============================================================================
// BLOCK_SUBSPACE STRATEGY
//==============================================================================

// ---------------------------------------------------------------------------
// Lock only a contiguous prefix of converged eigenpairs.  A small residual
// proves that a vector is an eigenpair, but it does not prove that it belongs
// to the requested lowest roots.  Keeping every band after the first
// unconverged root active lets an evolving lower root replace a later exact
// high-energy state.
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::lock_epairs(
    const std::vector<T>& residual,
    const std::vector<double>& ethr_band,
    std::vector<int>& active_cols) const
{
    active_cols.clear();
    active_cols.reserve(n_band_);
    std::vector<double> nrm2_all(n_band_, 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ * n_band_ > ppcg_openmp_work_threshold)
#endif
    for (int j = 0; j < n_band_; ++j)
    {
        double nrm2 = 0.0;
        for (int ig = 0; ig < n_dim_; ++ig)
        {
            nrm2 += double(std::norm(residual[idx(ig, j, ld_psi_)]));
        }
        nrm2_all[j] = nrm2;
    }
    reduce_pool_if_mpi_ready(nrm2_all.data(), n_band_);
    bool found_unconverged = false;
    for (int j = 0; j < n_band_; ++j)
    {
        const Real rnrm = std::sqrt(std::max(Real(nrm2_all[j]), Real(0)));
        const Real thr = std::max(Real(ethr_band[j]), diag_thr_);
        if (!std::isfinite(rnrm) || rnrm > thr)
        {
            found_unconverged = true;
        }
        if (found_unconverged)
        {
            active_cols.push_back(j);
        }
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::scale_to_unit_snorm(std::vector<T>& x,
                                                std::vector<T>& sx,
                                                std::vector<T>& hx,
                                                const int ncols) const
{
    std::vector<double> inverse_snorm(ncols, 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ * ncols > ppcg_openmp_work_threshold)
#endif
    for (int j = 0; j < ncols; ++j)
    {
        double snorm_squared = 0.0;
        for (int ig = 0; ig < n_dim_; ++ig)
        {
            snorm_squared += double(std::real(
                std::conj(x[idx(ig, j, ld_psi_)])
                * sx[idx(ig, j, ld_psi_)]));
        }
        inverse_snorm[j] = snorm_squared;
    }
    reduce_pool_if_mpi_ready(inverse_snorm.data(), ncols);
    for (int j = 0; j < ncols; ++j)
    {
        const Real snorm = std::sqrt(std::max(
            Real(inverse_snorm[j]), Real(ppcg_numerical_threshold)));
        inverse_snorm[j] = snorm > Real(ppcg_scaling_threshold)
                             ? double(Real(1) / snorm)
                             : 1.0;
    }
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) if (n_dim_ * ncols > ppcg_openmp_work_threshold)
#endif
    for (int j = 0; j < ncols; ++j)
    {
        for (int ig = 0; ig < n_dim_; ++ig)
        {
            const Real scale = Real(inverse_snorm[j]);
            x[idx(ig, j, ld_psi_)] *= scale;
            sx[idx(ig, j, ld_psi_)] *= scale;
            hx[idx(ig, j, ld_psi_)] *= scale;
        }
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::hermitize_projected(std::vector<T>& matrix,
                                                const int dim) const
{
    for (int j = 0; j < dim; ++j)
    {
        matrix[j + j * dim] = T(std::real(matrix[j + j * dim]), 0);
        for (int i = j + 1; i < dim; ++i)
        {
            const T value = (matrix[i + j * dim]
                             + std::conj(matrix[j + i * dim]))
                            * Real(0.5);
            matrix[i + j * dim] = value;
            matrix[j + i * dim] = std::conj(value);
        }
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::insert_startup_gram_block(
    const T* left,
    const T* right,
    const int row_offset,
    const int column_offset,
    std::vector<T>& matrix,
    std::vector<T>& workspace) const
{
    const int block_size = n_band_;
    const int matrix_dim = 2 * n_band_;
    gram(left, right, block_size, block_size, workspace, block_size);
    for (int j = 0; j < block_size; ++j)
    {
        for (int i = 0; i < block_size; ++i)
        {
            const int matrix_index = row_offset + i
                                     + (column_offset + j) * matrix_dim;
            matrix[matrix_index] = workspace[i + j * block_size];
        }
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::build_startup_global_subspace(
    const T* psi,
    SmallSubspace& subspace)
{
    scale_to_unit_snorm(w_, sw_, hw_, n_band_);

    const int local_dim = 2 * n_band_;
    const int matrix_size = local_dim * local_dim;
    subspace.k.assign(matrix_size, T(0));
    subspace.m.assign(matrix_size, T(0));
    subspace.eval.resize(local_dim);
    std::vector<T> gram_workspace;

    insert_startup_gram_block(psi,
                              hpsi_.data(),
                              0,
                              0,
                              subspace.k,
                              gram_workspace);
    insert_startup_gram_block(psi,
                              hw_.data(),
                              0,
                              n_band_,
                              subspace.k,
                              gram_workspace);
    insert_startup_gram_block(w_.data(),
                              hpsi_.data(),
                              n_band_,
                              0,
                              subspace.k,
                              gram_workspace);
    insert_startup_gram_block(w_.data(),
                              hw_.data(),
                              n_band_,
                              n_band_,
                              subspace.k,
                              gram_workspace);
    insert_startup_gram_block(psi,
                              spsi_.data(),
                              0,
                              0,
                              subspace.m,
                              gram_workspace);
    insert_startup_gram_block(psi,
                              sw_.data(),
                              0,
                              n_band_,
                              subspace.m,
                              gram_workspace);
    insert_startup_gram_block(w_.data(),
                              spsi_.data(),
                              n_band_,
                              0,
                              subspace.m,
                              gram_workspace);
    insert_startup_gram_block(w_.data(),
                              sw_.data(),
                              n_band_,
                              n_band_,
                              subspace.m,
                              gram_workspace);

    hermitize_projected(subspace.k, local_dim);
    hermitize_projected(subspace.m, local_dim);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::update_startup_global_subspace(
    T* psi,
    const SmallSubspace& subspace)
{
    const int local_dim = 2 * n_band_;
    const T* state_coefficients = subspace.k.data();
    const T* correction_coefficients = state_coefficients + n_band_;
    const T one = T(1);
    const T zero = T(0);

    set_zero(rr_psi_);
    set_zero(rr_spsi_);
    set_zero(rr_hpsi_);
    set_zero(p_);
    set_zero(sp_);
    set_zero(hp_);

    auto combine_state = [&](const T* state,
                             const T* correction,
                             T* output)
    {
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         n_dim_,
                                         n_band_,
                                         n_band_,
                                         &one,
                                         state,
                                         ld_psi_,
                                         state_coefficients,
                                         local_dim,
                                         &zero,
                                         output,
                                         ld_psi_);
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         n_dim_,
                                         n_band_,
                                         n_band_,
                                         &one,
                                         correction,
                                         ld_psi_,
                                         correction_coefficients,
                                         local_dim,
                                         &one,
                                         output,
                                         ld_psi_);
    };
    auto combine_direction = [&](const T* correction, T* output)
    {
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         n_dim_,
                                         n_band_,
                                         n_band_,
                                         &one,
                                         correction,
                                         ld_psi_,
                                         correction_coefficients,
                                         local_dim,
                                         &zero,
                                         output,
                                         ld_psi_);
    };

    combine_state(psi, w_.data(), rr_psi_.data());
    combine_state(spsi_.data(), sw_.data(), rr_spsi_.data());
    combine_state(hpsi_.data(), hw_.data(), rr_hpsi_.data());
    combine_direction(w_.data(), p_.data());
    combine_direction(sw_.data(), sp_.data());
    combine_direction(hw_.data(), hp_.data());

    std::copy(rr_psi_.begin(), rr_psi_.end(), psi);
    std::copy(rr_spsi_.begin(), rr_spsi_.end(), spsi_.begin());
    std::copy(rr_hpsi_.begin(), rr_hpsi_.end(), hpsi_.begin());
}

// ---------------------------------------------------------------------------
// Build K = V^H H V and M = V^H S V where V = [psi, w]
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::build_small_subspace(
    const T* psi,
    const std::vector<int>& cols,
    int nblk,
    SmallSubspace& subspace) const
{
    const int l = int(cols.size());
    const int dim = nblk * l;
    subspace.k.resize(dim * dim);
    subspace.m.resize(dim * dim);
    subspace.eval.resize(dim);

    copy_cols(psi, cols, subspace.psi_l);
    copy_cols(spsi_.data(), cols, subspace.spsi_l);
    copy_cols(hpsi_.data(), cols, subspace.hpsi_l);
    copy_cols(w_.data(), cols, subspace.w_l);
    copy_cols(sw_.data(), cols, subspace.sw_l);
    copy_cols(hw_.data(), cols, subspace.hw_l);
    if (nblk >= 3)
    {
        copy_cols(p_.data(), cols, subspace.p_l);
        copy_cols(sp_.data(), cols, subspace.sp_l);
        copy_cols(hp_.data(), cols, subspace.hp_l);
    }

    // ---------------------------------------------------------------------------
    // Normalize w columns to unit S-norm for numerical stability.
    //
    // The w block of the Gram matrix M has entries O(||w||^2) which become
    // tiny when residuals are small, making M nearly singular and causing
    // sygvd to produce garbage eigenvectors.
    //
    // Scaling to unit S-norm keeps M well-conditioned (diagonal ~1) without
    // changing the subspace. The same scaled basis is reused in update_one_block.
    // ---------------------------------------------------------------------------
    scale_to_unit_snorm(subspace.w_l,
                        subspace.sw_l,
                        subspace.hw_l,
                        l);
    if (nblk >= 3)
    {
        scale_to_unit_snorm(subspace.p_l,
                            subspace.sp_l,
                            subspace.hp_l,
                            l);
    }

    auto copy_block = [&](const std::vector<T>& src,
                          const int col0,
                          std::vector<T>& dst)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ld_psi_ * l > ppcg_openmp_work_threshold)
#endif
        for (int j = 0; j < l; ++j)
        {
            std::copy(src.begin() + j * ld_psi_,
                      src.begin() + (j + 1) * ld_psi_,
                      dst.begin() + (col0 + j) * ld_psi_);
        }
    };

    subspace.basis.resize(ld_psi_ * dim);
    subspace.hbasis.resize(ld_psi_ * dim);
    subspace.sbasis.resize(ld_psi_ * dim);
    copy_block(subspace.psi_l, 0, subspace.basis);
    copy_block(subspace.hpsi_l, 0, subspace.hbasis);
    copy_block(subspace.spsi_l, 0, subspace.sbasis);
    copy_block(subspace.w_l, l, subspace.basis);
    copy_block(subspace.hw_l, l, subspace.hbasis);
    copy_block(subspace.sw_l, l, subspace.sbasis);
    if (nblk >= 3)
    {
        copy_block(subspace.p_l, 2 * l, subspace.basis);
        copy_block(subspace.hp_l, 2 * l, subspace.hbasis);
        copy_block(subspace.sp_l, 2 * l, subspace.sbasis);
    }

    gram(subspace.basis.data(), subspace.hbasis.data(), dim, dim, subspace.k, dim);
    gram(subspace.basis.data(), subspace.sbasis.data(), dim, dim, subspace.m, dim);
    hermitize_projected(subspace.k, dim);
    hermitize_projected(subspace.m, dim);
}

// ---------------------------------------------------------------------------
// Solve K v = λ M v (small generalized eigenvalue problem)
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::solve_small_generalized(
    const int dim,
    const int nstates,
    SmallSubspace& subspace) const
{
    const std::vector<T> k0 = subspace.k;
    const std::vector<T> m0 = subspace.m;

    // A nearly dependent local basis makes the generalized solve unstable
    // even when LAPACK accepts M as positive definite.  Work in the retained
    // eigenspace of M and map the Ritz vectors back to the original basis.
    std::vector<T> overlap_eigenvectors = m0;
    std::vector<Real> overlap_eigenvalues(dim);
    bool local_rank_truncation_succeeded = false;
    try
    {
        container::kernels::lapack_heevd<T, container::DEVICE_CPU>()(
            dim,
            overlap_eigenvectors.data(),
            dim,
            overlap_eigenvalues.data());

        const Real largest_overlap_eigenvalue = std::max(overlap_eigenvalues.back(), Real(0));
        const Real precision_floor = Real(10) * std::numeric_limits<Real>::epsilon();
        const Real relative_rank_threshold = std::max(
            Real(ppcg_subspace_rank_threshold),
            precision_floor);
        const Real overlap_cutoff = largest_overlap_eigenvalue
                                    * relative_rank_threshold;
        std::vector<int> retained_indices;
        retained_indices.reserve(dim);
        for (int i = 0; i < dim; ++i)
        {
            if (overlap_eigenvalues[i] > overlap_cutoff)
            {
                retained_indices.push_back(i);
            }
        }

        const int retained_dim = int(retained_indices.size());
        if (retained_dim >= nstates && retained_dim < dim)
        {
            std::vector<T> retained_basis(dim * retained_dim, T(0));
            for (int j = 0; j < retained_dim; ++j)
            {
                const int source_index = retained_indices[j];
                const Real scale = Real(1)
                                   / std::sqrt(overlap_eigenvalues[source_index]);
                for (int row = 0; row < dim; ++row)
                {
                    retained_basis[row + j * dim] = overlap_eigenvectors[row + source_index * dim]
                                                      * scale;
                }
            }

            // Q = U_r Lambda_r^(-1/2) orthonormalizes the retained overlap
            // subspace.  Use BLAS for Q^H K Q rather than a four-index
            // contraction: the latter becomes prohibitively expensive when
            // the first global correction contains many nearly dependent
            // residual directions.
            const T one = T(1);
            const T zero = T(0);
            std::vector<T> projected_hamiltonian(dim * retained_dim, T(0));
            ModuleBase::gemm_op<T, Device>()('N',
                                             'N',
                                             dim,
                                             retained_dim,
                                             dim,
                                             &one,
                                             k0.data(),
                                             dim,
                                             retained_basis.data(),
                                             dim,
                                             &zero,
                                             projected_hamiltonian.data(),
                                             dim);
            std::vector<T> reduced_hamiltonian(retained_dim * retained_dim, T(0));
            ModuleBase::gemm_op<T, Device>()('C',
                                             'N',
                                             retained_dim,
                                             retained_dim,
                                             dim,
                                             &one,
                                             retained_basis.data(),
                                             dim,
                                             projected_hamiltonian.data(),
                                             dim,
                                             &zero,
                                             reduced_hamiltonian.data(),
                                             retained_dim);
            for (int j = 0; j < retained_dim; ++j)
            {
                const Real diagonal_value = std::real(
                    reduced_hamiltonian[j + j * retained_dim]);
                reduced_hamiltonian[j + j * retained_dim] = T(diagonal_value, Real(0));
                for (int i = 0; i < j; ++i)
                {
                    const T value = Real(0.5) * (
                        reduced_hamiltonian[i + j * retained_dim]
                        + std::conj(reduced_hamiltonian[j + i * retained_dim]));
                    reduced_hamiltonian[i + j * retained_dim] = value;
                    reduced_hamiltonian[j + i * retained_dim] = std::conj(value);
                }
            }

            std::vector<Real> reduced_eigenvalues(retained_dim);
            container::kernels::lapack_heevd<T, container::DEVICE_CPU>()(
                retained_dim,
                reduced_hamiltonian.data(),
                retained_dim,
                reduced_eigenvalues.data());

            std::fill(subspace.k.begin(), subspace.k.end(), T(0));
            for (int state = 0; state < nstates; ++state)
            {
                subspace.eval[state] = reduced_eigenvalues[state];
            }
            ModuleBase::gemm_op<T, Device>()('N',
                                             'N',
                                             dim,
                                             nstates,
                                             retained_dim,
                                             &one,
                                             retained_basis.data(),
                                             dim,
                                             reduced_hamiltonian.data(),
                                             retained_dim,
                                             &zero,
                                             subspace.k.data(),
                                             dim);
            local_rank_truncation_succeeded = true;
        }
    }
    catch (const std::runtime_error&)
    {
        // The shifted generalized solve below remains the fallback.
    }
    if (all_pool_operations_succeeded(local_rank_truncation_succeeded))
    {
        return;
    }

    // Try with increasing diagonal shifts; fall back to identity (no update)
    // if the subspace is too ill-conditioned.  sygvd modifies both matrices
    // in-place before it may fail.
    const Real shifts[] = {Real(ppcg_subspace_shifts[0]),
                           Real(ppcg_subspace_shifts[1]),
                           Real(ppcg_subspace_shifts[2]),
                           Real(ppcg_subspace_shifts[3])};
    for (const Real shift : shifts)
    {
        subspace.k = k0;
        subspace.m = m0;
        for (int i = 0; i < dim; ++i)
        {
            subspace.m[i + i * dim] += T(shift);
        }

        bool local_sygvd_succeeded = false;
        try
        {
            HermitianLapack<T>::sygvd(dim, subspace.k.data(),
                                      subspace.m.data(),
                                      subspace.eval.data());
            local_sygvd_succeeded = true;
        }
        catch (const std::runtime_error&)
        {
            // Try the next diagonal shift.
        }
        if (all_pool_operations_succeeded(local_sygvd_succeeded))
        {
            return;
        }
    }
    // All attempts failed — set eigenvectors to identity (no update).
    std::fill(subspace.k.begin(), subspace.k.end(), T(0));
    for (int i = 0; i < dim; ++i)
    {
        subspace.k[i + i * dim] = T(1);
        subspace.eval[i] = Real(std::real(k0[i + i * dim]))
                         / std::max(Real(std::real(m0[i + i * dim])),
                                    Real(ppcg_numerical_threshold));
    }
}

// ---------------------------------------------------------------------------
// Update wavefunctions from small subspace eigenvectors
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::update_one_block(
    T* psi,
    const std::vector<int>& cols,
    int l,
    int nblk,
    SmallSubspace& subspace)
{
    const int dim = nblk * l;
    const T* eigvec = subspace.k.data();

    subspace.psi_new.assign(ld_psi_ * l, T(0));
    subspace.spsi_new.assign(ld_psi_ * l, T(0));
    subspace.hpsi_new.assign(ld_psi_ * l, T(0));
    subspace.p_new.assign(ld_psi_ * l, T(0));
    subspace.sp_new.assign(ld_psi_ * l, T(0));
    subspace.hp_new.assign(ld_psi_ * l, T(0));

    // coeff_state: full Ritz-vector rows [psi, w, p] -> new iterate.
    // coeff_p:     the [w, p] rows (psi rows zero) -> new search direction.
    subspace.coeff_state.resize(dim * l);
    subspace.coeff_p.resize(dim * l);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (l * l > ppcg_openmp_work_threshold)
#endif
    for (int j = 0; j < l; ++j)
    {
        for (int i = 0; i < l; ++i)
        {
            const T c_psi = eigvec[i + j * dim];
            const T c_w = eigvec[(l + i) + j * dim];
            subspace.coeff_state[i + j * dim] = c_psi;
            subspace.coeff_state[(l + i) + j * dim] = c_w;
            subspace.coeff_p[i + j * dim] = T(0);
            subspace.coeff_p[(l + i) + j * dim] = c_w;
            if (nblk >= 3)
            {
                const T c_p = eigvec[(2 * l + i) + j * dim];
                subspace.coeff_state[(2 * l + i) + j * dim] = c_p;
                subspace.coeff_p[(2 * l + i) + j * dim] = c_p;
            }
        }
    }

    auto fill_basis = [&](const std::vector<T>& a,
                          const std::vector<T>& b,
                          const std::vector<T>& c,
                          std::vector<T>& basis)
    {
        basis.resize(ld_psi_ * dim);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ld_psi_ * l > ppcg_openmp_work_threshold)
#endif
        for (int j = 0; j < l; ++j)
        {
            std::copy(a.begin() + j * ld_psi_,
                      a.begin() + (j + 1) * ld_psi_,
                      basis.begin() + j * ld_psi_);
            std::copy(b.begin() + j * ld_psi_,
                      b.begin() + (j + 1) * ld_psi_,
                      basis.begin() + (l + j) * ld_psi_);
            if (nblk >= 3)
            {
                std::copy(c.begin() + j * ld_psi_,
                          c.begin() + (j + 1) * ld_psi_,
                          basis.begin() + (2 * l + j) * ld_psi_);
            }
        }
    };

    auto combine = [&](const std::vector<T>& basis,
                       const std::vector<T>& coeff,
                       std::vector<T>& out)
    {
        const T one = T(1);
        const T zero = T(0);
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         n_dim_,
                                         l,
                                         dim,
                                         &one,
                                         basis.data(),
                                         ld_psi_,
                                         coeff.data(),
                                         dim,
                                         &zero,
                                         out.data(),
                                         ld_psi_);
    };

    fill_basis(subspace.psi_l, subspace.w_l, subspace.p_l, subspace.basis);
    fill_basis(subspace.spsi_l, subspace.sw_l, subspace.sp_l, subspace.sbasis);
    fill_basis(subspace.hpsi_l, subspace.hw_l, subspace.hp_l, subspace.hbasis);

    combine(subspace.basis, subspace.coeff_state, subspace.psi_new);
    combine(subspace.sbasis, subspace.coeff_state, subspace.spsi_new);
    combine(subspace.hbasis, subspace.coeff_state, subspace.hpsi_new);
    combine(subspace.basis, subspace.coeff_p, subspace.p_new);
    combine(subspace.sbasis, subspace.coeff_p, subspace.sp_new);
    combine(subspace.hbasis, subspace.coeff_p, subspace.hp_new);

    scatter_cols(psi, cols, subspace.psi_new);
    scatter_cols(spsi_.data(), cols, subspace.spsi_new);
    scatter_cols(hpsi_.data(), cols, subspace.hpsi_new);
    scatter_cols(p_.data(), cols, subspace.p_new);
    scatter_cols(sp_.data(), cols, subspace.sp_new);
    scatter_cols(hp_.data(), cols, subspace.hp_new);
}

// ---------------------------------------------------------------------------
// Rayleigh-Ritz: full subspace diagonalization + residual computation
// ---------------------------------------------------------------------------
template <typename T, typename Device>
bool DiagoPPCG<T, Device>::rayleigh_ritz(
    T* psi, Real* eigenvalue,
    std::vector<int>& active_cols,
    const std::vector<double>& ethr_band)
{
    gram(psi, hpsi_.data(), n_band_, n_band_, rr_hsub_, n_band_);
    gram(psi, spsi_.data(), n_band_, n_band_, rr_ssub_, n_band_);

    bool local_sygvd_ok = false;
    try
    {
        HermitianLapack<T>::sygvd(n_band_, rr_hsub_.data(), rr_ssub_.data(),
                                  rr_eval_.data());
        local_sygvd_ok = true;
    }
    catch (const std::runtime_error&)
    {
        // The collective result below selects one common fallback on all ranks.
    }
    const bool sygvd_ok = all_pool_operations_succeeded(local_sygvd_ok);
    if (!sygvd_ok)
    {
        // LAPACK may overwrite the projected matrices before failing.  Re-form
        // them on every rank so all ranks enter the same fallback path.
        gram(psi, hpsi_.data(), n_band_, n_band_, rr_hsub_, n_band_);
        gram(psi, spsi_.data(), n_band_, n_band_, rr_ssub_, n_band_);
        for (int ii = 0; ii < n_band_; ++ii)
        {
            rr_eval_[ii] = Real(std::real(rr_hsub_[ii + ii * n_band_]))
                     / std::max(Real(
                                    std::real(rr_ssub_[ii + ii * n_band_])),
                                Real(ppcg_numerical_threshold));
        }
    }

    if (sygvd_ok)
    {
        const int sz = ld_psi_ * n_band_;
        std::copy(psi, psi + sz, rr_psi_.begin());
        std::copy(spsi_.begin(), spsi_.end(), rr_spsi_.begin());
        std::copy(hpsi_.begin(), hpsi_.end(), rr_hpsi_.begin());

        std::fill(psi, psi + ld_psi_ * n_band_, T(0));
        set_zero(spsi_);
        set_zero(hpsi_);

        const T one = T(1);
        const T zero = T(0);
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         n_dim_,
                                         n_band_,
                                         n_band_,
                                         &one,
                                         rr_psi_.data(),
                                         ld_psi_,
                                         rr_hsub_.data(),
                                         n_band_,
                                         &zero,
                                         psi,
                                         ld_psi_);
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         n_dim_,
                                         n_band_,
                                         n_band_,
                                         &one,
                                         rr_spsi_.data(),
                                         ld_psi_,
                                         rr_hsub_.data(),
                                         n_band_,
                                         &zero,
                                         spsi_.data(),
                                         ld_psi_);
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         n_dim_,
                                         n_band_,
                                         n_band_,
                                         &one,
                                         rr_hpsi_.data(),
                                         ld_psi_,
                                         rr_hsub_.data(),
                                         n_band_,
                                         &zero,
                                         hpsi_.data(),
                                         ld_psi_);

        for (int j = 0; j < n_band_; ++j)
        {
            eigenvalue[j] = rr_eval_[j];
        }
    }
    else
    {
        // No rotation: just update eigenvalues with Rayleigh quotients.
        for (int j = 0; j < n_band_; ++j)
        {
            eigenvalue[j] = rr_eval_[j];
        }
    }

    // Compute residual: w_i = H|psi_i> - eps_i * S|psi_i>
    set_zero(w_);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) if (n_dim_ * n_band_ > ppcg_openmp_work_threshold)
#endif
    for (int j = 0; j < n_band_; ++j)
    {
        for (int ig = 0; ig < n_dim_; ++ig)
        {
            w_[idx(ig, j, ld_psi_)] = hpsi_[idx(ig, j, ld_psi_)]
                                    - spsi_[idx(ig, j, ld_psi_)] * eigenvalue[j];
        }
    }

    if (sygvd_ok)
    {
        lock_epairs(w_, ethr_band, active_cols);
    }
    else
    {
        // Without a successful global Rayleigh-Ritz solve, the columns are
        // not known to be ordered by eigenvalue.  Do not hard-lock an
        // arbitrary column based on its residual alone.
        active_cols.resize(n_band_);
        std::iota(active_cols.begin(), active_cols.end(), 0);
    }
    return sygvd_ok;
}

} // namespace hsolver


namespace hsolver {

//==============================================================================
// MAIN DIAGONALIZATION ROUTINE
//==============================================================================
template <typename T, typename Device>
double DiagoPPCG<T, Device>::diag(const HPsiFunc& hpsi_func,
                                  const SPsiFunc& spsi_func,
                                  int ld_psi,
                                  int nband,
                                  int dim,
                                  T* psi_in,
                                  Real* eigenvalue_in,
                                  const std::vector<double>& ethr_band,
                                  const Real* prec)
{
    ld_psi_ = ld_psi;
    n_band_ = nband;
    n_dim_ = dim;
    converged_ = false;
    active_band_count_ = 0;
    iteration_limit_reached_ = false;

    validate_input(hpsi_func, psi_in, eigenvalue_in, ethr_band, prec);
    spsi_func_ = spsi_func;

    // Allocate working storage.
    const int ncol = n_band_;
    const int sz = ld_psi_ * ncol;

    hpsi_.assign(sz, T(0));
    spsi_.assign(sz, T(0));
    w_.assign(sz, T(0));
    sw_.assign(sz, T(0));
    hw_.assign(sz, T(0));
    p_.assign(sz, T(0));
    sp_.assign(sz, T(0));
    hp_.assign(sz, T(0));
    rr_psi_.resize(sz);
    rr_spsi_.resize(sz);
    rr_hpsi_.resize(sz);
    rr_hsub_.resize(ncol * ncol);
    rr_ssub_.resize(ncol * ncol);
    rr_eval_.resize(ncol);

    std::vector<int> all_cols(ncol);
    std::iota(all_cols.begin(), all_cols.end(), 0);

    force_g0_real(psi_in, ncol);
    apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
    apply_s_current(psi_in, spsi_.data(), ncol);

    double avg_iter = 1.0;
    int iter = 1;
    std::vector<int> active_cols;
    active_cols.reserve(ncol);

    std::ofstream residual_trace;
    if (const char* path = std::getenv("ABACUS_PPCG_RESIDUAL_TRACE"))
    {
        // Optional debug trace for plotting PPCG convergence curves.
        residual_trace.open(path);
        if (residual_trace)
        {
            residual_trace << std::setprecision(std::numeric_limits<Real>::max_digits10);
            residual_trace << "iteration,stage,max_residual";
            // Test-side tools compare these Ritz values with analytical or
            // LAPACK reference eigenvalues without coupling them to the solver.
            for (int ib = 0; ib < ncol; ++ib)
            {
                residual_trace << ",eigenvalue_" << ib;
            }
            residual_trace << ",local_overlap_min,local_overlap_max,local_overlap_condition\n";
        }
    }
    Real local_overlap_min = std::numeric_limits<Real>::quiet_NaN();
    Real local_overlap_max = std::numeric_limits<Real>::quiet_NaN();
    Real local_overlap_condition = std::numeric_limits<Real>::quiet_NaN();
    auto record_residual = [&](int iteration, const char* stage) {
        if (!residual_trace)
        {
            return;
        }
        const Real max_residual = max_generalized_residual(hpsi_.data(),
                                                            spsi_.data(),
                                                            eigenvalue_in,
                                                            ld_psi_,
                                                            n_dim_,
                                                            ncol);
        residual_trace << iteration << ',' << stage << ',' << max_residual;
        for (int ib = 0; ib < ncol; ++ib)
        {
            residual_trace << ',' << eigenvalue_in[ib];
        }
        residual_trace << ',' << local_overlap_min
                       << ',' << local_overlap_max
                       << ',' << local_overlap_condition
                       << '\n';
    };

    auto measure_local_overlap = [&](const SmallSubspace& local_subspace,
                                     const int local_dim) {
        if (!residual_trace)
        {
            return;
        }
        std::vector<T> overlap = local_subspace.m;
        std::vector<Real> overlap_eigenvalues(local_dim);
        container::kernels::lapack_heevd<T, container::DEVICE_CPU>()(
            local_dim,
            overlap.data(),
            local_dim,
            overlap_eigenvalues.data());
        const Real min_eigenvalue = overlap_eigenvalues.front();
        const Real max_eigenvalue = overlap_eigenvalues.back();
        const Real denominator = std::max(std::abs(min_eigenvalue),
                                          Real(ppcg_numerical_threshold));
        const Real condition = std::abs(max_eigenvalue) / denominator;
        if (!std::isfinite(local_overlap_condition) || condition > local_overlap_condition)
        {
            local_overlap_min = min_eigenvalue;
            local_overlap_max = max_eigenvalue;
            local_overlap_condition = condition;
        }
    };

    // Initialize with Rayleigh-Ritz.
    bool last_rr_succeeded = rayleigh_ritz(
        psi_in, eigenvalue_in, active_cols, ethr_band);
    const bool needs_startup_global_correction = !active_cols.empty()
                                                 && sbsize_ < ncol;
    if (needs_startup_global_correction)
    {
        // A local block can discard a lower-root residual direction before
        // the global Rayleigh-Ritz step can compare it with other blocks.
        // Keep every residual direction during the first correction, so the
        // requested lowest roots are selected in one common subspace.
        active_cols = all_cols;
    }
    // Recompute to keep hpsi/spi consistent with rotated psi.
    apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
    apply_s_current(psi_in, spsi_.data(), ncol);
    record_residual(0, "initial_rr");

    std::vector<T> w_active;
    std::vector<T> sw_active;
    std::vector<T> hw_active;
    std::vector<T> p_active;
    std::vector<T> sp_active;
    std::vector<T> hp_active;
    w_active.reserve(sz);
    sw_active.reserve(sz);
    hw_active.reserve(sz);
    p_active.reserve(sz);
    sp_active.reserve(sz);
    hp_active.reserve(sz);
    std::vector<int> cols;
    cols.reserve(std::min(sbsize_, ncol));
    SmallSubspace subspace;
    bool use_p = false;  // previous search direction becomes available after
                         // the first block update.
    Real best_res = std::numeric_limits<Real>::max();  // best residual seen
    int no_improve = 0;  // iterations since the residual last decreased

    while (!active_cols.empty() && iter <= maxiter_)
    {
        local_overlap_min = std::numeric_limits<Real>::quiet_NaN();
        local_overlap_max = std::numeric_limits<Real>::quiet_NaN();
        local_overlap_condition = std::numeric_limits<Real>::quiet_NaN();
        const int nact = int(active_cols.size());
        const int nsb = std::max(1, (nact + sbsize_ - 1) / sbsize_);

        // Precondition the residual.
        divide_by_preconditioner(active_cols, prec, w_);
        copy_cols(w_.data(), active_cols, w_active);
        sw_active.assign(ld_psi_ * nact, T(0));
        apply_s_current(w_active.data(), sw_active.data(), nact);
        scatter_cols(sw_.data(), active_cols, sw_active);
        project_against(psi_in, spsi_.data(), all_cols, w_, sw_, active_cols);

        // Apply H to the search direction.
        copy_cols(w_.data(), active_cols, w_active);
        force_g0_real(w_active.data(), nact);
        hw_active.assign(ld_psi_ * nact, T(0));
        sw_active.assign(ld_psi_ * nact, T(0));
        scatter_cols(w_.data(), active_cols, w_active);
        apply_h(hpsi_func, w_active.data(), hw_active.data(), nact);
        apply_s_current(w_active.data(), sw_active.data(), nact);
        scatter_cols(hw_.data(), active_cols, hw_active);
        scatter_cols(sw_.data(), active_cols, sw_active);

        // S-orthogonalize the previous search direction p against psi, then
        // re-apply H/S.  The full Rayleigh-Ritz rotation re-mixes psi columns
        // every step, which would otherwise let p drift into psi's span and
        // destabilize the [psi, w, p] block subspace.
        if (use_p)
        {
            copy_cols(p_.data(), active_cols, p_active);
            sp_active.assign(ld_psi_ * nact, T(0));
            apply_s_current(p_active.data(), sp_active.data(), nact);
            scatter_cols(sp_.data(), active_cols, sp_active);
            project_against(psi_in, spsi_.data(), all_cols, p_, sp_, active_cols);

            copy_cols(p_.data(), active_cols, p_active);
            force_g0_real(p_active.data(), nact);
            hp_active.assign(ld_psi_ * nact, T(0));
            sp_active.assign(ld_psi_ * nact, T(0));
            scatter_cols(p_.data(), active_cols, p_active);
            apply_h(hpsi_func, p_active.data(), hp_active.data(), nact);
            apply_s_current(p_active.data(), sp_active.data(), nact);
            scatter_cols(hp_.data(), active_cols, hp_active);
            scatter_cols(sp_.data(), active_cols, sp_active);
        }

        avg_iter += double(nact) / double(ncol);

        // LOBPCG-style block subspace.  On the first sweep only [psi, w] is
        // available; afterwards the previous search direction p is added as a
        // third block, which restores the conjugate-gradient acceleration.
        // The w/p blocks are normalized to unit S-norm before building the
        // Gram matrix (see build_small_subspace), keeping M well-conditioned.

        // The startup correction needs a global [X, W] comparison to retain
        // the requested lowest roots.  Build its projected matrices directly
        // from the persistent blocks to avoid materializing three additional
        // ld_psi-by-nband bases.  Other sweeps retain the generic local path.
        if (iter == 1 && needs_startup_global_correction)
        {
            const int local_dim = 2 * ncol;
            build_startup_global_subspace(psi_in, subspace);
            measure_local_overlap(subspace, local_dim);
            solve_small_generalized(local_dim, ncol, subspace);
            update_startup_global_subspace(psi_in, subspace);
        }
        else
        {
            const int nblk = use_p ? 3 : 2;
            for (int isb = 0; isb < nsb; ++isb)
            {
                const int i0 = isb * sbsize_;
                const int l = std::min(sbsize_, nact - i0);
                cols.assign(active_cols.begin() + i0,
                            active_cols.begin() + i0 + l);

                build_small_subspace(psi_in, cols, nblk, subspace);
                const int local_dim = nblk * l;
                measure_local_overlap(subspace, local_dim);
                solve_small_generalized(local_dim, l, subspace);
                update_one_block(psi_in, cols, l, nblk, subspace);
            }
        }
        use_p = true;

        // Rayleigh-Ritz after each block update keeps the global subspace
        // synchronized with the updated active vectors.  The block update
        // can otherwise drift into an ill-conditioned basis before the next
        // Ritz rotation.
        last_rr_succeeded = rayleigh_ritz(
            psi_in, eigenvalue_in, active_cols, ethr_band);
        // Restart the search direction if the residual keeps rising.  With a
        // poor preconditioner the LOBPCG recurrence can stagnate (or slowly
        // diverge) instead of reducing the residual.  Requiring several
        // consecutive rises avoids resetting on a transient bump, and a reset
        // falls back to a steepest-descent step to recover the low eigenpairs.
        {
            const Real cur_res = max_generalized_residual(hpsi_.data(),
                                                          spsi_.data(),
                                                          eigenvalue_in,
                                                          ld_psi_,
                                                          n_dim_,
                                                          ncol);
            if (cur_res < best_res)
            {
                best_res = cur_res;
                no_improve = 0;
            }
            else
            {
                ++no_improve;
            }
            if (no_improve >= ppcg_stagnation_restart_interval)
            {
                std::fill(p_.begin(), p_.end(), T(0));
                std::fill(sp_.begin(), sp_.end(), T(0));
                std::fill(hp_.begin(), hp_.end(), T(0));
                use_p = false;
                no_improve = 0;
                best_res = cur_res;
            }
        }
        // The Rayleigh-Ritz rotation already keeps hpsi_/spsi_ consistent
        // with the rotated psi up to rounding; re-applying H/S exactly is
        // only needed every rr_step_ iterations to reset the accumulated
        // rounding drift.
        if ((iter % rr_step_) == 0)
        {
            apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
            apply_s_current(psi_in, spsi_.data(), ncol);
            apply_h(hpsi_func, p_.data(), hp_.data(), ncol);
            apply_s_current(p_.data(), sp_.data(), ncol);
        }
        record_residual(iter, "rayleigh_ritz");

        ++iter;
    }

    // Final consistency: ensure H|psi> and S|psi> match the returned vectors,
    // then classify convergence from those final residuals rather than from a
    // cached Rayleigh-Ritz state.
    apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
    apply_s_current(psi_in, spsi_.data(), ncol);
    set_zero(w_);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) if (n_dim_ * n_band_ > ppcg_openmp_work_threshold)
#endif
    for (int j = 0; j < n_band_; ++j)
    {
        for (int ig = 0; ig < n_dim_; ++ig)
        {
            w_[idx(ig, j, ld_psi_)] = hpsi_[idx(ig, j, ld_psi_)]
                                      - spsi_[idx(ig, j, ld_psi_)] * eigenvalue_in[j];
        }
    }
    if (last_rr_succeeded)
    {
        lock_epairs(w_, ethr_band, active_cols);
    }
    else
    {
        active_cols = all_cols;
    }
    active_band_count_ = int(active_cols.size());
    converged_ = last_rr_succeeded && active_cols.empty();
    iteration_limit_reached_ = !converged_ && iter > maxiter_;
    if (!active_cols.empty())
    {
        const char* failure_stage = iteration_limit_reached_
                                    ? "max_iterations"
                                    : "final_residual_failure";
        record_residual(iter - 1, failure_stage);
    }
    record_residual(iter - 1, "final");
    return avg_iter;
}

} // namespace hsolver

namespace hsolver {

template class DiagoPPCG<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_CPU>;

} // namespace hsolver
