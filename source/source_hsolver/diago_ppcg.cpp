#include "diago_ppcg.h"
#include "source_base/parallel_reduce.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <ATen/kernels/lapack.h>
#include <stdexcept>
#include "source_base/kernels/math_kernel_op.h"
#include <limits>
#include <cstdlib>
#include <fstream>
#include <numeric>

namespace hsolver {
namespace {

const int ppcg_openmp_work_threshold = 4096;
const int ppcg_openmp_column_threshold = 16;
const double ppcg_minimum_diagonalization_threshold = 1.0e-14;
const double ppcg_preconditioner_threshold = 1.0e-12;
const double ppcg_numerical_threshold = 1.0e-30;
const double ppcg_scaling_threshold = 1.0e-15;

// Increasing diagonal shifts used to regularize an ill-conditioned Gram matrix
// when a Cholesky factorization or a small projected generalized eigenproblem
// fails numerically. The ladder is tried from no shift up to a unit shift.
const double ppcg_cholesky_shifts[] = {0.0,   1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6,
                                       1.0e-4, 1.0e-3,  1.0e-2,  1.0e-1, 1.0};
// Subset of the shift ladder used by the small projected eigensolve fallback.
const double ppcg_subspace_shifts[] = {0.0, 1.0e-10, 1.0e-8, 1.0e-6};

// Orthogonality check tolerance expressed as a multiple of machine epsilon.
const double ppcg_orthogonality_tolerance_factor = 10.0;
// Line-search root-selection tolerance expressed as a multiple of machine epsilon.
const double ppcg_line_search_tolerance_factor = 100.0;
// Quadratic-formula coefficients in the line-search root solve (b^2 - 4ac and 2a).
const double ppcg_quadratic_discriminant_coefficient = 4.0;
const double ppcg_quadratic_root_denominator_coefficient = 2.0;

} // namespace
} // namespace hsolver



namespace hsolver {
namespace {

template <typename Value>
void reduce_pool_if_mpi_ready(Value& value)
{
#ifdef __MPI
    int initialized = 0;
    int finalized = 0;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    if (initialized && !finalized)
    {
        Parallel_Reduce::reduce_pool(value);
    }
#endif
}

template <typename Value>
void reduce_pool_if_mpi_ready(Value* value, const int n)
{
#ifdef __MPI
    int initialized = 0;
    int finalized = 0;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    if (initialized && !finalized)
    {
        Parallel_Reduce::reduce_pool(value, n);
    }
#endif
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
        max_res = std::max(max_res, std::sqrt(Real(nrm2_all[j])));
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

    static void potrf(int n, Scalar* a)
    {
        Real diag_max = 0;
        for (int i = 0; i < n; ++i)
        {
            diag_max = std::max(diag_max, std::abs(a[i + i * n]));
        }
        std::vector<Scalar> a0(a, a + n * n);

        for (const double shift : ppcg_cholesky_shifts)
        {
            std::copy(a0.begin(), a0.end(), a);
            if (shift > 0.0)
            {
                for (int i = 0; i < n; ++i)
                {
                    a[i + i * n] += Scalar(Real(shift) * std::max(diag_max, Real(1.0)), 0.0);
                }
            }
            try
            {
                container::kernels::lapack_potrf<Scalar, Device>()('U', n, a, n);
                return;
            }
            catch (const std::runtime_error&)
            {
                // Try the next diagonal shift.
            }
        }
        throw std::runtime_error("PPCG: potrf failed.");
    }

    static void trtri(int n, Scalar* a)
    {
        container::kernels::lapack_trtri<Scalar, Device>()('U', 'N', n, a, n);
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

template <typename T, typename Device>
T DiagoPPCG<T, Device>::complex_dot(const T* x, const T* y) const
{
    T acc = T(0);
    for (int i = 0; i < n_dim_; ++i)
    {
        acc += std::conj(x[i]) * y[i];
    }
    reduce_pool_if_mpi_ready(&acc, 1);
    return acc;
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
// Lock converged eigenpairs: bands whose eigenvalue stops changing between
// successive Rayleigh-Ritz steps are considered converged.  This matches the
// convergence criterion used by CG and Davidson (eigenvalue change < ethr).
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::lock_epairs(
    const Real* eigenvalue_prev,
    const Real* eigenvalue,
    const std::vector<double>& ethr_band,
    std::vector<int>& active_cols) const
{
    active_cols.clear();
    active_cols.reserve(n_band_);
    for (int j = 0; j < n_band_; ++j)
    {
        const Real thr = std::max(Real(ethr_band[j]), diag_thr_);
        const Real delta = std::abs(eigenvalue[j] - eigenvalue_prev[j]);
        if (delta > thr)
        {
            active_cols.push_back(j);
        }
    }
}

// ---------------------------------------------------------------------------
// Build K = V^H H V and M = V^H S V where V = [psi, w]
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::build_small_subspace(
    const T* psi,
    const std::vector<int>& cols,
    SmallSubspace& subspace) const
{
    const int l = int(cols.size());
    const int dim = 2 * l;
    subspace.k.resize(dim * dim);
    subspace.m.resize(dim * dim);
    subspace.eval.resize(dim);

    copy_cols(psi, cols, subspace.psi_l);
    copy_cols(spsi_.data(), cols, subspace.spsi_l);
    copy_cols(hpsi_.data(), cols, subspace.hpsi_l);
    copy_cols(w_.data(), cols, subspace.w_l);
    copy_cols(sw_.data(), cols, subspace.sw_l);
    copy_cols(hw_.data(), cols, subspace.hw_l);

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
    auto scale_to_unit_snorm = [this](std::vector<T>& x,
                                      std::vector<T>& sx,
                                      std::vector<T>& hx,
                                      int lcols) {
        std::vector<double> sn_scale_all(lcols, 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ * lcols > ppcg_openmp_work_threshold)
#endif
        for (int j = 0; j < lcols; ++j) {
            double sn2 = 0.0;
            for (int ig = 0; ig < n_dim_; ++ig)
            {
                sn2 += double(std::real(std::conj(x[idx(ig, j, ld_psi_)])
                                                     * sx[idx(ig, j, ld_psi_)]));
            }
            sn_scale_all[j] = sn2;
        }
        reduce_pool_if_mpi_ready(sn_scale_all.data(), lcols);
        for (int j = 0; j < lcols; ++j) {
            Real sn = std::sqrt(std::max(Real(sn_scale_all[j]),
                                         Real(ppcg_numerical_threshold)));
            // Only scale if the norm is non-negligible; a near-zero
            // column is a converged band whose contribution is harmless.
            sn_scale_all[j] = (sn > Real(ppcg_scaling_threshold))
                            ? double(Real(1) / sn)
                            : 1.0;
        }
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) if (n_dim_ * lcols > ppcg_openmp_work_threshold)
#endif
        for (int j = 0; j < lcols; ++j) {
            for (int ig = 0; ig < n_dim_; ++ig) {
                const Real scale = Real(sn_scale_all[j]);
                x[ idx(ig, j, ld_psi_)] *= scale;
                sx[idx(ig, j, ld_psi_)] *= scale;
                hx[idx(ig, j, ld_psi_)] *= scale;
            }
        }
    };
    scale_to_unit_snorm(subspace.w_l,
                        subspace.sw_l,
                        subspace.hw_l,
                        l);

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

    auto hermitize = [&](std::vector<T>& mat)
    {
        for (int j = 0; j < dim; ++j)
        {
            mat[j + j * dim] = T(std::real(mat[j + j * dim]), 0);
            for (int i = j + 1; i < dim; ++i)
            {
                const T avg = (mat[i + j * dim] + std::conj(mat[j + i * dim]))
                            * Real(0.5);
                mat[i + j * dim] = avg;
                mat[j + i * dim] = std::conj(avg);
            }
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

    gram(subspace.basis.data(), subspace.hbasis.data(), dim, dim, subspace.k, dim);
    gram(subspace.basis.data(), subspace.sbasis.data(), dim, dim, subspace.m, dim);
    hermitize(subspace.k);
    hermitize(subspace.m);
}

// ---------------------------------------------------------------------------
// Solve K v = λ M v (small generalized eigenvalue problem)
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::solve_small_generalized(
    int dim, SmallSubspace& subspace) const
{
    // Try with increasing diagonal shifts; fall back to identity (no update)
    // if the subspace is too ill-conditioned.
    // Save originals; sygvd modifies both matrices in-place before it may
    // fail.
    const std::vector<T> k0 = subspace.k;
    const std::vector<T> m0 = subspace.m;
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

        try
        {
            HermitianLapack<T>::sygvd(dim, subspace.k.data(),
                                      subspace.m.data(),
                                      subspace.eval.data());
            return;
        }
        catch (const std::runtime_error&)
        {
            // Try the next diagonal shift.
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
    SmallSubspace& subspace)
{
    const int dim = 2 * l;
    const T* eigvec = subspace.k.data();

    subspace.psi_new.assign(ld_psi_ * l, T(0));
    subspace.spsi_new.assign(ld_psi_ * l, T(0));
    subspace.hpsi_new.assign(ld_psi_ * l, T(0));

    subspace.coeff_state.resize(dim * l);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (l * l > ppcg_openmp_work_threshold)
#endif
    for (int j = 0; j < l; ++j)
    {
        for (int i = 0; i < l; ++i)
        {
            subspace.coeff_state[i + j * dim] = eigvec[i + j * dim];
            subspace.coeff_state[(l + i) + j * dim] = eigvec[(l + i) + j * dim];
        }
    }

    auto fill_basis = [&](const std::vector<T>& a,
                          const std::vector<T>& b,
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

    fill_basis(subspace.psi_l, subspace.w_l, subspace.basis);
    fill_basis(subspace.spsi_l, subspace.sw_l, subspace.sbasis);
    fill_basis(subspace.hpsi_l, subspace.hw_l, subspace.hbasis);

    combine(subspace.basis, subspace.coeff_state, subspace.psi_new);
    combine(subspace.sbasis, subspace.coeff_state, subspace.spsi_new);
    combine(subspace.hbasis, subspace.coeff_state, subspace.hpsi_new);

    scatter_cols(psi, cols, subspace.psi_new);
    scatter_cols(spsi_.data(), cols, subspace.spsi_new);
    scatter_cols(hpsi_.data(), cols, subspace.hpsi_new);
}

} // namespace hsolver


namespace hsolver {

// ---------------------------------------------------------------------------
// Check S-orthonormality of a column block.
// ---------------------------------------------------------------------------
template <typename T, typename Device>
bool DiagoPPCG<T, Device>::is_s_orthonormal(
    const T* psi, const T* spsi, int ncol) const
{
    const Real orth_tol = Real(ppcg_orthogonality_tolerance_factor)
                        * std::sqrt(std::numeric_limits<Real>::epsilon());
    std::vector<T> gram_s;
    gram(psi, spsi, ncol, ncol, gram_s, ncol);
    for (int j = 0; j < ncol; ++j)
    {
        for (int i = 0; i < ncol; ++i)
        {
            const T sij = gram_s[i + j * ncol];
            const T target = (i == j) ? T(1) : T(0);
            if (std::abs(sij - target) > orth_tol)
            {
                return false;
            }
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// Iterative S-Gram-Schmidt fallback with one reorthogonalization pass.
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::s_gram_schmidt(
    T* psi, T* hpsi, T* spsi, int ncol) const
{
    for (int j = 0; j < ncol; ++j)
    {
        for (int pass = 0; pass < 2; ++pass)
        {
            apply_s_current(psi + j * ld_psi_, spsi + j * ld_psi_, 1);
            for (int k = 0; k < j; ++k)
            {
                T coeff = complex_dot(psi + k * ld_psi_,
                                      spsi + j * ld_psi_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ > ppcg_openmp_work_threshold)
#endif
                for (int ig = 0; ig < n_dim_; ++ig)
                {
                    psi [idx(ig, j, ld_psi_)] -= coeff * psi [idx(ig, k, ld_psi_)];
                    hpsi[idx(ig, j, ld_psi_)] -= coeff * hpsi[idx(ig, k, ld_psi_)];
                    spsi[idx(ig, j, ld_psi_)] -= coeff * spsi[idx(ig, k, ld_psi_)];
                }
            }
        }
        apply_s_current(psi + j * ld_psi_, spsi + j * ld_psi_, 1);
        Real nrm = std::sqrt(std::max(
            gamma_dot(psi + j * ld_psi_, spsi + j * ld_psi_),
            Real(ppcg_numerical_threshold)));
        Real inv_nrm = Real(1) / nrm;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ > ppcg_openmp_work_threshold)
#endif
        for (int ig = 0; ig < n_dim_; ++ig)
        {
            psi [idx(ig, j, ld_psi_)] *= inv_nrm;
            hpsi[idx(ig, j, ld_psi_)] *= inv_nrm;
            spsi[idx(ig, j, ld_psi_)] *= inv_nrm;
        }
    }
}

// ---------------------------------------------------------------------------
// Rayleigh-Ritz: full subspace diagonalization + residual computation
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::rayleigh_ritz(
    T* psi, Real* eigenvalue,
    std::vector<int>& active_cols,
    const std::vector<double>& ethr_band)
{
    // Remember the eigenvalues of the previous step; convergence is measured
    // as the eigenvalue change between successive Rayleigh-Ritz steps.
    eval_prev_.resize(n_band_);
    std::copy(eigenvalue, eigenvalue + n_band_, eval_prev_.begin());

    gram(psi, hpsi_.data(), n_band_, n_band_, rr_hsub_, n_band_);
    gram(psi, spsi_.data(), n_band_, n_band_, rr_ssub_, n_band_);

    bool sygvd_ok = false;
    try
    {
        HermitianLapack<T>::sygvd(n_band_, rr_hsub_.data(), rr_ssub_.data(),
                                  rr_eval_.data());
        sygvd_ok = true;
    }
    catch (const std::runtime_error&)
    {
        // Fallback: diagonal Rayleigh quotients.
        // hsub and ssub may be corrupted by sygvd; re-form them.
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

    lock_epairs(eval_prev_.data(), eigenvalue, ethr_band, active_cols);
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
        residual_trace << "iteration,stage,max_residual\n";
    }
    }
    auto record_residual = [&](int iteration, const char* stage) {
    if (!residual_trace)
    {
        return;
    }
    residual_trace
        << iteration << ','
        << stage << ','
        << max_generalized_residual(hpsi_.data(),
                                    spsi_.data(),
                                    eigenvalue_in,
                                    ld_psi_,
                                    n_dim_,
                                    ncol)
        << '\n';
    };

    // Initialize with Rayleigh-Ritz.
    rayleigh_ritz(psi_in, eigenvalue_in, active_cols, ethr_band);
    // Recompute to keep hpsi/spi consistent with rotated psi.
    apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
    apply_s_current(psi_in, spsi_.data(), ncol);
    record_residual(0, "initial_rr");

    std::vector<T> w_active;
    std::vector<T> sw_active;
    std::vector<T> hw_active;
    w_active.reserve(sz);
    sw_active.reserve(sz);
    hw_active.reserve(sz);
    std::vector<int> cols;
    cols.reserve(std::min(sbsize_, ncol));
    SmallSubspace subspace;

    while (!active_cols.empty() && iter <= maxiter_)
    {
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

        avg_iter += double(nact) / double(ncol);

        // Use the stable 2-block [psi, w] projected subspace.  The
        // preconditioned residual w is normalized to unit S-norm before
        // building the Gram matrix (see build_small_subspace), which
        // keeps M well-conditioned even when residuals are small.

        // Block subspace solve.
        for (int isb = 0; isb < nsb; ++isb)
        {
            const int i0 = isb * sbsize_;
            const int l = std::min(sbsize_, nact - i0);
            cols.assign(active_cols.begin() + i0,
                        active_cols.begin() + i0 + l);

            build_small_subspace(psi_in, cols, subspace);
            solve_small_generalized(2 * l, subspace);
            update_one_block(psi_in, cols, l, subspace);
        }

        // Rayleigh-Ritz after each block update keeps the global subspace
        // synchronized with the updated active vectors.  The block update
        // can otherwise drift into an ill-conditioned basis before the next
        // Ritz rotation.
        rayleigh_ritz(psi_in, eigenvalue_in, active_cols, ethr_band);
        // The Rayleigh-Ritz rotation already keeps hpsi_/spsi_ consistent
        // with the rotated psi up to rounding; re-applying H/S exactly is
        // only needed every rr_step_ iterations to reset the accumulated
        // rounding drift.
        if ((iter % rr_step_) == 0)
        {
            apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
            apply_s_current(psi_in, spsi_.data(), ncol);
        }
        record_residual(iter, "rayleigh_ritz");

        ++iter;
    }

    // Final consistency: ensure hpsi/spi match the converged psi.
    apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
    apply_s_current(psi_in, spsi_.data(), ncol);
    record_residual(iter - 1, "final");
    return avg_iter;
}

} // namespace hsolver

namespace hsolver {

template class DiagoPPCG<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_CPU>;

} // namespace hsolver
