#include "diago_ppcg.h"
#include "hs_operator.h"
#include "kernels/ppcg_kernel_op.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_cell.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include <ATen/kernels/lapack.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace hsolver
{
namespace
{

#ifdef _OPENMP
// Shared scheduling cutoff; an enumerator has no mutable or allocated storage.
enum
{
    openmp_work_threshold = 4096
};
#endif

// A whole-function timer also closes correctly on early return or exception.
class PPCGFunctionTimer
{
  public:
    explicit PPCGFunctionTimer(const char *function) : function_(function)
    {
        ModuleBase::timer::start("DiagoPPCG", function_);
    }

    ~PPCGFunctionTimer() { ModuleBase::timer::end("DiagoPPCG", function_); }

  private:
    const char *function_;
};

} // namespace
} // namespace hsolver

namespace hsolver
{
namespace
{

template <typename Value>
void reduce_pool_values(Value &value)
{
    Parallel_Reduce::reduce_pool(value);
}

template <typename Value>
void reduce_pool_values(Value *value, const int n)
{
    Parallel_Reduce::reduce_pool(value, n);
}

bool all_pool_operations_succeeded(const bool local_success)
{
    int failure_count = local_success ? 0 : 1;
    reduce_pool_values(&failure_count, 1);
    return failure_count == 0;
}

template <typename T, typename Device, typename Real>
Real max_generalized_residual(
    const T *hpsi,
    const T *spsi,
    const Real *eigenvalue,
    int ld,
    int n_dim,
    int ncol)
{
    const PPCGFunctionTimer timer("max_generalized_residual");
    typename PPCGWorkspaceType<T, Device>::type residual(ld * ncol, T(0));
    typename PPCGWorkspaceType<Real, Device>::type values(ncol, Real(0));
    base_device::memory::synchronize_memory_op<Real, Device, base_device::DEVICE_CPU>()(
        values.data(), eigenvalue, ncol);
    PPCGVectorOps<T, Device>::residual(n_dim, ld, ncol, hpsi, spsi, values.data(), residual.data());
    std::vector<double> nrm2_all(ncol, 0.0);
    PPCGVectorOps<T, Device>::column_products(n_dim, ld, ncol, residual.data(), residual.data(), nrm2_all.data());
    reduce_pool_values(nrm2_all.data(), ncol);
    Real max_res = 0;
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

template <typename T, typename Device>
void copy_values(const T *source, const int count, T *destination)
{
    if (count > 0)
    {
        base_device::memory::synchronize_memory_op<T, Device, Device>()(destination, source, count);
    }
}

template <typename T, typename Device>
void zero_values(T *destination, const int count)
{
    if (count > 0)
    {
        base_device::memory::set_memory_op<T, Device>()(destination, 0, count);
    }
}

template <typename T, typename Device, typename Buffer>
void set_zero(Buffer &values)
{
    zero_values<T, Device>(values.data(), int(values.size()));
}

template <typename T, typename Device>
typename PPCGWorkspaceType<T, Device>::type upload_values(const T *values, const int count)
{
    typename PPCGWorkspaceType<T, Device>::type result(count, T(0));
    if (count > 0)
    {
        base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(result.data(), values, count);
    }
    return result;
}

} // anonymous namespace
} // namespace hsolver

namespace hsolver
{
namespace
{

template <typename Scalar>
struct HermitianLapack
{
    using Real = typename container::GetTypeReal<Scalar>::type;
    using Device = container::DEVICE_CPU;

    static void sygvd(
        int n,
        Scalar *a,
        Scalar *b,
        Real *w)
    {
        std::vector<Scalar> eigenvectors(n * n);
        container::kernels::lapack_hegvd<Scalar, Device>()(n, n, a, b, w, eigenvectors.data());
        std::copy(eigenvectors.begin(), eigenvectors.end(), a);
    }
};

} // anonymous namespace
} // namespace hsolver

namespace hsolver
{
namespace
{

bool ppcg_contiguous_cols(const std::vector<int> &cols, int &first)
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

bool ppcg_workspace_dimensions_valid(
    const int ld_psi,
    const int nband,
    const int dim,
    const int sbsize)
{
    if (ld_psi <= 0 || nband <= 0 || dim <= 0 || dim > ld_psi)
    {
        return false;
    }
    const long long limit = std::numeric_limits<int>::max();
    const int bounded_sbsize = std::min(nband, sbsize);
    const int block_size = std::max(1, bounded_sbsize);
    if (static_cast<long long>(nband) > limit / ld_psi || block_size > limit / 3 || nband > limit / 3)
    {
        return false;
    }
    // The largest local basis holds [X, W, P]. Startup correction can hold
    // [X, W], while discarded directions can add up to 2*nband columns during
    // root recovery. The expanded recovery solve therefore reaches 3*nband.
    const long long local_dim = 3LL * block_size;
    const long long recovery_dim = 3LL * nband;
    if (local_dim > limit / ld_psi || recovery_dim > limit / ld_psi)
    {
        return false;
    }
    long long projected_dim = std::max(local_dim, recovery_dim);
    if (block_size < nband)
    {
        const long long startup_dim = 2LL * nband;
        projected_dim = std::max(projected_dim, startup_dim);
    }
    const long long matrix_size = projected_dim * projected_dim;
    const long long linear_workspace = 1 + 6 * projected_dim;
    if (matrix_size > limit || linear_workspace > limit)
    {
        return false;
    }
    // lapack_heevd/hegvd additionally form 1+6*d+2*d*d in int arithmetic.
    return matrix_size <= (limit - linear_workspace) / 2;
}

// =============================================================================
// Constructor
// =============================================================================
template <typename T, typename Device>
DiagoPPCG<T, Device>::DiagoPPCG(
    const Real &diag_thr,
    const int &diag_iter_max,
    const int &sbsize,
    const int &rr_step)
    : maxiter_(diag_iter_max), sbsize_(sbsize), rr_step_(rr_step), diag_thr_(diag_thr)
{
    sbsize_ = std::max(1, sbsize_);
    rr_step_ = std::max(1, rr_step_);
    const Real minimum_diagonalization_threshold = Real(1.0e-14);
    diag_thr_ = std::max(diag_thr_, minimum_diagonalization_threshold);
}

// =============================================================================
// Input validation
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::validate_input(
    const T *psi_in,
    const Real *eigenvalue_in,
    const std::vector<double> &ethr_band,
    const Real *prec) const
{
    std::string local_error;
    try
    {
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
        if (!ppcg_workspace_dimensions_valid(ld_psi_, n_band_, n_dim_, sbsize_))
        {
            throw std::invalid_argument("PPCG: workspace dimensions exceed integer indexing or LAPACK limits.");
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
    catch (const std::invalid_argument &error)
    {
        local_error = error.what();
    }
    // All pool ranks must reject the call before any rank enters an operator
    // callback or Gram-matrix collective, even if only one input is invalid.
    const bool local_success = local_error.empty();
    if (!all_pool_operations_succeeded(local_success))
    {
        if (local_error.empty())
        {
            local_error = "PPCG: invalid input on another rank in the MPI pool.";
        }
        throw std::invalid_argument(local_error);
    }
}

// =============================================================================
// Operator application
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::apply_h(
    const HSOperator<T, Device> &op,
    T *psi_in,
    T *hpsi_out,
    int ncol) const
{
    const PPCGFunctionTimer timer("apply_h");
    op.hpsi(psi_in, hpsi_out, ld_psi_, ncol);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::apply_s(
    const HSOperator<T, Device> &op,
    T *psi_in,
    T *spsi_out,
    int ncol) const
{
    const PPCGFunctionTimer timer("apply_s");
    op.spsi(psi_in, spsi_out, ld_psi_, ncol);
}

// =============================================================================
// Inner product <x|y> (real part only, for Hermitian operators)
// =============================================================================
template <typename T, typename Device>
typename DiagoPPCG<T, Device>::Real DiagoPPCG<T, Device>::gamma_dot(const T *x, const T *y) const
{
    Real result = ModuleBase::dot_real_op<T, Device>()(n_dim_, x, y, false);
    reduce_pool_values(result);
    return result;
}

// =============================================================================
// Gram matrix: out[i, j] = <a_i | b_j>
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::gram(
    const T *mat_a,
    const T *mat_b,
    int ncol_a,
    int ncol_b,
    std::vector<T> &out,
    int ld_out) const
{
    const int reduction_count = ld_out * ncol_b;
    out.resize(reduction_count);
    Workspace projected(reduction_count, T(0));
    const T one = T(1);
    const T zero = T(0);
    ModuleBase::gemm_op<T, Device>()(
        'C', 'N', ncol_a, ncol_b, n_dim_, &one, mat_a, ld_psi_, mat_b, ld_psi_, &zero, projected.data(), ld_out);
    base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>()(
        out.data(), projected.data(), reduction_count);
    reduce_pool_values(out.data(), reduction_count);
}

// =============================================================================
// Column gather: extract selected columns into contiguous storage
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::copy_cols(const T *src, const std::vector<int> &cols, Workspace &dst) const
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
        const T *source = src + first * ld_psi_;
        const int count = ld_psi_ * ncols;
        copy_values<T, Device>(source, count, dst.data());
        return;
    }
    for (int j = 0; j < ncols; ++j)
    {
        const T *source = src + cols[j] * ld_psi_;
        T *destination = dst.data() + j * ld_psi_;
        copy_values<T, Device>(source, ld_psi_, destination);
    }
}

// =============================================================================
// Column scatter: write contiguous storage back into selected columns
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::scatter_cols(T *dst, const std::vector<int> &cols, const Workspace &src) const
{
    const int ncols = int(cols.size());
    if (ncols == 0)
    {
        return;
    }
    int first = 0;
    if (ppcg_contiguous_cols(cols, first))
    {
        T *destination = dst + first * ld_psi_;
        const int count = ld_psi_ * ncols;
        copy_values<T, Device>(src.data(), count, destination);
        return;
    }
    for (int j = 0; j < ncols; ++j)
    {
        const T *source = src.data() + j * ld_psi_;
        T *destination = dst + cols[j] * ld_psi_;
        copy_values<T, Device>(source, ld_psi_, destination);
    }
}

// =============================================================================
// Project x onto vectors orthogonal to S-orthonormal basis
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::project_against(
    const T *basis,
    const T *sbasis,
    const T *hbasis,
    const std::vector<int> &basis_cols,
    Workspace &x,
    Workspace &sx,
    Workspace *hx,
    const std::vector<int> &x_cols) const
{
    const PPCGFunctionTimer timer("project_against");
    if (basis_cols.empty() || x_cols.empty())
    {
        return;
    }

    const int nbasis = int(basis_cols.size());
    const int nx = int(x_cols.size());

    int x_first = 0;
    const bool contiguous_x = ppcg_contiguous_cols(x_cols, x_first);

    Workspace x_l;
    Workspace sx_l;
    Workspace hx_l;
    T *x_data = x.data() + x_first * ld_psi_;
    T *sx_data = sx.data() + x_first * ld_psi_;
    T *hx_data = hx == nullptr ? nullptr : hx->data() + x_first * ld_psi_;
    if (!contiguous_x)
    {
        x_l.reserve(ld_psi_ * nx);
        sx_l.reserve(ld_psi_ * nx);
        copy_cols(x.data(), x_cols, x_l);
        copy_cols(sx.data(), x_cols, sx_l);
        x_data = x_l.data();
        sx_data = sx_l.data();
        if (hx != nullptr)
        {
            copy_cols(hx->data(), x_cols, hx_l);
            hx_data = hx_l.data();
        }
    }

    int basis_first = 0;
    const bool contiguous_basis = ppcg_contiguous_cols(basis_cols, basis_first);

    Workspace basis_l;
    Workspace sbasis_l;
    Workspace hbasis_l;
    const T *basis_data = basis + basis_first * ld_psi_;
    const T *sbasis_data = sbasis + basis_first * ld_psi_;
    const T *hbasis_data = hbasis == nullptr ? nullptr : hbasis + basis_first * ld_psi_;
    if (!contiguous_basis)
    {
        basis_l.reserve(ld_psi_ * nbasis);
        sbasis_l.reserve(ld_psi_ * nbasis);
        copy_cols(basis, basis_cols, basis_l);
        copy_cols(sbasis, basis_cols, sbasis_l);
        basis_data = basis_l.data();
        sbasis_data = sbasis_l.data();
        if (hx != nullptr)
        {
            copy_cols(hbasis, basis_cols, hbasis_l);
            hbasis_data = hbasis_l.data();
        }
    }

    std::vector<T> coeff(nbasis * nx, T(0));
    gram(basis_data, sx_data, nbasis, nx, coeff, nbasis);

    const int coefficient_count = nbasis * nx;
    auto device_coeff = upload_values<T, Device>(coeff.data(), coefficient_count);
    const T minus_one = T(-1);
    const T one = T(1);
    ModuleBase::gemm_op<T, Device>()(
        'N', 'N', n_dim_, nx, nbasis, &minus_one, basis_data, ld_psi_, device_coeff.data(), nbasis, &one, x_data,
        ld_psi_);
    ModuleBase::gemm_op<T, Device>()(
        'N', 'N', n_dim_, nx, nbasis, &minus_one, sbasis_data, ld_psi_, device_coeff.data(), nbasis, &one, sx_data,
        ld_psi_);
    if (hx != nullptr)
    {
        // H(x-X*C) = Hx-HX*C; the same coefficients update all three images.
        ModuleBase::gemm_op<T, Device>()(
            'N', 'N', n_dim_, nx, nbasis, &minus_one, hbasis_data, ld_psi_, device_coeff.data(), nbasis, &one, hx_data,
            ld_psi_);
    }

    if (!contiguous_x)
    {
        scatter_cols(x.data(), x_cols, x_l);
        scatter_cols(sx.data(), x_cols, sx_l);
        if (hx != nullptr)
        {
            scatter_cols(hx->data(), x_cols, hx_l);
        }
    }
}

// =============================================================================
// Project P and refresh its operator images when subtraction loses precision.
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::project_search_directions(
    const HSOperator<T, Device> &op,
    const T *psi,
    const std::vector<int> &all_cols,
    const std::vector<int> &active_cols)
{
    const PPCGFunctionTimer timer("project_search_directions");
    const int nactive = int(active_cols.size());
    std::vector<double> norms(2 * n_band_, 0.0);
    PPCGVectorOps<T, Device>::column_products(n_dim_, ld_psi_, n_band_, p_.data(), p_.data(), norms.data());
    project_against(psi, spsi_.data(), hpsi_.data(), all_cols, p_, sp_, &hp_, active_cols);
    double *projected_norms = norms.data() + n_band_;
    PPCGVectorOps<T, Device>::column_products(n_dim_, ld_psi_, n_band_, p_.data(), p_.data(), projected_norms);
    const int norm_count = int(norms.size());
    reduce_pool_values(norms.data(), norm_count);
    std::vector<int> refresh_cols;
    for (int j = 0; j < nactive; ++j)
    {
        const int column = active_cols[j];
        if (norms[column] == 0.0 || norms[n_band_ + column] == 0.0)
        {
            const int offset = active_cols[j] * ld_psi_;
            T *p_column = p_.data() + offset;
            zero_values<T, Device>(p_column, ld_psi_);
            T *sp_column = sp_.data() + offset;
            zero_values<T, Device>(sp_column, ld_psi_);
            T *hp_column = hp_.data() + offset;
            zero_values<T, Device>(hp_column, ld_psi_);
            continue;
        }
        // A strongly reduced P suffers cancellation in HP-HX*C and SP-SX*C.
        // Normalizing that tiny P would amplify the absolute image error.
        // Refresh only these columns; reduction gives every rank one selection.
        if (norms[n_band_ + column] < 0.5 * norms[column])
        {
            refresh_cols.push_back(active_cols[j]);
        }
    }
    if (!refresh_cols.empty())
    {
        Workspace projected;
        copy_cols(p_.data(), refresh_cols, projected);
        Workspace hprojected(projected.size(), T(0));
        Workspace sprojected(projected.size(), T(0));
        const int refresh_count = int(refresh_cols.size());
        apply_h(op, projected.data(), hprojected.data(), refresh_count);
        apply_s(op, projected.data(), sprojected.data(), refresh_count);
        scatter_cols(hp_.data(), refresh_cols, hprojected);
        scatter_cols(sp_.data(), refresh_cols, sprojected);
    }
}

// =============================================================================
// Preconditioner: x[c] /= max(prec, eps) for each active column c
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::divide_by_preconditioner(
    const std::vector<int> &active_cols, const Real *prec, Workspace &x) const
{
    const int ncols = int(active_cols.size());
    if (ncols == 0)
    {
        return;
    }
    int first = 0;
    if (ppcg_contiguous_cols(active_cols, first))
    {
        T *values = x.data() + first * ld_psi_;
        PPCGVectorOps<T, Device>::precondition(n_dim_, ld_psi_, ncols, prec, values);
        return;
    }
    for (const int column : active_cols)
    {
        T *values = x.data() + column * ld_psi_;
        PPCGVectorOps<T, Device>::precondition(n_dim_, ld_psi_, 1, prec, values);
    }
}

} // namespace hsolver

namespace hsolver
{

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
    const Workspace &residual, const std::vector<double> &ethr_band, std::vector<int> &active_cols)
{
    active_cols.clear();
    active_cols.reserve(n_band_);
    std::vector<double> nrm2_all(n_band_, 0.0);
    PPCGVectorOps<T, Device>::column_products(
        n_dim_, ld_psi_, n_band_, residual.data(), residual.data(), nrm2_all.data());
    reduce_pool_values(nrm2_all.data(), n_band_);
    bool found_unconverged = false;
    max_residual_ = Real(0);
    for (int j = 0; j < n_band_; ++j)
    {
        const Real rnrm = std::sqrt(std::max(Real(nrm2_all[j]), Real(0)));
        const Real thr = std::max(Real(ethr_band[j]), diag_thr_);
        if (std::isfinite(rnrm))
        {
            max_residual_ = std::max(max_residual_, rnrm);
        }
        else
        {
            max_residual_ = std::numeric_limits<Real>::infinity();
        }
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
void DiagoPPCG<T, Device>::scale_to_unit_snorm(
    Workspace &x,
    Workspace &sx,
    Workspace &hx,
    const int ncols) const
{
    const Real squared_snorm_floor = Real(1.0e-30);
    const Real minimum_snorm = Real(1.0e-15);
    std::vector<double> inverse_snorm(ncols, 0.0);
    PPCGVectorOps<T, Device>::column_products(n_dim_, ld_psi_, ncols, x.data(), sx.data(), inverse_snorm.data());
    reduce_pool_values(inverse_snorm.data(), ncols);
    for (int j = 0; j < ncols; ++j)
    {
        const Real snorm = std::sqrt(std::max(Real(inverse_snorm[j]), squared_snorm_floor));
        inverse_snorm[j] = snorm > minimum_snorm ? double(Real(1) / snorm) : 1.0;
    }
    auto scales = upload_values<double, Device>(inverse_snorm.data(), ncols);
    PPCGVectorOps<T, Device>::scale_columns(n_dim_, ld_psi_, ncols, scales.data(), x.data(), sx.data(), hx.data());
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::hermitize_projected(std::vector<T> &matrix, const int dim) const
{
    for (int j = 0; j < dim; ++j)
    {
        matrix[j + j * dim] = T(std::real(matrix[j + j * dim]), 0);
        for (int i = j + 1; i < dim; ++i)
        {
            const T value = (matrix[i + j * dim] + std::conj(matrix[j + i * dim])) * Real(0.5);
            matrix[i + j * dim] = value;
            matrix[j + i * dim] = std::conj(value);
        }
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::insert_startup_gram_block(
    const T *left,
    const T *right,
    const int row_count,
    const int column_count,
    const int row_offset,
    const int column_offset,
    const int matrix_dim,
    std::vector<T> &matrix,
    std::vector<T> &workspace) const
{
    gram(left, right, row_count, column_count, workspace, row_count);
    for (int j = 0; j < column_count; ++j)
    {
        for (int i = 0; i < row_count; ++i)
        {
            const int matrix_index = row_offset + i + (column_offset + j) * matrix_dim;
            matrix[matrix_index] = workspace[i + j * row_count];
        }
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::build_startup_global_subspace(
    const T *psi,
    Workspace &w_active,
    Workspace &sw_active,
    Workspace &hw_active,
    const int active_count,
    SmallSubspace &subspace)
{
    scale_to_unit_snorm(w_active, sw_active, hw_active, active_count);

    const int state_count = n_band_;
    const int correction_offset = state_count;
    const int local_dim = state_count + active_count;
    const int matrix_size = local_dim * local_dim;
    subspace.k.assign(matrix_size, T(0));
    subspace.m.assign(matrix_size, T(0));
    subspace.eval.resize(local_dim);
    std::vector<T> gram_workspace;

    insert_startup_gram_block(psi, hpsi_.data(), state_count, state_count, 0, 0, local_dim, subspace.k, gram_workspace);
    insert_startup_gram_block(
        psi, hw_active.data(), state_count, active_count, 0, correction_offset, local_dim, subspace.k, gram_workspace);
    insert_startup_gram_block(
        w_active.data(), hpsi_.data(), active_count, state_count, correction_offset, 0, local_dim, subspace.k,
        gram_workspace);
    insert_startup_gram_block(
        w_active.data(), hw_active.data(), active_count, active_count, correction_offset, correction_offset, local_dim,
        subspace.k, gram_workspace);
    insert_startup_gram_block(psi, spsi_.data(), state_count, state_count, 0, 0, local_dim, subspace.m, gram_workspace);
    insert_startup_gram_block(
        psi, sw_active.data(), state_count, active_count, 0, correction_offset, local_dim, subspace.m, gram_workspace);
    insert_startup_gram_block(
        w_active.data(), spsi_.data(), active_count, state_count, correction_offset, 0, local_dim, subspace.m,
        gram_workspace);
    insert_startup_gram_block(
        w_active.data(), sw_active.data(), active_count, active_count, correction_offset, correction_offset, local_dim,
        subspace.m, gram_workspace);

    hermitize_projected(subspace.k, local_dim);
    hermitize_projected(subspace.m, local_dim);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::update_startup_global_subspace(
    T *psi,
    const Workspace &w_active,
    const Workspace &sw_active,
    const Workspace &hw_active,
    const int active_count,
    const SmallSubspace &subspace)
{
    const int local_dim = n_band_ + active_count;
    const int coefficient_count = local_dim * n_band_;
    auto device_coeff = upload_values<T, Device>(subspace.k.data(), coefficient_count);
    const T *state_coefficients = device_coeff.data();
    const T *correction_coefficients = state_coefficients + n_band_;
    const T one = T(1);
    const T zero = T(0);

    set_zero<T, Device>(w_);
    set_zero<T, Device>(sw_);
    set_zero<T, Device>(hw_);
    set_zero<T, Device>(p_);
    set_zero<T, Device>(sp_);
    set_zero<T, Device>(hp_);

    auto combine_state = [&](const T *state, const T *correction, T *output)
    {
        ModuleBase::gemm_op<T, Device>()(
            'N', 'N', n_dim_, n_band_, n_band_, &one, state, ld_psi_, state_coefficients, local_dim, &zero, output,
            ld_psi_);
        ModuleBase::gemm_op<T, Device>()(
            'N', 'N', n_dim_, n_band_, active_count, &one, correction, ld_psi_, correction_coefficients, local_dim,
            &one, output, ld_psi_);
    };
    auto combine_direction = [&](const T *correction, T *output)
    {
        ModuleBase::gemm_op<T, Device>()(
            'N', 'N', n_dim_, n_band_, active_count, &one, correction, ld_psi_, correction_coefficients, local_dim,
            &zero, output, ld_psi_);
    };

    combine_state(psi, w_active.data(), w_.data());
    combine_state(spsi_.data(), sw_active.data(), sw_.data());
    combine_state(hpsi_.data(), hw_active.data(), hw_.data());
    combine_direction(w_active.data(), p_.data());
    combine_direction(sw_active.data(), sp_.data());
    combine_direction(hw_active.data(), hp_.data());

    copy_values<T, Device>(w_.data(), int(w_.size()), psi);
    copy_values<T, Device>(sw_.data(), int(sw_.size()), spsi_.data());
    copy_values<T, Device>(hw_.data(), int(hw_.size()), hpsi_.data());
}

// ---------------------------------------------------------------------------
// Build K = V^H H V and M = V^H S V where V = [psi, w]
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::build_small_subspace(
    const T *psi,
    const std::vector<int> &cols,
    int nblk,
    SmallSubspace &subspace) const
{
    const PPCGFunctionTimer timer("build_small_subspace");
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
    scale_to_unit_snorm(subspace.w_l, subspace.sw_l, subspace.hw_l, l);
    if (nblk >= 3)
    {
        scale_to_unit_snorm(subspace.p_l, subspace.sp_l, subspace.hp_l, l);
    }

    auto copy_block = [&](const Workspace &src, const int col0, Workspace &dst)
    {
        T *destination = dst.data() + col0 * ld_psi_;
        const int count = ld_psi_ * l;
        copy_values<T, Device>(src.data(), count, destination);
    };

    const int basis_size = ld_psi_ * dim;
    subspace.basis.resize(basis_size);
    subspace.hbasis.resize(basis_size);
    subspace.sbasis.resize(basis_size);
    copy_block(subspace.psi_l, 0, subspace.basis);
    copy_block(subspace.hpsi_l, 0, subspace.hbasis);
    copy_block(subspace.spsi_l, 0, subspace.sbasis);
    copy_block(subspace.w_l, l, subspace.basis);
    copy_block(subspace.hw_l, l, subspace.hbasis);
    copy_block(subspace.sw_l, l, subspace.sbasis);
    if (nblk >= 3)
    {
        const int p_block_offset = 2 * l;
        copy_block(subspace.p_l, p_block_offset, subspace.basis);
        copy_block(subspace.hp_l, p_block_offset, subspace.hbasis);
        copy_block(subspace.sp_l, p_block_offset, subspace.sbasis);
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
int DiagoPPCG<T, Device>::solve_small_generalized(const int dim, const int nstates, SmallSubspace &subspace) const
{
    const PPCGFunctionTimer timer("solve_small_generalized");
    const std::vector<T> k0 = subspace.k;
    const std::vector<T> m0 = subspace.m;

    // A nearly dependent local basis makes the generalized solve unstable
    // even when LAPACK accepts M as positive definite.  Work in the retained
    // eigenspace of M and map the Ritz vectors back to the original basis.
    std::vector<T> overlap_eigenvectors = m0;
    std::vector<Real> overlap_eigenvalues(dim);
    bool local_rank_truncation_succeeded = false;
    int retained_count = 0;
    try
    {
        container::kernels::lapack_heevd<T, container::DEVICE_CPU>()(
            dim, overlap_eigenvectors.data(), dim, overlap_eigenvalues.data());

        const Real largest_overlap_eigenvalue = std::max(overlap_eigenvalues.back(), Real(0));
        const Real precision_floor = Real(10) * std::numeric_limits<Real>::epsilon();
        // Discard numerically dependent directions in the local overlap spectrum.
        const Real minimum_relative_overlap = Real(1.0e-8);
        const Real relative_rank_threshold = std::max(minimum_relative_overlap, precision_floor);
        const Real overlap_cutoff = largest_overlap_eigenvalue * relative_rank_threshold;
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
                const Real scale = Real(1) / std::sqrt(overlap_eigenvalues[source_index]);
                for (int row = 0; row < dim; ++row)
                {
                    retained_basis[row + j * dim] = overlap_eigenvectors[row + source_index * dim] * scale;
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
            ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()(
                'N', 'N', dim, retained_dim, dim, &one, k0.data(), dim, retained_basis.data(), dim, &zero,
                projected_hamiltonian.data(), dim);
            std::vector<T> reduced_hamiltonian(retained_dim * retained_dim, T(0));
            ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()(
                'C', 'N', retained_dim, retained_dim, dim, &one, retained_basis.data(), dim,
                projected_hamiltonian.data(), dim, &zero, reduced_hamiltonian.data(), retained_dim);
            for (int j = 0; j < retained_dim; ++j)
            {
                const Real diagonal_value = std::real(reduced_hamiltonian[j + j * retained_dim]);
                reduced_hamiltonian[j + j * retained_dim] = T(diagonal_value, Real(0));
                for (int i = 0; i < j; ++i)
                {
                    const T value = Real(0.5) * (reduced_hamiltonian[i + j * retained_dim] +
                                                    std::conj(reduced_hamiltonian[j + i * retained_dim]));
                    reduced_hamiltonian[i + j * retained_dim] = value;
                    reduced_hamiltonian[j + i * retained_dim] = std::conj(value);
                }
            }

            std::vector<Real> reduced_eigenvalues(retained_dim);
            container::kernels::lapack_heevd<T, container::DEVICE_CPU>()(
                retained_dim, reduced_hamiltonian.data(), retained_dim, reduced_eigenvalues.data());

            std::fill(subspace.k.begin(), subspace.k.end(), T(0));
            for (int state = 0; state < retained_dim; ++state)
            {
                subspace.eval[state] = reduced_eigenvalues[state];
            }
            ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()(
                'N', 'N', dim, retained_dim, retained_dim, &one, retained_basis.data(), dim, reduced_hamiltonian.data(),
                retained_dim, &zero, subspace.k.data(), dim);
            local_rank_truncation_succeeded = true;
            retained_count = retained_dim;
        }
    }
    catch (const std::runtime_error &)
    {
        // The shifted generalized solve below remains the fallback.
    }
    if (all_pool_operations_succeeded(local_rank_truncation_succeeded))
    {
        return retained_count;
    }

    // Try with increasing diagonal shifts; fall back to identity (no update)
    // if the subspace is too ill-conditioned.  sygvd modifies both matrices
    // in-place before it may fail.
    const Real shifts[] = {Real(0.0), Real(1.0e-10), Real(1.0e-8), Real(1.0e-6)};
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
            HermitianLapack<T>::sygvd(dim, subspace.k.data(), subspace.m.data(), subspace.eval.data());
            local_sygvd_succeeded = true;
        }
        catch (const std::runtime_error &)
        {
            // Try the next diagonal shift.
        }
        if (all_pool_operations_succeeded(local_sygvd_succeeded))
        {
            return dim;
        }
    }
    // All attempts failed — set eigenvectors to identity (no update).
    const Real minimum_overlap = Real(1.0e-30);
    std::fill(subspace.k.begin(), subspace.k.end(), T(0));
    for (int i = 0; i < dim; ++i)
    {
        subspace.k[i + i * dim] = T(1);
        subspace.eval[i] =
            Real(std::real(k0[i + i * dim])) / std::max(Real(std::real(m0[i + i * dim])), minimum_overlap);
    }
    return 0;
}

// Keep valid local Ritz directions that a block-size restriction would discard.
template <typename T, typename Device>
void DiagoPPCG<T, Device>::collect_discarded_ritz(
    const SmallSubspace &subspace,
    const int kept_count,
    const int valid_count,
    OperatorBlock &candidates) const
{
    const PPCGFunctionTimer timer("collect_discarded_ritz");
    const int count = valid_count - kept_count;
    if (count <= 0)
    {
        return;
    }
    const int dim = int(subspace.eval.size());
    const int offset = int(candidates.vectors.size());
    const T *host_coefficients = subspace.k.data() + kept_count * dim;
    const int coefficient_count = dim * count;
    auto device_coeff = upload_values<T, Device>(host_coefficients, coefficient_count);
    const T *coefficients = device_coeff.data();
    const T one(1);
    const T zero(0);
    const auto append = [&](const Workspace &basis, Workspace &output)
    {
        output.resize(offset + ld_psi_ * count, T(0));
        T *destination = output.data() + offset;
        ModuleBase::gemm_op<T, Device>()(
            'N', 'N', n_dim_, count, dim, &one, basis.data(), ld_psi_, coefficients, dim, &zero, destination, ld_psi_);
    };
    append(subspace.basis, candidates.vectors);
    append(subspace.hbasis, candidates.hvectors);
    append(subspace.sbasis, candidates.svectors);
}

// Preserve X exactly as a subspace: rank decisions apply only to its complement.
template <typename T, typename Device>
void DiagoPPCG<T, Device>::project_candidate_complement(
    const HSOperator<T, Device> &op, const T *psi, OperatorBlock &candidates) const
{
    const PPCGFunctionTimer timer("project_candidate_complement");
    const int count = int(candidates.vectors.size()) / ld_psi_;
    std::vector<int> basis_cols(n_band_);
    std::iota(basis_cols.begin(), basis_cols.end(), 0);
    std::vector<int> candidate_cols(count);
    std::iota(candidate_cols.begin(), candidate_cols.end(), 0);
    std::vector<double> norms(2 * count, 0.0);
    PPCGVectorOps<T, Device>::column_products(
        n_dim_, ld_psi_, count, candidates.vectors.data(), candidates.vectors.data(), norms.data());
    // Two passes suppress overlap left by cancellation against the retained X.
    for (int pass = 0; pass < 2; ++pass)
    {
        project_against(
            psi, spsi_.data(), hpsi_.data(), basis_cols, candidates.vectors, candidates.svectors, &candidates.hvectors,
            candidate_cols);
    }
    double *projected_norms = norms.data() + count;
    PPCGVectorOps<T, Device>::column_products(
        n_dim_, ld_psi_, count, candidates.vectors.data(), candidates.vectors.data(), projected_norms);
    const int norm_count = int(norms.size());
    reduce_pool_values(norms.data(), norm_count);
    const Real precision_floor = Real(10) * std::numeric_limits<Real>::epsilon();
    const Real minimum_retained_norm_ratio = Real(1.0e-8);
    const Real rank_threshold = std::max(minimum_retained_norm_ratio, precision_floor);
    std::vector<int> retained;
    std::vector<int> refresh;
    for (int column = 0; column < count; ++column)
    {
        const double before = norms[column];
        const double after = norms[count + column];
        if (!std::isfinite(before) || !std::isfinite(after) || after <= rank_threshold * before)
        {
            continue;
        }
        retained.push_back(column);
        if (after < 0.5 * before)
        {
            refresh.push_back(column);
        }
    }
    // Bound each operator batch by the number of target bands.
    for (int first = 0; first < int(refresh.size()); first += n_band_)
    {
        const int last = std::min(first + n_band_, int(refresh.size()));
        const std::vector<int> columns(refresh.begin() + first, refresh.begin() + last);
        Workspace vectors;
        copy_cols(candidates.vectors.data(), columns, vectors);
        Workspace hvectors(vectors.size(), T(0));
        Workspace svectors(vectors.size(), T(0));
        const int batch_count = int(columns.size());
        apply_h(op, vectors.data(), hvectors.data(), batch_count);
        apply_s(op, vectors.data(), svectors.data(), batch_count);
        scatter_cols(candidates.hvectors.data(), columns, hvectors);
        scatter_cols(candidates.svectors.data(), columns, svectors);
    }
    OperatorBlock projected;
    copy_cols(candidates.vectors.data(), retained, projected.vectors);
    copy_cols(candidates.hvectors.data(), retained, projected.hvectors);
    copy_cols(candidates.svectors.data(), retained, projected.svectors);
    candidates = std::move(projected);
}

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::recover_competing_ritz(
    const HSOperator<T, Device> &op,
    T *psi,
    Real *eigenvalue,
    const bool ordered_ritz,
    std::vector<int> &active_cols,
    const std::vector<double> &ethr_band,
    OperatorBlock &candidates)
{
    const PPCGFunctionTimer timer("recover_competing_ritz");
    if (candidates.vectors.empty())
    {
        return false;
    }
    Real upper_bound = eigenvalue[n_band_ - 1];
    if (!ordered_ritz)
    {
        // Between Ritz steps a column Rayleigh quotient is not an upper bound
        // on the largest target Ritz value. Use the Hermitian row-sum bound.
        upper_bound = std::numeric_limits<Real>::lowest();
        for (int row = 0; row < n_band_; ++row)
        {
            Real row_bound = std::real(rr_hsub_[row + row * n_band_]);
            for (int column = 0; column < n_band_; ++column)
            {
                if (column != row)
                {
                    row_bound += std::abs(rr_hsub_[row + column * n_band_]);
                }
            }
            upper_bound = std::max(upper_bound, row_bound);
        }
    }
    project_candidate_complement(op, psi, candidates);
    const int count = int(candidates.vectors.size()) / ld_psi_;
    if (count == 0)
    {
        return false;
    }
    std::vector<double> moments(2 * count, 0.0);
    PPCGVectorOps<T, Device>::column_products(
        n_dim_, ld_psi_, count, candidates.vectors.data(), candidates.svectors.data(), moments.data());
    double *energies = moments.data() + count;
    PPCGVectorOps<T, Device>::column_products(
        n_dim_, ld_psi_, count, candidates.vectors.data(), candidates.hvectors.data(), energies);
    const int moment_count = int(moments.size());
    reduce_pool_values(moments.data(), moment_count);
    const Real scale = std::max(Real(1), std::abs(upper_bound));
    const Real margin = Real(32) * std::numeric_limits<Real>::epsilon() * scale;
    std::vector<int> competing;
    for (int column = 0; column < count; ++column)
    {
        const double norm = moments[column];
        const double energy = moments[count + column];
        if (!std::isfinite(norm) || !std::isfinite(energy) || norm <= 0.0)
        {
            continue;
        }
        // Use the physical quotient, not an eigenvalue of the shifted-M
        // fallback; they need not represent the same projected problem.
        if (energy / norm < upper_bound - margin)
        {
            competing.push_back(column);
        }
    }
    if (competing.empty())
    {
        return false;
    }
    OperatorBlock selected;
    copy_cols(candidates.vectors.data(), competing, selected.vectors);
    copy_cols(candidates.hvectors.data(), competing, selected.hvectors);
    copy_cols(candidates.svectors.data(), competing, selected.svectors);
    // The selected copies own the surviving directions; release the larger
    // candidate buffers before constructing their orthogonal basis.
    candidates = OperatorBlock{};
    const int selected_count = int(competing.size());
    scale_to_unit_snorm(selected.vectors, selected.svectors, selected.hvectors, selected_count);
    std::vector<T> overlap;
    gram(selected.vectors.data(), selected.svectors.data(), selected_count, selected_count, overlap, selected_count);
    hermitize_projected(overlap, selected_count);
    std::vector<Real> norms(selected_count);
    bool local_success = false;
    try
    {
        container::kernels::lapack_heevd<T, container::DEVICE_CPU>()(
            selected_count, overlap.data(), selected_count, norms.data());
        local_success = true;
    }
    catch (const std::runtime_error &)
    {
    }
    if (!all_pool_operations_succeeded(local_success))
    {
        throw std::runtime_error("PPCG: competing-direction overlap diagonalization failed.");
    }
    const Real precision_floor = Real(10) * std::numeric_limits<Real>::epsilon();
    const Real minimum_relative_overlap = Real(1.0e-8);
    const Real rank_threshold = std::max(minimum_relative_overlap, precision_floor);
    const Real cutoff = norms.back() * rank_threshold;
    int retained_count = 0;
    for (int column = 0; column < selected_count; ++column)
    {
        if (norms[column] <= cutoff)
        {
            continue;
        }
        const Real inverse_norm = Real(1) / std::sqrt(norms[column]);
        for (int row = 0; row < selected_count; ++row)
        {
            overlap[row + retained_count * selected_count] = overlap[row + column * selected_count] * inverse_norm;
        }
        ++retained_count;
    }
    if (retained_count == 0)
    {
        throw std::runtime_error("PPCG: competing-direction complement has zero numerical rank.");
    }
    const int coefficient_count = selected_count * retained_count;
    auto device_coeff = upload_values<T, Device>(overlap.data(), coefficient_count);
    OperatorBlock orthogonal;
    const T one(1);
    const T zero(0);
    const auto transform = [&](const Workspace &input, Workspace &output)
    {
        output.assign(ld_psi_ * retained_count, T(0));
        ModuleBase::gemm_op<T, Device>()(
            'N', 'N', n_dim_, retained_count, selected_count, &one, input.data(), ld_psi_, device_coeff.data(),
            selected_count, &zero, output.data(), ld_psi_);
    };
    transform(selected.vectors, orthogonal.vectors);
    transform(selected.hvectors, orthogonal.hvectors);
    transform(selected.svectors, orthogonal.svectors);
    selected = OperatorBlock{};
    SmallSubspace expanded;
    build_startup_global_subspace(
        psi, orthogonal.vectors, orthogonal.svectors, orthogonal.hvectors, retained_count, expanded);
    const int expanded_dim = n_band_ + retained_count;
    local_success = false;
    try
    {
        // No rank truncation of [X,Q]: X is retained intact and Q has already
        // been rank-revealed independently in its S-orthogonal complement.
        HermitianLapack<T>::sygvd(expanded_dim, expanded.k.data(), expanded.m.data(), expanded.eval.data());
        local_success = true;
    }
    catch (const std::runtime_error &)
    {
    }
    if (!all_pool_operations_succeeded(local_success))
    {
        throw std::runtime_error("PPCG: competing-direction Rayleigh-Ritz solve failed.");
    }
    update_startup_global_subspace(
        psi, orthogonal.vectors, orthogonal.svectors, orthogonal.hvectors, retained_count, expanded);
    std::copy(expanded.eval.begin(), expanded.eval.begin() + n_band_, eigenvalue);
    set_zero<T, Device>(p_);
    set_zero<T, Device>(sp_);
    set_zero<T, Device>(hp_);
    compute_eigenpair_residual(eigenvalue);
    lock_epairs(w_, ethr_band, active_cols);
    ++rayleigh_ritz_count_;
    ++root_recovery_count_;
    return true;
}

// ---------------------------------------------------------------------------
// Update wavefunctions from small subspace eigenvectors
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::update_one_block(
    T *psi,
    const std::vector<int> &cols,
    int l,
    int nblk,
    SmallSubspace &subspace)
{
    const PPCGFunctionTimer timer("update_one_block");
    const int dim = nblk * l;
    const T *eigvec = subspace.k.data();

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
#pragma omp parallel for schedule(static) if (l * l > openmp_work_threshold)
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

    auto fill_basis = [&](const Workspace &a, const Workspace &b, const Workspace &c, Workspace &basis)
    {
        const int count = ld_psi_ * l;
        basis.resize(ld_psi_ * dim);
        T *second = basis.data() + count;
        copy_values<T, Device>(a.data(), count, basis.data());
        copy_values<T, Device>(b.data(), count, second);
        if (nblk >= 3)
        {
            T *third = basis.data() + 2 * count;
            copy_values<T, Device>(c.data(), count, third);
        }
    };

    auto combine = [&](const Workspace &basis, const Workspace &coeff, Workspace &out)
    {
        const T one = T(1);
        const T zero = T(0);
        ModuleBase::gemm_op<T, Device>()(
            'N', 'N', n_dim_, l, dim, &one, basis.data(), ld_psi_, coeff.data(), dim, &zero, out.data(), ld_psi_);
    };

    fill_basis(subspace.psi_l, subspace.w_l, subspace.p_l, subspace.basis);
    fill_basis(subspace.spsi_l, subspace.sw_l, subspace.sp_l, subspace.sbasis);
    fill_basis(subspace.hpsi_l, subspace.hw_l, subspace.hp_l, subspace.hbasis);

    const int coefficient_count = dim * l;
    auto state_coeff = upload_values<T, Device>(subspace.coeff_state.data(), coefficient_count);
    auto direction_coeff = upload_values<T, Device>(subspace.coeff_p.data(), coefficient_count);

    combine(subspace.basis, state_coeff, subspace.psi_new);
    combine(subspace.sbasis, state_coeff, subspace.spsi_new);
    combine(subspace.hbasis, state_coeff, subspace.hpsi_new);
    combine(subspace.basis, direction_coeff, subspace.p_new);
    combine(subspace.sbasis, direction_coeff, subspace.sp_new);
    combine(subspace.hbasis, direction_coeff, subspace.hp_new);

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
    T *psi,
    Real *eigenvalue,
    std::vector<int> &active_cols,
    const std::vector<double> &ethr_band)
{
    const PPCGFunctionTimer timer("rayleigh_ritz");
    ++rayleigh_ritz_count_;
    gram(psi, hpsi_.data(), n_band_, n_band_, rr_hsub_, n_band_);
    gram(psi, spsi_.data(), n_band_, n_band_, rr_ssub_, n_band_);

    bool local_sygvd_ok = false;
    try
    {
        HermitianLapack<T>::sygvd(n_band_, rr_hsub_.data(), rr_ssub_.data(), rr_eval_.data());
        local_sygvd_ok = true;
    }
    catch (const std::runtime_error &)
    {
        // The collective result below selects one common fallback on all ranks.
    }
    const bool sygvd_ok = all_pool_operations_succeeded(local_sygvd_ok);
    if (!sygvd_ok)
    {
        // LAPACK may overwrite the projected matrices before failing.  Re-form
        // them on every rank so all ranks enter the same fallback path.
        const Real minimum_overlap = Real(1.0e-30);
        gram(psi, hpsi_.data(), n_band_, n_band_, rr_hsub_, n_band_);
        gram(psi, spsi_.data(), n_band_, n_band_, rr_ssub_, n_band_);
        for (int ii = 0; ii < n_band_; ++ii)
        {
            rr_eval_[ii] = Real(std::real(rr_hsub_[ii + ii * n_band_])) /
                           std::max(Real(std::real(rr_ssub_[ii + ii * n_band_])), minimum_overlap);
        }
    }

    if (sygvd_ok)
    {
        const int coefficient_count = n_band_ * n_band_;
        auto device_coeff = upload_values<T, Device>(rr_hsub_.data(), coefficient_count);
        const int sz = ld_psi_ * n_band_;
        copy_values<T, Device>(psi, sz, w_.data());
        copy_values<T, Device>(spsi_.data(), int(spsi_.size()), sw_.data());
        copy_values<T, Device>(hpsi_.data(), int(hpsi_.size()), hw_.data());

        zero_values<T, Device>(psi, sz);
        set_zero<T, Device>(spsi_);
        set_zero<T, Device>(hpsi_);

        const T one = T(1);
        const T zero = T(0);
        ModuleBase::gemm_op<T, Device>()(
            'N', 'N', n_dim_, n_band_, n_band_, &one, w_.data(), ld_psi_, device_coeff.data(), n_band_, &zero, psi,
            ld_psi_);
        ModuleBase::gemm_op<T, Device>()(
            'N', 'N', n_dim_, n_band_, n_band_, &one, sw_.data(), ld_psi_, device_coeff.data(), n_band_, &zero,
            spsi_.data(), ld_psi_);
        ModuleBase::gemm_op<T, Device>()(
            'N', 'N', n_dim_, n_band_, n_band_, &one, hw_.data(), ld_psi_, device_coeff.data(), n_band_, &zero,
            hpsi_.data(), ld_psi_);

        for (int j = 0; j < n_band_; ++j)
        {
            eigenvalue[j] = rr_eval_[j];
        }
        // Keep each conjugate direction attached to the same Ritz column
        // when the global solve mixes or permutes the retained states.
        const auto rotate_direction = [&](Workspace &values)
        {
            copy_values<T, Device>(values.data(), int(values.size()), w_.data());
            set_zero<T, Device>(values);
            ModuleBase::gemm_op<T, Device>()(
                'N', 'N', n_dim_, n_band_, n_band_, &one, w_.data(), ld_psi_, device_coeff.data(), n_band_, &zero,
                values.data(), ld_psi_);
        };
        rotate_direction(p_);
        rotate_direction(sp_);
        rotate_direction(hp_);
    }
    else
    {
        // No rotation: just update eigenvalues with Rayleigh quotients.
        for (int j = 0; j < n_band_; ++j)
        {
            eigenvalue[j] = rr_eval_[j];
        }
    }

    compute_eigenpair_residual(eigenvalue);

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
        max_residual_ =
            max_generalized_residual<T, Device>(hpsi_.data(), spsi_.data(), eigenvalue, ld_psi_, n_dim_, n_band_);
    }
    return sygvd_ok;
}

// Refine one retained pair during stagnation without rotating the other states.
// The projected correction preserves S-orthogonality to the other columns.
// A residual-minimizing step and eigenvalue correction solve a two-variable
// least-squares problem. Only a normalized candidate with a smaller freshly
// evaluated residual and no resolved energy increase replaces the original.
template <typename T, typename Device>
bool DiagoPPCG<T, Device>::refine_stalled_residual(
    const HSOperator<T, Device> &op,
    T *psi_in,
    Real *eigenvalue_in,
    const std::vector<double> &ethr_band,
    const std::vector<int> &all_cols,
    const int band_hint)
{
    const PPCGFunctionTimer timer("refine_stalled_residual");
    const int ncol = n_band_;
    const int sz = ld_psi_ * ncol;
    std::vector<double> old_norms(ncol, 0);
    PPCGVectorOps<T, Device>::column_products(n_dim_, ld_psi_, ncol, w_.data(), w_.data(), old_norms.data());
    reduce_pool_values(old_norms.data(), ncol);
    int selected = -1;
    double largest_relative_residual = 1.0;
    for (int band = 0; band < ncol; ++band)
    {
        const double threshold = std::max(double(diag_thr_), ethr_band[band]);
        const double relative_residual = old_norms[band] / (threshold * threshold);
        if (std::isfinite(relative_residual) && relative_residual > largest_relative_residual)
        {
            selected = band;
            largest_relative_residual = relative_residual;
        }
    }
    if (band_hint >= 0 && band_hint < ncol)
    {
        selected = band_hint;
    }
    if (selected < 0 || !std::isfinite(old_norms[selected]) || old_norms[selected] == 0.0)
    {
        return false;
    }
    const std::vector<int> selected_cols(1, selected);
    Workspace original;
    copy_cols(psi_in, selected_cols, original);
    const Real original_value = eigenvalue_in[selected];
    const Real margin = Real(16) * std::numeric_limits<Real>::epsilon() * std::max(Real(1), std::abs(original_value));
    const auto refine_value = [&]()
    {
        // An exact state can retain a few ulps of eigenvalue error. Its residual
        // is parallel to Sx, so projecting a vector correction would remove it.
        const int value_moment_count = 2 * ncol;
        std::vector<double> value_moments(value_moment_count, 0.0);
        double *value_overlap_norms = value_moments.data() + ncol;
        PPCGVectorOps<T, Device>::column_products(n_dim_, ld_psi_, ncol, w_.data(), spsi_.data(), value_moments.data());
        PPCGVectorOps<T, Device>::column_products(
            n_dim_, ld_psi_, ncol, spsi_.data(), spsi_.data(), value_overlap_norms);
        reduce_pool_values(value_moments.data(), value_moment_count);
        const double value_overlap_norm = value_overlap_norms[selected];
        if (value_overlap_norm > 0.0 && std::isfinite(value_overlap_norm))
        {
            const Real correction = Real(value_moments[selected] / value_overlap_norm);
            const Real candidate_value = original_value + correction;
            const bool above_previous = selected == 0 || candidate_value >= eigenvalue_in[selected - 1] - margin;
            const bool below_next = selected + 1 == ncol || candidate_value <= eigenvalue_in[selected + 1] + margin;
            if (std::isfinite(correction) && std::abs(correction) <= margin && candidate_value != original_value &&
                above_previous && below_next)
            {
                eigenvalue_in[selected] = candidate_value;
                compute_eigenpair_residual(eigenvalue_in);
                std::vector<double> value_residual_norms(ncol, 0.0);
                PPCGVectorOps<T, Device>::column_products(
                    n_dim_, ld_psi_, ncol, w_.data(), w_.data(), value_residual_norms.data());
                reduce_pool_values(value_residual_norms.data(), ncol);
                if (std::isfinite(value_residual_norms[selected]) &&
                    value_residual_norms[selected] < old_norms[selected])
                {
                    return true;
                }
                eigenvalue_in[selected] = original_value;
                compute_eigenpair_residual(eigenvalue_in);
            }
        }
        return false;
    };
    Workspace direction(sz, T(0));
    copy_values<T, Device>(w_.data(), sz, direction.data());
    divide_by_preconditioner(all_cols, device_prec_.data(), direction);
    Workspace h_direction(sz, T(0));
    Workspace s_direction(sz, T(0));
    Workspace response(sz, T(0));
    apply_s(op, direction.data(), s_direction.data(), ncol);
    project_against(psi_in, spsi_.data(), nullptr, all_cols, direction, s_direction, nullptr, all_cols);
    apply_h(op, direction.data(), h_direction.data(), ncol);
    auto device_values = upload_values<Real, Device>(eigenvalue_in, ncol);
    PPCGVectorOps<T, Device>::residual(
        n_dim_, ld_psi_, ncol, h_direction.data(), s_direction.data(), device_values.data(), response.data());
    const int moment_count = 5 * ncol;
    std::vector<double> moments(moment_count, 0.0);
    double *response_norms = moments.data() + ncol;
    double *residual_overlaps = moments.data() + 2 * ncol;
    double *response_overlaps = moments.data() + 3 * ncol;
    double *overlap_norms = moments.data() + 4 * ncol;
    PPCGVectorOps<T, Device>::column_products(n_dim_, ld_psi_, ncol, w_.data(), response.data(), moments.data());
    PPCGVectorOps<T, Device>::column_products(n_dim_, ld_psi_, ncol, response.data(), response.data(), response_norms);
    PPCGVectorOps<T, Device>::column_products(n_dim_, ld_psi_, ncol, w_.data(), spsi_.data(), residual_overlaps);
    PPCGVectorOps<T, Device>::column_products(n_dim_, ld_psi_, ncol, response.data(), spsi_.data(), response_overlaps);
    PPCGVectorOps<T, Device>::column_products(n_dim_, ld_psi_, ncol, spsi_.data(), spsi_.data(), overlap_norms);
    reduce_pool_values(moments.data(), moment_count);
    std::vector<double> steps(ncol, 0);
    const double ss = overlap_norms[selected];
    const double vs = response_overlaps[selected];
    const double rs = residual_overlaps[selected];
    if (!(ss > 0.0) || !std::isfinite(ss))
    {
        return refine_value();
    }
    // Eliminate the eigenvalue correction from
    // min ||r + step * (H - lambda S) d - delta_lambda * Sx||^2.
    const double denominator = response_norms[selected] - vs * vs / ss;
    if (!(denominator > 0.0) || !std::isfinite(denominator))
    {
        return refine_value();
    }
    steps[selected] = -(moments[selected] - rs * vs / ss) / denominator;
    if (!std::isfinite(steps[selected]))
    {
        return refine_value();
    }
    auto device_steps = upload_values<double, Device>(steps.data(), ncol);
    PPCGVectorOps<T, Device>::scale_columns(
        n_dim_, ld_psi_, ncol, device_steps.data(), direction.data(), s_direction.data(), h_direction.data());
    const T one(1);
    const int selected_offset = selected * ld_psi_;
    const T *selected_direction = direction.data() + selected_offset;
    T *selected_state = psi_in + selected_offset;
    ModuleBase::axpy_op<T, Device>()(ld_psi_, &one, selected_direction, 1, selected_state, 1);
    apply_h(op, psi_in, hpsi_.data(), ncol);
    apply_s(op, psi_in, spsi_.data(), ncol);
    const int normalization_count = 4 * ncol;
    std::vector<double> normalization(normalization_count, 0.0);
    double *state_norms = normalization.data() + ncol;
    double *candidate_residual_overlaps = normalization.data() + 2 * ncol;
    double *candidate_overlap_norms = normalization.data() + 3 * ncol;
    // Evaluate the Rayleigh correction from the small residual instead of
    // subtracting nearly equal Rayleigh quotients after division.
    PPCGVectorOps<T, Device>::residual(
        n_dim_, ld_psi_, ncol, hpsi_.data(), spsi_.data(), device_values.data(), response.data());
    PPCGVectorOps<T, Device>::column_products(n_dim_, ld_psi_, ncol, psi_in, response.data(), normalization.data());
    PPCGVectorOps<T, Device>::column_products(n_dim_, ld_psi_, ncol, psi_in, spsi_.data(), state_norms);
    PPCGVectorOps<T, Device>::column_products(
        n_dim_, ld_psi_, ncol, spsi_.data(), response.data(), candidate_residual_overlaps);
    PPCGVectorOps<T, Device>::column_products(
        n_dim_, ld_psi_, ncol, spsi_.data(), spsi_.data(), candidate_overlap_norms);
    reduce_pool_values(normalization.data(), normalization_count);
    std::vector<double> scales(ncol, 1);
    const double snorm = state_norms[selected];
    const double overlap_image_norm = candidate_overlap_norms[selected];
    Real candidate_energy = std::numeric_limits<Real>::quiet_NaN();
    if (snorm > 0 && overlap_image_norm > 0)
    {
        scales[selected] = 1 / std::sqrt(snorm);
        candidate_energy = original_value + Real(normalization[selected] / snorm);
        // Minimize the unscaled generalized residual for the candidate vector.
        // Keep its Rayleigh energy separately for the variational safeguard.
        eigenvalue_in[selected] = original_value + Real(candidate_residual_overlaps[selected] / overlap_image_norm);
    }
    auto device_scales = upload_values<double, Device>(scales.data(), ncol);
    PPCGVectorOps<T, Device>::scale_columns(
        n_dim_, ld_psi_, ncol, device_scales.data(), psi_in, spsi_.data(), hpsi_.data());
    apply_h(op, psi_in, hpsi_.data(), ncol);
    apply_s(op, psi_in, spsi_.data(), ncol);
    compute_eigenpair_residual(eigenvalue_in);
    std::vector<double> new_norms(ncol, 0);
    PPCGVectorOps<T, Device>::column_products(n_dim_, ld_psi_, ncol, w_.data(), w_.data(), new_norms.data());
    reduce_pool_values(new_norms.data(), ncol);
    bool ordered = true;
    if (selected > 0)
    {
        ordered = eigenvalue_in[selected] >= eigenvalue_in[selected - 1] - margin;
    }
    if (selected + 1 < ncol)
    {
        ordered = ordered && eigenvalue_in[selected] <= eigenvalue_in[selected + 1] + margin;
    }
    const bool accepted = ordered && snorm > 0.0 && std::isfinite(snorm) && std::isfinite(new_norms[selected]) &&
                          new_norms[selected] < old_norms[selected] &&
                          eigenvalue_in[selected] <= original_value + margin &&
                          candidate_energy <= original_value + margin;
    if (!accepted)
    {
        scatter_cols(psi_in, selected_cols, original);
        eigenvalue_in[selected] = original_value;
        apply_h(op, psi_in, hpsi_.data(), ncol);
        apply_s(op, psi_in, spsi_.data(), ncol);
        compute_eigenpair_residual(eigenvalue_in);
        return refine_value();
    }
    return accepted;
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::compute_eigenpair_residual(const Real *eigenvalue)
{
    const PPCGFunctionTimer timer("compute_eigenpair_residual");
    set_zero<T, Device>(w_);
    auto values = upload_values<Real, Device>(eigenvalue, n_band_);
    PPCGVectorOps<T, Device>::residual(n_dim_, ld_psi_, n_band_, hpsi_.data(), spsi_.data(), values.data(), w_.data());
}

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::orthonormalize(T *psi)
{
    const PPCGFunctionTimer timer("orthonormalize");
    gram(psi, spsi_.data(), n_band_, n_band_, rr_ssub_, n_band_);
    bool local_success = false;
    try
    {
        container::kernels::lapack_potrf<T, container::DEVICE_CPU>()('U', n_band_, rr_ssub_.data(), n_band_);
        container::kernels::lapack_trtri<T, container::DEVICE_CPU>()('U', 'N', n_band_, rr_ssub_.data(), n_band_);
        local_success = true;
    }
    catch (const std::runtime_error &)
    {
        // The caller performs a common Rayleigh-Ritz fallback on all ranks.
    }
    if (!all_pool_operations_succeeded(local_success))
    {
        return false;
    }
    for (int j = 0; j < n_band_; ++j)
    {
        for (int i = j + 1; i < n_band_; ++i)
        {
            rr_ssub_[i + j * n_band_] = T(0);
        }
    }
    // X^H SX = R^H R. Right multiplication by the upper triangular R^{-1}
    // keeps a previously locked prefix independent of the active columns.
    const int coefficient_count = n_band_ * n_band_;
    auto device_coeff = upload_values<T, Device>(rr_ssub_.data(), coefficient_count);
    const int size = ld_psi_ * n_band_;
    const T one = T(1);
    const T zero = T(0);
    const auto rotate = [&](T *values)
    {
        copy_values<T, Device>(values, size, w_.data());
        zero_values<T, Device>(values, size);
        ModuleBase::gemm_op<T, Device>()(
            'N', 'N', n_dim_, n_band_, n_band_, &one, w_.data(), ld_psi_, device_coeff.data(), n_band_, &zero, values,
            ld_psi_);
    };
    rotate(psi);
    rotate(hpsi_.data());
    rotate(spsi_.data());
    // The non-Ritz sweep changes the column basis too. Transform the
    // conjugate directions and both cached images with the same coefficients.
    rotate(p_.data());
    rotate(sp_.data());
    rotate(hp_.data());
    return true;
}

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::update_block_gradient(const T *psi, Real *eigenvalue, const std::vector<double> &ethr_band)
{
    const PPCGFunctionTimer timer("update_block_gradient");
    gram(psi, hpsi_.data(), n_band_, n_band_, rr_hsub_, n_band_);
    hermitize_projected(rr_hsub_, n_band_);
    for (int j = 0; j < n_band_; ++j)
    {
        // These are Rayleigh quotients for tracing, not ordered Ritz values.
        eigenvalue[j] = std::real(rr_hsub_[j + j * n_band_]);
    }
    copy_values<T, Device>(hpsi_.data(), int(hpsi_.size()), w_.data());
    const int coefficient_count = n_band_ * n_band_;
    auto device_coeff = upload_values<T, Device>(rr_hsub_.data(), coefficient_count);
    const T minus_one = T(-1);
    const T one = T(1);
    ModuleBase::gemm_op<T, Device>()(
        'N', 'N', n_dim_, n_band_, n_band_, &minus_one, spsi_.data(), ld_psi_, device_coeff.data(), n_band_, &one,
        w_.data(), ld_psi_);
    std::vector<int> candidates;
    lock_epairs(w_, ethr_band, candidates);
    // Do not change the actual locked prefix without a successful Ritz solve.
    return candidates.empty();
}

} // namespace hsolver

namespace hsolver
{

//==============================================================================
// MAIN DIAGONALIZATION ROUTINE
//==============================================================================
template <typename T, typename Device>
double DiagoPPCG<T, Device>::diag(
    const HSOperator<T, Device> &op,
    int ld_psi,
    int nband,
    int dim,
    T *psi_in,
    Real *eigenvalue_in,
    const std::vector<double> &ethr_band,
    const Real *prec,
    std::ostream &log,
    const PPCGTraceContext &trace_context)
{
    const PPCGFunctionTimer timer("diag");
    ld_psi_ = ld_psi;
    n_band_ = nband;
    n_dim_ = dim;
    converged_ = false;
    active_band_count_ = 0;
    iteration_limit_reached_ = false;
    max_residual_ = Real(0);
    first_active_band_ = -1;
    first_active_threshold_ = Real(0);
    rayleigh_ritz_count_ = 0;
    root_recovery_count_ = 0;

    validate_input(psi_in, eigenvalue_in, ethr_band, prec);

    device_prec_ = upload_values<Real, Device>(prec, n_dim_);

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
    rr_hsub_.resize(ncol * ncol);
    rr_ssub_.resize(ncol * ncol);
    rr_eval_.resize(ncol);

    std::vector<int> all_cols(ncol);
    std::iota(all_cols.begin(), all_cols.end(), 0);

    apply_h(op, psi_in, hpsi_.data(), ncol);
    apply_s(op, psi_in, spsi_.data(), ncol);

    double avg_iter = 1.0;
    int iter = 1;
    std::vector<int> active_cols;
    active_cols.reserve(ncol);

    std::ofstream residual_trace;
    std::streamoff trace_id = 0;
    const char *residual_trace_path = std::getenv("ABACUS_PPCG_RESIDUAL_TRACE");
    bool trace_writer = false;
    if (residual_trace_path != nullptr)
    {
        const auto world = ModuleBase::world_comm_domain();
        trace_writer = world.rank() == 0;
    }
    if (trace_writer)
    {
        // A single writer avoids file corruption when several MPI pools solve
        // different k-points concurrently.
        residual_trace.open(residual_trace_path, std::ios::out | std::ios::app | std::ios::ate);
        if (residual_trace.is_open())
        {
            trace_id = residual_trace.tellp();
        }
        if (residual_trace)
        {
#ifdef _WIN32
            char *full_path = ::_fullpath(nullptr, residual_trace_path, 0);
#else
            char *full_path = ::realpath(residual_trace_path, nullptr);
#endif
            if (full_path != nullptr)
            {
                const std::string trace_filename(full_path);
                std::free(full_path);
                log << "PPCG residual trace written: " << trace_filename << std::endl;
            }
            else
            {
                log << "WARNING: PPCG residual trace opened, but the full path could not be resolved: "
                    << residual_trace_path << std::endl;
            }
            residual_trace << std::setprecision(std::numeric_limits<Real>::max_digits10);
            // Preserve earlier k-points and SCF/ionic steps. The append offset
            // identifies this solve without a global counter or persistent state.
            if (trace_id == 0)
            {
                residual_trace << "iteration,stage,max_residual";
                for (int ib = 0; ib < ncol; ++ib)
                {
                    residual_trace << ",eigenvalue_" << ib;
                }
                residual_trace << ",local_overlap_min,local_overlap_max,local_overlap_condition"
                               << ",scf_iteration,kpoint,solve_id\n";
            }
        }
    }
    const bool write_residual_trace = residual_trace.is_open() && residual_trace.good();
    int residual_trace_writer_count = write_residual_trace ? 1 : 0;
    reduce_pool_values(&residual_trace_writer_count, 1);
    const bool residual_trace_enabled = residual_trace_writer_count > 0;
    Real local_overlap_min = std::numeric_limits<Real>::quiet_NaN();
    Real local_overlap_max = std::numeric_limits<Real>::quiet_NaN();
    Real local_overlap_condition = std::numeric_limits<Real>::quiet_NaN();
    auto record_residual = [&](int iteration, const char *stage)
    {
        if (!residual_trace_enabled)
        {
            return;
        }
        const Real max_residual =
            max_generalized_residual<T, Device>(hpsi_.data(), spsi_.data(), eigenvalue_in, ld_psi_, n_dim_, ncol);
        if (!write_residual_trace)
        {
            return;
        }
        residual_trace << iteration << ',' << stage << ',' << max_residual;
        for (int ib = 0; ib < ncol; ++ib)
        {
            residual_trace << ',' << eigenvalue_in[ib];
        }
        residual_trace << ',' << local_overlap_min << ',' << local_overlap_max << ',' << local_overlap_condition << ','
                       << trace_context.scf_iteration << ',' << trace_context.kpoint << ',' << trace_id << '\n';
    };

    auto measure_local_overlap = [&](const SmallSubspace &local_subspace, const int local_dim)
    {
        if (!write_residual_trace)
        {
            return;
        }
        std::vector<T> overlap = local_subspace.m;
        std::vector<Real> overlap_eigenvalues(local_dim);
        try
        {
            container::kernels::lapack_heevd<T, container::DEVICE_CPU>()(
                local_dim, overlap.data(), local_dim, overlap_eigenvalues.data());
        }
        catch (const std::exception &)
        {
            // Trace diagnostics must not change the solver's MPI control flow.
            local_overlap_min = std::numeric_limits<Real>::quiet_NaN();
            local_overlap_max = std::numeric_limits<Real>::quiet_NaN();
            local_overlap_condition = std::numeric_limits<Real>::infinity();
            return;
        }
        const Real min_eigenvalue = overlap_eigenvalues.front();
        const Real max_eigenvalue = overlap_eigenvalues.back();
        const Real minimum_overlap_magnitude = Real(1.0e-30);
        const Real denominator = std::max(std::abs(min_eigenvalue), minimum_overlap_magnitude);
        const Real condition = std::abs(max_eigenvalue) / denominator;
        if (std::isnan(local_overlap_condition) || condition > local_overlap_condition)
        {
            local_overlap_min = min_eigenvalue;
            local_overlap_max = max_eigenvalue;
            local_overlap_condition = condition;
        }
    };

    // Initialize with Rayleigh-Ritz.
    bool last_rr_succeeded = rayleigh_ritz(psi_in, eigenvalue_in, active_cols, ethr_band);
    const bool needs_startup_global_correction = !active_cols.empty() && sbsize_ < ncol;
    // The Ritz transformation updates X, HX and SX with identical coefficients.
    record_residual(0, "initial_rr");

    Workspace active_vector;
    Workspace s_active_vector;
    Workspace h_active_vector;
    std::vector<int> cols;
    cols.reserve(std::min(sbsize_, ncol));
    SmallSubspace subspace;
    bool use_p = false;                               // previous search direction becomes available after
                                                      // the first block update.
    Real best_res = std::numeric_limits<Real>::max(); // best residual seen
    int no_improve = 0;                               // iterations since the residual last decreased
    int final_refinement_band = -1;
    int next_final_refinement_band = 0;

    while (!active_cols.empty() && iter <= maxiter_)
    {
        OperatorBlock discarded;
        // Near the floating-point cancellation scale, form a reliable
        // residual before using it as a new search direction. This safeguard
        // also covers cancellation in the W/P recurrence, not only in the
        // explicit projection of P. The final fresh residual remains the
        // convergence certificate; this scale is not an error bound.
        Real spectral_scale = Real(1);
        for (int ib = 0; ib < ncol; ++ib)
        {
            spectral_scale = std::max(spectral_scale, std::abs(eigenvalue_in[ib]));
        }
        const Real reliable_update_threshold = std::sqrt(std::numeric_limits<Real>::epsilon()) * spectral_scale;
        const bool near_precision = max_residual_ <= reliable_update_threshold;
        if (near_precision)
        {
            apply_h(op, psi_in, hpsi_.data(), ncol);
            apply_s(op, psi_in, spsi_.data(), ncol);
            apply_h(op, p_.data(), hp_.data(), ncol);
            apply_s(op, p_.data(), sp_.data(), ncol);
            if (last_rr_succeeded)
            {
                compute_eigenpair_residual(eigenvalue_in);
            }
            else
            {
                update_block_gradient(psi_in, eigenvalue_in, ethr_band);
            }
        }
        // At stagnation within a few stopping tolerances, another global
        // Ritz rotation can perturb already accurate vectors. Try a single
        // normalized residual correction before resuming the block iteration.
        // Do not wait past the remaining budget: reserve the last few steps
        // for the worst residual rather than allowing another Ritz rotation
        // to perturb near-converged states before refinement can run.
        const bool startup_correction = iter == 1 && needs_startup_global_correction;
        const Real refinement_threshold = Real(4) * diag_thr_;
        const int refinement_stagnation = 5;
        // One correction updates one band. Leave room for a full retained
        // set as well as the stagnation delay within the existing budget.
        const int refinement_budget = ncol + refinement_stagnation;
        const bool final_refinement_window = maxiter_ - iter < refinement_budget;
        const bool refinement_available = !final_refinement_window || final_refinement_band < ncol;
        if (!startup_correction && last_rr_succeeded && near_precision && max_residual_ <= refinement_threshold &&
            refinement_available && (no_improve >= refinement_stagnation || final_refinement_window))
        {
            // Cycle through the retained states: errors in an already locked
            // neighbor can limit the attainable residual of the remaining pair.
            const int refinement_band = final_refinement_window ? final_refinement_band : iter % ncol;
            const bool accepted =
                refine_stalled_residual(op, psi_in, eigenvalue_in, ethr_band, all_cols, refinement_band);
            if (final_refinement_window)
            {
                // A rejected worst-pair step cannot improve if retried with
                // identical inputs. Visit its retained neighbors, then allow
                // a block sweep if the entire cycle fails to make progress.
                if (accepted)
                {
                    final_refinement_band = -1;
                }
                else
                {
                    final_refinement_band = next_final_refinement_band;
                    ++next_final_refinement_band;
                }
            }
            if (!accepted && refinement_band == 0)
            {
                no_improve = 0;
            }
            lock_epairs(w_, ethr_band, active_cols);
            record_residual(iter, "residual_refinement");
            set_zero<T, Device>(p_);
            set_zero<T, Device>(sp_);
            set_zero<T, Device>(hp_);
            use_p = false;
            avg_iter += 1.0 / double(ncol);
            ++iter;
            continue;
        }
        final_refinement_band = -1;
        next_final_refinement_band = 0;
        local_overlap_min = std::numeric_limits<Real>::quiet_NaN();
        local_overlap_max = std::numeric_limits<Real>::quiet_NaN();
        local_overlap_condition = std::numeric_limits<Real>::quiet_NaN();
        const int nact = int(active_cols.size());
        const int nsb = 1 + (nact - 1) / sbsize_;
        // Keep the same number of local blocks and the same size limit, but
        // distribute the remainder evenly.  When the limit permits it, this
        // avoids a trailing one-column block that cannot exchange search
        // directions with a neighboring degenerate Ritz vector.
        const int base_block_size = nact / nsb;
        const int larger_block_count = nact % nsb;

        // Precondition the residual.
        divide_by_preconditioner(active_cols, device_prec_.data(), w_);
        copy_cols(w_.data(), active_cols, active_vector);
        s_active_vector.assign(ld_psi_ * nact, T(0));
        apply_s(op, active_vector.data(), s_active_vector.data(), nact);
        scatter_cols(sw_.data(), active_cols, s_active_vector);
        project_against(psi_in, spsi_.data(), nullptr, all_cols, w_, sw_, nullptr, active_cols);

        // Apply H to the search direction.
        copy_cols(w_.data(), active_cols, active_vector);
        h_active_vector.assign(ld_psi_ * nact, T(0));
        copy_cols(sw_.data(), active_cols, s_active_vector);
        apply_h(op, active_vector.data(), h_active_vector.data(), nact);
        scatter_cols(hw_.data(), active_cols, h_active_vector);

        // Project the cached P, SP and HP together. Linearity preserves their
        // consistency without another Hamiltonian/overlap application except
        // when strong cancellation requires a fresh image of selected columns.
        if (use_p)
        {
            project_search_directions(op, psi_in, all_cols, active_cols);
        }

        avg_iter += double(nact) / double(ncol);

        // LOBPCG-style block subspace.  On the first sweep only [psi, w] is
        // available; afterwards the previous search direction p is added as a
        // third block, which restores the conjugate-gradient acceleration.
        // The w/p blocks are normalized to unit S-norm before building the
        // Gram matrix (see build_small_subspace), keeping M well-conditioned.

        // The startup correction needs a global [X, W_active] comparison to
        // retain the requested lowest roots.  All current Ritz vectors remain
        // available, while a converged prefix contributes no correction
        // vectors and stays outside the active iteration count.
        if (startup_correction)
        {
            const int local_dim = ncol + nact;
            build_startup_global_subspace(psi_in, active_vector, s_active_vector, h_active_vector, nact, subspace);
            measure_local_overlap(subspace, local_dim);
            solve_small_generalized(local_dim, ncol, subspace);
            update_startup_global_subspace(psi_in, active_vector, s_active_vector, h_active_vector, nact, subspace);
        }
        else
        {
            const int nblk = use_p ? 3 : 2;
            int i0 = 0;
            for (int isb = 0; isb < nsb; ++isb)
            {
                int l = base_block_size;
                if (isb < larger_block_count)
                {
                    ++l;
                }
                cols.assign(active_cols.begin() + i0, active_cols.begin() + i0 + l);

                build_small_subspace(psi_in, cols, nblk, subspace);
                const int local_dim = nblk * l;
                measure_local_overlap(subspace, local_dim);
                const int valid_ritz_count = solve_small_generalized(local_dim, l, subspace);
                collect_discarded_ritz(subspace, l, valid_ritz_count, discarded);
                update_one_block(psi_in, cols, l, nblk, subspace);
                i0 += l;
            }
        }
        use_p = true;

        // Refresh cached images independently of the full Ritz cadence, so
        // rr_step=1 does not add two full H/S applications on every sweep.
        // Intermediate sweeps maintain S-orthonormality without diagonalizing
        // the full projected Hamiltonian.
        const bool scheduled_rr = (iter % rr_step_) == 0;
        const int operator_refresh_interval = 16;
        const bool refresh_operators = (iter % operator_refresh_interval) == 0 && !near_precision;
        if (refresh_operators)
        {
            apply_h(op, psi_in, hpsi_.data(), ncol);
            apply_s(op, psi_in, spsi_.data(), ncol);
            apply_h(op, p_.data(), hp_.data(), ncol);
            apply_s(op, p_.data(), sp_.data(), ncol);
        }
        // The startup subspace can insert a newly discovered lower root into
        // an old locked prefix. Reclassify it before freezing any column.
        bool perform_rr = scheduled_rr || startup_correction || iter == maxiter_;
        if (!perform_rr)
        {
            perform_rr = !orthonormalize(psi_in);
        }
        if (!perform_rr)
        {
            perform_rr = update_block_gradient(psi_in, eigenvalue_in, ethr_band);
        }
        last_rr_succeeded = false;
        if (perform_rr)
        {
            last_rr_succeeded = rayleigh_ritz(psi_in, eigenvalue_in, active_cols, ethr_band);
        }
        // A residual only certifies the retained eigenpairs, not that no lower
        // state was discarded by another local block. Inspect those available
        // directions before accepting a newly converged set of columns.
        if (!perform_rr || last_rr_succeeded)
        {
            const bool recovered =
                recover_competing_ritz(op, psi_in, eigenvalue_in, last_rr_succeeded, active_cols, ethr_band, discarded);
            if (recovered)
            {
                last_rr_succeeded = true;
                perform_rr = true;
                use_p = false;
            }
        }
        if (last_rr_succeeded && active_cols.empty())
        {
            // A small cached residual is only tentative. Refresh the
            // operator images before stopping, and keep iterating if
            // roundoff had made the cached criterion overly optimistic.
            apply_h(op, psi_in, hpsi_.data(), ncol);
            apply_s(op, psi_in, spsi_.data(), ncol);
            compute_eigenpair_residual(eigenvalue_in);
            lock_epairs(w_, ethr_band, active_cols);
            if (!active_cols.empty())
            {
                set_zero<T, Device>(p_);
                set_zero<T, Device>(sp_);
                set_zero<T, Device>(hp_);
                use_p = false;
            }
        }
        // Restart the search direction if the residual keeps rising.  With a
        // poor preconditioner the LOBPCG recurrence can stagnate (or slowly
        // diverge) instead of reducing the residual.  Requiring several
        // consecutive rises avoids resetting on a transient bump, and a reset
        // falls back to a steepest-descent step to recover the low eigenpairs.
        {
            const int stagnation_restart_interval = 15;
            const Real cur_res = max_residual_;
            if (cur_res < best_res)
            {
                best_res = cur_res;
                no_improve = 0;
            }
            else
            {
                ++no_improve;
            }
            if (no_improve >= stagnation_restart_interval)
            {
                set_zero<T, Device>(p_);
                set_zero<T, Device>(sp_);
                set_zero<T, Device>(hp_);
                use_p = false;
                no_improve = 0;
                best_res = cur_res;
            }
        }
        const char *stage = perform_rr ? "rayleigh_ritz" : "block_sweep";
        record_residual(iter, stage);

        ++iter;
    }

    // Final consistency: ensure H|psi> and S|psi> match the returned vectors,
    // then classify convergence from those final residuals rather than from a
    // cached Rayleigh-Ritz state.
    apply_h(op, psi_in, hpsi_.data(), ncol);
    apply_s(op, psi_in, spsi_.data(), ncol);
    compute_eigenpair_residual(eigenvalue_in);
    if (last_rr_succeeded)
    {
        lock_epairs(w_, ethr_band, active_cols);
    }
    else
    {
        active_cols = all_cols;
        max_residual_ =
            max_generalized_residual<T, Device>(hpsi_.data(), spsi_.data(), eigenvalue_in, ld_psi_, n_dim_, ncol);
    }
    active_band_count_ = int(active_cols.size());
    if (!active_cols.empty())
    {
        first_active_band_ = active_cols.front();
        first_active_threshold_ = std::max(Real(ethr_band[first_active_band_]), diag_thr_);
    }
    converged_ = last_rr_succeeded && active_cols.empty();
    iteration_limit_reached_ = !converged_ && iter > maxiter_;
    if (!active_cols.empty())
    {
        const char *failure_stage = iteration_limit_reached_ ? "max_iterations" : "final_residual_failure";
        record_residual(iter - 1, failure_stage);
    }
    record_residual(iter - 1, "final");
    return avg_iter;
}

} // namespace hsolver

namespace hsolver
{

template class DiagoPPCG<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_CPU>;
#ifdef __CUDA
template class DiagoPPCG<std::complex<float>, base_device::DEVICE_GPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hsolver
