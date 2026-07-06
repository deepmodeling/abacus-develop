#include "source_base/kernels/math_kernel_op.h"
namespace hsolver {
namespace {

inline bool ppcg_contiguous_cols(const std::vector<int>& cols, int& first)
{
    if (cols.empty())
        return false;

    first = cols.front();
    for (int j = 0; j < static_cast<int>(cols.size()); ++j)
    {
        if (cols[j] != first + j)
            return false;
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
                                 const bool gamma_g0_real,
                                 const PpcgStrategy strategy)
    : maxiter_(diag_iter_max),
      sbsize_(std::max(1, sbsize)),
      rr_step_(std::max(1, rr_step)),
      diag_thr_(std::max(diag_thr, static_cast<Real>(1.0e-14))),
      gamma_g0_real_(gamma_g0_real),
      strategy_(strategy)
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
        throw std::invalid_argument("PPCG: H operator is empty.");
    if (psi_in == nullptr || eigenvalue_in == nullptr)
        throw std::invalid_argument("PPCG: psi/eigenvalue pointer is null.");
    if (prec == nullptr)
        throw std::invalid_argument("PPCG: preconditioner pointer is null.");
    if (ld_psi_ <= 0 || n_band_ <= 0 || n_dim_ <= 0)
        throw std::invalid_argument("PPCG: invalid dimensions.");
    if (n_dim_ > ld_psi_)
        throw std::invalid_argument("PPCG: dim must not exceed ld_psi.");
    if (ethr_band.size() < static_cast<size_t>(n_band_))
        throw std::invalid_argument("PPCG: ethr_band size is smaller than nband.");
    for (int i = 0; i < n_band_; ++i)
        if (!std::isfinite(ethr_band[i]))
            throw std::invalid_argument("PPCG: ethr_band contains non-finite value.");
    for (int i = 0; i < n_dim_; ++i)
        if (!std::isfinite(prec[i]))
            throw std::invalid_argument("PPCG: preconditioner contains non-finite value.");
}

// =============================================================================
// Gamma-point symmetry: enforce real-valued first element
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::force_g0_real(T* x, int ncol) const
{
    if (!gamma_g0_real_ || n_dim_ <= 0)
        return;
    for (int j = 0; j < ncol; ++j)
        x[idx(0, j, ld_psi_)] = T(std::real(x[idx(0, j, ld_psi_)]), 0.0);
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
        spsi_func(psi_in, spsi_out, ld_psi_, ncol);
    else
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ld_psi_ * ncol > 4096)
#endif
        for (int j = 0; j < ncol; ++j)
            std::copy(psi_in + j * ld_psi_, psi_in + (j + 1) * ld_psi_,
                      spsi_out + j * ld_psi_);
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
        acc += std::conj(x[i]) * y[i];
    reduce_pool_if_mpi_ready(&acc, 1);
    return acc;
}

// =============================================================================
// Gram matrix: out[i, j] = <a_i | b_j>
// =============================================================================
template <typename T, typename Device>
void DiagoPPCG<T, Device>::gram(const T* a, const T* b,
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
                                     a,
                                     ld_psi_,
                                     b,
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
    const int ncols = static_cast<int>(cols.size());
    dst.resize(ld_psi_ * ncols);
    if (ncols == 0)
        return;

    int first = 0;
    if (ppcg_contiguous_cols(cols, first))
    {
        std::copy(src + first * ld_psi_,
                  src + (first + ncols) * ld_psi_,
                  dst.begin());
        return;
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ld_psi_ * ncols > 4096)
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
    const int ncols = static_cast<int>(cols.size());
    if (ncols == 0)
        return;

    int first = 0;
    if (ppcg_contiguous_cols(cols, first))
    {
        std::copy(src.begin(),
                  src.begin() + ld_psi_ * ncols,
                  dst + first * ld_psi_);
        return;
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ld_psi_ * ncols > 4096)
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
        return;

    const int nbasis = static_cast<int>(basis_cols.size());
    const int nx = static_cast<int>(x_cols.size());

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
    const int ncols = static_cast<int>(active_cols.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ * ncols > 4096)
#endif
    for (int j = 0; j < ncols; ++j)
    {
        const int c = active_cols[j];
        for (int ig = 0; ig < n_dim_; ++ig)
            x[idx(ig, c, ld_psi_)] /=
                std::max(prec[ig], static_cast<Real>(1.0e-12));
    }
}

} // namespace hsolver
