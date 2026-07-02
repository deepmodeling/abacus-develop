#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_reduce.h"

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
        Parallel_Reduce::reduce_pool(value);
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
        Parallel_Reduce::reduce_pool(value, n);
#endif
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
    const T* psi_in,
    const Real* eigenvalue_in,
    const std::vector<double>& ethr_band,
    const Real* prec) const
{
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
    out.assign(ld_out * ncol_b, T(0));
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
    dst.assign(ld_psi_ * cols.size(), T(0));
    for (int j = 0; j < static_cast<int>(cols.size()); ++j)
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
    for (int j = 0; j < static_cast<int>(cols.size()); ++j)
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

    for (const int c : x_cols)
    {
        for (const int bc : basis_cols)
        {
            // Full complex inner product <basis_bc | sx_c>
            T coeff = 0;
            const T* bb = basis + bc * ld_psi_;
            const T* sc = sx.data() + c * ld_psi_;
            for (int ig = 0; ig < n_dim_; ++ig)
                coeff += std::conj(bb[ig]) * sc[ig];
            if (std::abs(coeff) <= std::numeric_limits<Real>::epsilon())
                continue;
            const T* sb = sbasis + bc * ld_psi_;
            T* xc = x.data() + c * ld_psi_;
            T* sxc = sx.data() + c * ld_psi_;
            for (int ig = 0; ig < n_dim_; ++ig)
            {
                xc[ig] -= bb[ig] * coeff;
                sxc[ig] -= sb[ig] * coeff;
            }
        }
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
    for (const int c : active_cols)
        for (int ig = 0; ig < n_dim_; ++ig)
            x[idx(ig, c, ld_psi_)] /=
                std::max(prec[ig], static_cast<Real>(1.0e-12));
}

} // namespace hsolver
