#include "diago_ppcg.h"

#include <cstdint>

#include "ATen/kernels/lapack.h"

#include <cstdlib>
#include <fstream>

namespace hsolver {

// =============================================================================
// Small dense eigensolver helpers
// =============================================================================
namespace {

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
    for (int j = 0; j < ncol; ++j)
    {
        Real nrm2 = 0;
        for (int ig = 0; ig < n_dim; ++ig)
        {
            const T r = hpsi[ig + j * ld] - T(eigenvalue[j]) * spsi[ig + j * ld];
            nrm2 += static_cast<Real>(std::norm(r));
        }
        max_res = std::max(max_res, std::sqrt(nrm2));
    }
    return max_res;
}

template <typename T>
struct PpcgLapack
{
    using Real = typename ct::kernels::lapack_hegvx<T, base_device::DEVICE_CPU>::Real;

    static void heevd(int n, T* a, Real* w)
    {
        ct::kernels::lapack_heevd<T, base_device::DEVICE_CPU>()(n, a, n, w);
    }

    static void hegvd(int n, T* a, T* b, Real* w)
    {
        std::vector<T> eigen_vec(n * n, T(0));
        ct::kernels::lapack_hegvx<T, base_device::DEVICE_CPU>()(
            n, n, a, b, n, w, eigen_vec.data());
        std::copy(eigen_vec.begin(), eigen_vec.end(), a);
    }

    static void potrf(int n, T* a)
    {
        const char uplo = 'U';
        const int lda = n;
        Real diag_max = 0;
        for (int i = 0; i < n; ++i)
            diag_max = std::max(diag_max, static_cast<Real>(std::abs(a[i + i * lda])));
        const std::vector<T> a0(a, a + n * lda);

        const Real shifts[] = {static_cast<Real>(0),
                               static_cast<Real>(1e-12),
                               static_cast<Real>(1e-10),
                               static_cast<Real>(1e-8),
                               static_cast<Real>(1e-6),
                               static_cast<Real>(1e-4),
                               static_cast<Real>(1e-3),
                               static_cast<Real>(1e-2),
                               static_cast<Real>(1e-1),
                               static_cast<Real>(1)};
        for (const Real shift : shifts)
        {
            std::copy(a0.begin(), a0.end(), a);
            if (shift > 0)
            {
                const Real scaled_shift = shift * std::max(diag_max, static_cast<Real>(1));
                for (int i = 0; i < n; ++i)
                    a[i + i * lda] += T(scaled_shift);
            }
            try
            {
                ct::kernels::lapack_potrf<T, base_device::DEVICE_CPU>()(
                    uplo, n, a, lda);
                return;
            }
            catch (const std::runtime_error&)
            {
                // Try the next diagonal shift.
            }
        }
        throw std::runtime_error("PPCG: lapack_potrf failed.");
    }

    static void trtri(int n, T* a)
    {
        const char uplo = 'U';
        const char diag = 'N';
        const int lda = n;
        ct::kernels::lapack_trtri<T, base_device::DEVICE_CPU>()(
            uplo, diag, n, a, lda);
    }
};

template <typename T>
inline void set_zero(std::vector<T>& x)
{
    std::fill(x.begin(), x.end(), T(0));
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
    Real acc = 0;
    for (int i = 0; i < n_dim_; ++i)
        acc += static_cast<Real>(std::real(std::conj(x[i]) * y[i]));
    return acc;
}

template <typename T, typename Device>
T DiagoPPCG<T, Device>::complex_dot(const T* x, const T* y) const
{
    T acc = T(0);
    for (int i = 0; i < n_dim_; ++i)
        acc += std::conj(x[i]) * y[i];
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
    for (int jb = 0; jb < ncol_b; ++jb)
        for (int ia = 0; ia < ncol_a; ++ia)
            out[ia + jb * ld_out] = complex_dot(a + ia * ld_psi_,
                                                 b + jb * ld_psi_);
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

//==============================================================================
// BLOCK_SUBSPACE STRATEGY
//==============================================================================

// ---------------------------------------------------------------------------
// Lock converged eigenpairs: columns with residual below threshold
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::lock_epairs(
    const std::vector<T>& residual,
    const std::vector<double>& ethr_band,
    std::vector<int>& active_cols) const
{
    active_cols.clear();
    for (int j = 0; j < n_band_; ++j)
    {
        Real nrm2 = 0;
        for (int ig = 0; ig < n_dim_; ++ig)
            nrm2 += static_cast<Real>(std::norm(residual[idx(ig, j, ld_psi_)]));
        const Real rnrm = std::sqrt(std::max(nrm2, static_cast<Real>(0)));
        const Real thr = std::max(static_cast<Real>(ethr_band[j]), diag_thr_);
        if (rnrm > thr)
            active_cols.push_back(j);
    }
}

// ---------------------------------------------------------------------------
// Build K = V^H H V and M = V^H S V where V = [psi, w, p]
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::build_small_subspace(
    const T* psi,
    const std::vector<int>& cols,
    bool use_p,
    SmallSubspace& subspace) const
{
    const int l = static_cast<int>(cols.size());
    const int nblk = use_p ? 3 : 2;
    const int dim = nblk * l;
    subspace.k.assign(dim * dim, T(0));
    subspace.m.assign(dim * dim, T(0));
    subspace.eval.assign(dim, static_cast<Real>(0));

    std::vector<T> psi_l, spsi_l, hpsi_l;
    std::vector<T> w_l, sw_l, hw_l;
    std::vector<T> p_l, sp_l, hp_l;
    copy_cols(psi, cols, psi_l);
    copy_cols(spsi_.data(), cols, spsi_l);
    copy_cols(hpsi_.data(), cols, hpsi_l);
    copy_cols(w_.data(), cols, w_l);
    copy_cols(sw_.data(), cols, sw_l);
    copy_cols(hw_.data(), cols, hw_l);
    if (use_p)
    {
        copy_cols(p_.data(), cols, p_l);
        copy_cols(sp_.data(), cols, sp_l);
        copy_cols(hp_.data(), cols, hp_l);
    }

    // ---------------------------------------------------------------------------
    // Normalize w and p columns to unit S-norm for numerical stability.
    //
    // The [w, p] block of the Gram matrix M has entries O(||w||^2) which
    // become tiny when residuals are small, making M nearly singular and
    // causing sygvd to produce garbage eigenvectors.
    //
    // Scaling to unit S-norm keeps M well-conditioned (diagonal ~1) without
    // changing the subspace.  The Ritz values are identical and the Ritz
    // vector coefficients in update_one_block automatically compensate.
    // ---------------------------------------------------------------------------
    auto scale_to_unit_snorm = [this](std::vector<T>& x, std::vector<T>& sx,
                                       std::vector<T>& hx, int lcols) {
        for (int j = 0; j < lcols; ++j) {
            Real sn2 = 0;
            for (int ig = 0; ig < n_dim_; ++ig)
                sn2 += std::real(std::conj(x[idx(ig, j, ld_psi_)])
                                 * sx[idx(ig, j, ld_psi_)]);
            Real sn = std::sqrt(std::max(sn2, static_cast<Real>(1e-30)));
            // Only scale if the norm is non-negligible; a near-zero
            // column is a converged band whose contribution is harmless.
            if (sn > static_cast<Real>(1e-15)) {
                Real inv = static_cast<Real>(1) / sn;
                for (int ig = 0; ig < n_dim_; ++ig) {
                    x[ idx(ig, j, ld_psi_)]  *= inv;
                    sx[idx(ig, j, ld_psi_)] *= inv;
                    hx[idx(ig, j, ld_psi_)] *= inv;
                }
            }
        }
    };
    scale_to_unit_snorm(w_l, sw_l, hw_l, l);
    if (use_p)
        scale_to_unit_snorm(p_l, sp_l, hp_l, l);

    auto fill_sym = [&](const std::vector<T>& a, const std::vector<T>& b,
                        int r0, int c0, std::vector<T>& mat)
    {
        std::vector<T> g;
        gram(a.data(), b.data(), l, l, g, l);
        for (int j = 0; j < l; ++j)
            for (int i = 0; i < l; ++i)
            {
                mat[(r0 + i) + (c0 + j) * dim] = g[i + j * l];
                mat[(c0 + j) + (r0 + i) * dim] = std::conj(g[i + j * l]);
            }
    };

    fill_sym(psi_l, hpsi_l, 0,   0,   subspace.k);
    fill_sym(psi_l, spsi_l, 0,   0,   subspace.m);
    fill_sym(w_l,   hw_l,   l,   l,   subspace.k);
    fill_sym(w_l,   sw_l,   l,   l,   subspace.m);
    fill_sym(psi_l, hw_l,   0,   l,   subspace.k);
    fill_sym(psi_l, sw_l,   0,   l,   subspace.m);

    if (use_p)
    {
        fill_sym(p_l, hp_l, 2*l, 2*l, subspace.k);
        fill_sym(p_l, sp_l, 2*l, 2*l, subspace.m);
        fill_sym(psi_l, hp_l, 0,   2*l, subspace.k);
        fill_sym(psi_l, sp_l, 0,   2*l, subspace.m);
        fill_sym(w_l,   hp_l, l,   2*l, subspace.k);
        fill_sym(w_l,   sp_l, l,   2*l, subspace.m);
    }
}

// ---------------------------------------------------------------------------
// Solve K v = lambda M v (small generalized eigenvalue problem)
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
    const Real shifts[] = {static_cast<Real>(0),
                           static_cast<Real>(1e-10),
                           static_cast<Real>(1e-8),
                           static_cast<Real>(1e-6)};
    for (const Real shift : shifts)
    {
        subspace.k = k0;
        subspace.m = m0;
        for (int i = 0; i < dim; ++i)
            subspace.m[i + i * dim] += T(shift);

        try
        {
            PpcgLapack<T>::hegvd(dim, subspace.k.data(),
                                      subspace.m.data(),
                                      subspace.eval.data());
            return;
        }
        catch (const std::runtime_error&)
        {
            // Try the next diagonal shift.
        }
    }
    // All attempts failed; set eigenvectors to identity (no update).
    std::fill(subspace.k.begin(), subspace.k.end(), T(0));
    for (int i = 0; i < dim; ++i)
        subspace.k[i + i * dim] = T(1);
    std::fill(subspace.eval.begin(), subspace.eval.end(), static_cast<Real>(0));
}

// ---------------------------------------------------------------------------
// Update wavefunctions from small subspace eigenvectors
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::update_one_block(
    T* psi,
    const std::vector<int>& cols,
    int l,
    bool use_p,
    const SmallSubspace& subspace)
{
    const int dim = (use_p ? 3 : 2) * l;
    const T* eigvec = subspace.k.data();

    std::vector<T> psi_l, spsi_l, hpsi_l;
    std::vector<T> w_l, sw_l, hw_l;
    std::vector<T> p_l, sp_l, hp_l;
    copy_cols(psi, cols, psi_l);
    copy_cols(spsi_.data(), cols, spsi_l);
    copy_cols(hpsi_.data(), cols, hpsi_l);
    copy_cols(w_.data(), cols, w_l);
    copy_cols(sw_.data(), cols, sw_l);
    copy_cols(hw_.data(), cols, hw_l);
    if (use_p)
    {
        copy_cols(p_.data(), cols, p_l);
        copy_cols(sp_.data(), cols, sp_l);
        copy_cols(hp_.data(), cols, hp_l);
    }

    std::vector<T> psi_new(ld_psi_ * l, T(0));
    std::vector<T> spsi_new(ld_psi_ * l, T(0));
    std::vector<T> hpsi_new(ld_psi_ * l, T(0));
    std::vector<T> p_new(ld_psi_ * l, T(0));
    std::vector<T> sp_new(ld_psi_ * l, T(0));
    std::vector<T> hp_new(ld_psi_ * l, T(0));

    for (int j = 0; j < l; ++j)
    {
        for (int i = 0; i < l; ++i)
        {
            const T cpsi = eigvec[i + j * dim];
            const T cw   = eigvec[(l + i) + j * dim];

            for (int ig = 0; ig < n_dim_; ++ig)
            {
                psi_new[idx(ig, j, ld_psi_)]  += psi_l[idx(ig, i, ld_psi_)] * cpsi
                                                + w_l[ idx(ig, i, ld_psi_)] * cw;
                spsi_new[idx(ig, j, ld_psi_)] += spsi_l[idx(ig, i, ld_psi_)] * cpsi
                                                + sw_l[ idx(ig, i, ld_psi_)] * cw;
                hpsi_new[idx(ig, j, ld_psi_)] += hpsi_l[idx(ig, i, ld_psi_)] * cpsi
                                                + hw_l[ idx(ig, i, ld_psi_)] * cw;
                p_new[idx(ig, j, ld_psi_)]    += w_l[ idx(ig, i, ld_psi_)] * cw;
                sp_new[idx(ig, j, ld_psi_)]   += sw_l[ idx(ig, i, ld_psi_)] * cw;
                hp_new[idx(ig, j, ld_psi_)]   += hw_l[ idx(ig, i, ld_psi_)] * cw;
            }

            if (use_p)
            {
                const T cp = eigvec[(2*l + i) + j * dim];
                for (int ig = 0; ig < n_dim_; ++ig)
                {
                    psi_new[idx(ig, j, ld_psi_)]  += p_l[ idx(ig, i, ld_psi_)] * cp;
                    spsi_new[idx(ig, j, ld_psi_)] += sp_l[idx(ig, i, ld_psi_)] * cp;
                    hpsi_new[idx(ig, j, ld_psi_)] += hp_l[idx(ig, i, ld_psi_)] * cp;
                    p_new[idx(ig, j, ld_psi_)]    += p_l[ idx(ig, i, ld_psi_)] * cp;
                    sp_new[idx(ig, j, ld_psi_)]   += sp_l[idx(ig, i, ld_psi_)] * cp;
                    hp_new[idx(ig, j, ld_psi_)]   += hp_l[idx(ig, i, ld_psi_)] * cp;
                }
            }
        }
    }

    scatter_cols(psi, cols, psi_new);
    scatter_cols(spsi_.data(), cols, spsi_new);
    scatter_cols(hpsi_.data(), cols, hpsi_new);
    scatter_cols(p_.data(), cols, p_new);
    scatter_cols(sp_.data(), cols, sp_new);
    scatter_cols(hp_.data(), cols, hp_new);
}

// ---------------------------------------------------------------------------
// Back-substitute with upper triangular Cholesky factor: X *= R^{-1}
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::right_solve_upper(
    const std::vector<T>& r, int n, std::vector<T>& x) const
{
    std::vector<T> b = x;
    for (int row = 0; row < n_dim_; ++row)
    {
        for (int j = 0; j < n; ++j)
        {
            T v = b[idx(row, j, ld_psi_)];
            for (int k = 0; k < j; ++k)
                v -= x[idx(row, k, ld_psi_)] * r[k + j * n];
            x[idx(row, j, ld_psi_)] = v / r[j + j * n];
        }
    }
}

// ---------------------------------------------------------------------------
// Check S-orthonormality of a column block.
// ---------------------------------------------------------------------------
template <typename T, typename Device>
bool DiagoPPCG<T, Device>::is_s_orthonormal(
    const T* psi, const T* spsi, int ncol) const
{
    const Real orth_tol = static_cast<Real>(10)
                        * std::sqrt(std::numeric_limits<Real>::epsilon());
    for (int j = 0; j < ncol; ++j)
    {
        for (int i = 0; i < ncol; ++i)
        {
            const T sij = complex_dot(psi + i * ld_psi_,
                                      spsi + j * ld_psi_);
            const T target = (i == j) ? T(1) : T(0);
            if (std::abs(sij - target) > orth_tol)
                return false;
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
            static_cast<Real>(1e-30)));
        Real inv_nrm = static_cast<Real>(1) / nrm;
        for (int ig = 0; ig < n_dim_; ++ig)
        {
            psi [idx(ig, j, ld_psi_)] *= inv_nrm;
            hpsi[idx(ig, j, ld_psi_)] *= inv_nrm;
            spsi[idx(ig, j, ld_psi_)] *= inv_nrm;
        }
    }
}

// ---------------------------------------------------------------------------
// Cholesky QR: S-orthonormalize active columns via Cholesky on S-gram
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::chol_qr_active(
    T* psi, const std::vector<int>& active_cols)
{
    if (active_cols.empty())
        return;

    const int nact = static_cast<int>(active_cols.size());
    std::vector<T> psi_a, spsi_a, hpsi_a;
    copy_cols(psi, active_cols, psi_a);
    copy_cols(spsi_.data(), active_cols, spsi_a);
    copy_cols(hpsi_.data(), active_cols, hpsi_a);

    std::vector<T> s(nact * nact, T(0));
    gram(psi_a.data(), spsi_a.data(), nact, nact, s, nact);

    bool cholesky_ok = false;
    try
    {
        PpcgLapack<T>::potrf(nact, s.data());
        right_solve_upper(s, nact, psi_a);
        right_solve_upper(s, nact, spsi_a);
        right_solve_upper(s, nact, hpsi_a);
        cholesky_ok = is_s_orthonormal(psi_a.data(), spsi_a.data(), nact);
    }
    catch (const std::runtime_error&)
    {
        cholesky_ok = false;
    }

    if (!cholesky_ok)
        s_gram_schmidt(psi_a.data(), hpsi_a.data(), spsi_a.data(), nact);

    scatter_cols(psi, active_cols, psi_a);
    scatter_cols(spsi_.data(), active_cols, spsi_a);
    scatter_cols(hpsi_.data(), active_cols, hpsi_a);
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
    std::vector<T> hsub(n_band_ * n_band_, T(0));
    std::vector<T> ssub(n_band_ * n_band_, T(0));
    gram(psi, hpsi_.data(), n_band_, n_band_, hsub, n_band_);
    gram(psi, spsi_.data(), n_band_, n_band_, ssub, n_band_);

    std::vector<Real> eval(n_band_, static_cast<Real>(0));
    bool sygvd_ok = false;
    try
    {
        PpcgLapack<T>::hegvd(n_band_, hsub.data(), ssub.data(),
                                  eval.data());
        sygvd_ok = true;
    }
    catch (const std::runtime_error&)
    {
        // Fallback: diagonal Rayleigh quotients.
        // hsub and ssub may be corrupted by sygvd; re-form them.
        gram(psi, hpsi_.data(), n_band_, n_band_, hsub, n_band_);
        gram(psi, spsi_.data(), n_band_, n_band_, ssub, n_band_);
        for (int ii = 0; ii < n_band_; ++ii)
            eval[ii] = static_cast<Real>(std::real(hsub[ii + ii * n_band_]))
                     / std::max(static_cast<Real>(
                                    std::real(ssub[ii + ii * n_band_])),
                                static_cast<Real>(1e-30));
    }

    if (sygvd_ok)
    {
        std::vector<T> psi_old(psi, psi + ld_psi_ * n_band_);
        std::vector<T> spsi_old = spsi_;
        std::vector<T> hpsi_old = hpsi_;

        std::fill(psi, psi + ld_psi_ * n_band_, T(0));
        set_zero(spsi_);
        set_zero(hpsi_);

        for (int j = 0; j < n_band_; ++j)
        {
            for (int i = 0; i < n_band_; ++i)
            {
                const T c = hsub[i + j * n_band_];
                for (int ig = 0; ig < n_dim_; ++ig)
                {
                    psi[ idx(ig, j, ld_psi_)] += psi_old[ idx(ig, i, ld_psi_)] * c;
                    spsi_[idx(ig, j, ld_psi_)] += spsi_old[idx(ig, i, ld_psi_)] * c;
                    hpsi_[idx(ig, j, ld_psi_)] += hpsi_old[idx(ig, i, ld_psi_)] * c;
                }
            }
            eigenvalue[j] = eval[j];
        }
    }
    else
    {
        // No rotation: just update eigenvalues with Rayleigh quotients.
        for (int j = 0; j < n_band_; ++j)
            eigenvalue[j] = eval[j];
    }

    // Compute residual: w_i = H|psi_i> - eps_i * S|psi_i>
    set_zero(w_);
    for (int j = 0; j < n_band_; ++j)
        for (int ig = 0; ig < n_dim_; ++ig)
            w_[idx(ig, j, ld_psi_)] = hpsi_[idx(ig, j, ld_psi_)]
                                    - spsi_[idx(ig, j, ld_psi_)] * eigenvalue[j];

    lock_epairs(w_, ethr_band, active_cols);
}

// ---------------------------------------------------------------------------
// Trace of H|psi> within active columns
// ---------------------------------------------------------------------------
template <typename T, typename Device>
typename DiagoPPCG<T, Device>::Real
DiagoPPCG<T, Device>::trace_of_active_projected(
    const T* psi, const std::vector<int>& active_cols) const
{
    if (active_cols.empty())
        return static_cast<Real>(0);

    std::vector<T> psi_a, hpsi_a;
    copy_cols(psi, active_cols, psi_a);
    copy_cols(hpsi_.data(), active_cols, hpsi_a);

    const int nact = static_cast<int>(active_cols.size());
    std::vector<T> g(nact * nact, T(0));
    gram(psi_a.data(), hpsi_a.data(), nact, nact, g, nact);

    Real tr = 0;
    for (int i = 0; i < nact; ++i)
        tr += static_cast<Real>(std::real(g[i + i * nact]));
    return tr;
}

//==============================================================================
// CONJUGATE_GRADIENT STRATEGY
//==============================================================================

// ---------------------------------------------------------------------------
// Compute gradient: grad_i = H|psi_i> - eps_i * S|psi_i>
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_gradient(
    const Real* /*prec*/,
    const T* hpsi,
    const T* spsi,
    const T* /*psi*/,
    const Real* eigenvalue,
    std::vector<T>& grad) const
{
    grad.assign(ld_psi_ * n_band_, T(0));
    for (int j = 0; j < n_band_; ++j)
    {
        const Real ej = eigenvalue[j];
        for (int ig = 0; ig < n_dim_; ++ig)
            grad[idx(ig, j, ld_psi_)] = hpsi[idx(ig, j, ld_psi_)]
                                      - spsi[idx(ig, j, ld_psi_)] * ej;
    }
}

// ---------------------------------------------------------------------------
// Orthogonalize gradient: grad_j -= sum_i <psi_i|grad_j> * S|psi_i>
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::orth_gradient(
    const T* psi, const T* spsi,
    std::vector<T>& grad) const
{
    for (int j = 0; j < n_band_; ++j)
    {
        for (int i = 0; i < n_band_; ++i)
        {
            // Full complex inner product <psi_i | grad_j>
            T coeff = 0;
            const T* pi = psi + i * ld_psi_;
            const T* gj = grad.data() + j * ld_psi_;
            for (int ig = 0; ig < n_dim_; ++ig)
                coeff += std::conj(pi[ig]) * gj[ig];
            if (std::abs(coeff) <= std::numeric_limits<Real>::epsilon())
                continue;
            // grad_j -= S|psi_i> * coeff
            const T* si = spsi + i * ld_psi_;
            T* gj_out = grad.data() + j * ld_psi_;
            for (int ig = 0; ig < n_dim_; ++ig)
                gj_out[ig] -= si[ig] * coeff;
        }
    }
}

// ---------------------------------------------------------------------------
// Polak-Ribiere conjugate gradient update with preconditioning:
//   z_new = -P^{-1} * r_new
//   beta = max(0, <z_new, r_new - r_old> / <z_old, r_old>)
//   d_new = z_new + beta * d_old
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::update_polak_ribiere(
    const std::vector<T>& grad,
    std::vector<T>& p,
    std::vector<T>& grad_old,
    std::vector<T>& z_old,
    std::vector<Real>& beta_denom,
    const Real* prec) const
{
    const bool first_iter = p.empty();
    if (first_iter)
    {
        p.assign(ld_psi_ * n_band_, T(0));
        z_old.assign(ld_psi_ * n_band_, T(0));
        beta_denom.assign(n_band_, std::numeric_limits<Real>::infinity());
    }

    std::vector<T> z_new(ld_psi_ * n_band_, T(0));

    for (int j = 0; j < n_band_; ++j)
    {
        const T* g  = grad.data() + j * ld_psi_;
        T* pj  = p.data() + j * ld_psi_;
        T* zn  = z_new.data() + j * ld_psi_;
        T* zo  = z_old.data() + j * ld_psi_;

        Real beta_num_zr = 0;
        Real beta_num_zo = 0;

        for (int ig = 0; ig < n_dim_; ++ig)
        {
            // z_new = -P^{-1} * grad
            T z = -g[ig] / std::max(prec[ig], static_cast<Real>(1.0e-12));
            zn[ig] = z;

            // r_old = -P * z_old (recover old raw residual)
            T r_old = -prec[ig] * zo[ig];

            beta_num_zr += static_cast<Real>(std::real(z * std::conj(g[ig])));
            beta_num_zo += static_cast<Real>(std::real(z * std::conj(r_old)));
        }

        Real beta = 0;
        const Real denom = beta_denom[j];
        if (denom > static_cast<Real>(1.0e-30))
        {
            beta = (beta_num_zr - beta_num_zo) / denom;
            if (beta < 0)
                beta = 0;
        }

        // d_new = z_new + beta * d_old
        for (int ig = 0; ig < n_dim_; ++ig)
            pj[ig] = zn[ig] + beta * pj[ig];

        // Save <z_new, r_new> as denominator for next iteration.
        beta_denom[j] = beta_num_zr + static_cast<Real>(1.0e-30);
    }

    // Persist state for next iteration.
    z_old.swap(z_new);
    grad_old = grad;
}

// ---------------------------------------------------------------------------
// Line minimization along search direction:
//   For each band j: find optimal step alpha by minimizing the Rayleigh quotient
//   in the 2D subspace spanned by |psi_j> and |p_j>.
//
//   The Rayleigh quotient:
//     R(alpha) = (h_ii + 2 alpha h_ip + alpha^2 h_pp)
//              / (s_ii + 2 alpha s_ip + alpha^2 s_pp)
//
//   Setting dR/dalpha = 0 gives a QUADRATIC equation
//   A alpha^2 + B alpha + C = 0 with:
//     A = s_ip * h_pp - h_ip * s_pp
//     B = s_ii * h_pp - h_ii * s_pp
//     C = s_ii * h_ip - h_ii * s_ip
//
//   The linear approximation alpha = -C / B (dropping the alpha^2 term)
//   picks one of
//   the two stationary points more-or-less arbitrarily.  For bands far from
//   convergence this can select the MAXIMUM, driving psi toward high-energy
//   states.  We solve the full quadratic and explicitly pick the root with
//   the lower Rayleigh quotient.
//
//   Update: |psi>  += alpha |p>
//           H|psi> += alpha H|p>
//           S|psi> += alpha S|p>
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::line_minimize(
    T* psi, T* hpsi, T* spsi,
    const T* p, const T* hp, const T* sp,
    int ncol) const
{
    for (int j = 0; j < ncol; ++j)
    {
        const int off = j * ld_psi_;
        T* pj  = psi  + off;
        T* hj  = hpsi + off;
        T* sj  = spsi + off;
        const T* pp = p   + off;
        const T* hpp = hp + off;
        const T* spp = sp + off;

        Real h_ii = gamma_dot(pj, hj);
        Real s_ii = gamma_dot(pj, sj);
        const T h_ip_c = complex_dot(pj, hpp);
        const T s_ip_c = complex_dot(pj, spp);
        Real h_pp = gamma_dot(pp, hpp);
        Real s_pp = gamma_dot(pp, spp);

        // Rotate the search direction so the first-order Rayleigh quotient
        // derivative is real. The scalar alpha solve below stays unchanged for
        // real problems, while complex PW states can use a complex step.
        T phase = T(1);
        const Real lambda = h_ii / std::max(s_ii, static_cast<Real>(1e-30));
        const T q = h_ip_c - T(lambda) * s_ip_c;
        const Real q_abs = std::abs(q);
        if (q_abs > static_cast<Real>(1e-30))
            phase = std::conj(q) / q_abs;

        Real h_ip = static_cast<Real>(std::real(phase * h_ip_c));
        Real s_ip = static_cast<Real>(std::real(phase * s_ip_c));

        // Coefficients of A alpha^2 + B alpha + C = 0
        const Real A = s_ip * h_pp - h_ip * s_pp;
        const Real B = s_ii * h_pp - h_ii * s_pp;
        const Real C = s_ii * h_ip - h_ii * s_ip;

        auto ray_quot = [&](Real a) -> Real {
            return (h_ii + static_cast<Real>(2) * a * h_ip + a * a * h_pp)
                 / std::max(s_ii + static_cast<Real>(2) * a * s_ip + a * a * s_pp,
                            static_cast<Real>(1e-30));
        };

        Real alpha = 0;
        Real alpha_linear = (std::abs(B) > static_cast<Real>(1e-30))
                          ? -C / B : static_cast<Real>(0);

        const Real tol = std::numeric_limits<Real>::epsilon() * static_cast<Real>(100);
        if (std::abs(A) > tol * std::max(static_cast<Real>(1), std::abs(B)))
        {
            const Real disc = B * B - static_cast<Real>(4) * A * C;
            if (disc >= static_cast<Real>(0))
            {
                const Real sqrt_disc = std::sqrt(disc);
                const Real a1 = (-B + sqrt_disc) / (static_cast<Real>(2) * A);
                const Real a2 = (-B - sqrt_disc) / (static_cast<Real>(2) * A);

                const Real r1 = ray_quot(a1);
                const Real r2 = ray_quot(a2);
                const Real r_lin = ray_quot(alpha_linear);

                if (r1 < r2 && r1 < r_lin)
                    alpha = a1;
                else if (r2 < r1 && r2 < r_lin)
                    alpha = a2;
                else
                    alpha = alpha_linear;
            }
            else
            {
                alpha = alpha_linear;
            }
        }
        else
        {
            alpha = alpha_linear;
        }

        for (int ig = 0; ig < n_dim_; ++ig)
        {
            const T step = T(alpha) * phase;
            pj[ig] += step * pp[ig];
            hj[ig] += step * hpp[ig];
            sj[ig] += step * spp[ig];
        }
    }
}

// ---------------------------------------------------------------------------
// Cholesky orthonormalization (S-orthonormal):
//   1. Form S-gram matrix J = psi^H * S * psi
//   2. Cholesky: J = U^T * U  (upper)
//   3. Invert U: U^{-1}
//   4. psi *= U^{-1},  Hpsi *= U^{-1},  Spsi *= U^{-1}
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::orth_cholesky(
    T* psi, T* hpsi, T* spsi, int ncol) const
{
    // Save original vectors in case Cholesky fails numerically.
    std::vector<T> psi_orig(psi, psi + ld_psi_ * ncol);
    std::vector<T> hpsi_orig(hpsi, hpsi + ld_psi_ * ncol);
    std::vector<T> spsi_orig(spsi, spsi + ld_psi_ * ncol);

    // Gram matrix of S-orthonormality: J_{ij} = <psi_i | S | psi_j>
    std::vector<T> gram_s(ncol * ncol, T(0));
    for (int j = 0; j < ncol; ++j)
        for (int i = 0; i < ncol; ++i)
            gram_s[i + j * ncol] = complex_dot(psi + i * ld_psi_,
                                                spsi + j * ld_psi_);

    bool cholesky_ok = false;
    try
    {
        PpcgLapack<T>::potrf(ncol, gram_s.data());
        PpcgLapack<T>::trtri(ncol, gram_s.data());

        std::vector<T> tmp(ld_psi_ * ncol, T(0));
        for (int j = 0; j < ncol; ++j)
            for (int i = 0; i < ncol; ++i) {
                const T uinv = gram_s[i + j * ncol];
                for (int ig = 0; ig < n_dim_; ++ig)
                    tmp[idx(ig, j, ld_psi_)] += psi[idx(ig, i, ld_psi_)] * uinv;
            }
        std::copy(tmp.begin(), tmp.end(), psi);

        set_zero(tmp);
        for (int j = 0; j < ncol; ++j)
            for (int i = 0; i < ncol; ++i) {
                const T uinv = gram_s[i + j * ncol];
                for (int ig = 0; ig < n_dim_; ++ig)
                    tmp[idx(ig, j, ld_psi_)] += hpsi[idx(ig, i, ld_psi_)] * uinv;
            }
        std::copy(tmp.begin(), tmp.end(), hpsi);

        set_zero(tmp);
        for (int j = 0; j < ncol; ++j)
            for (int i = 0; i < ncol; ++i) {
                const T uinv = gram_s[i + j * ncol];
                for (int ig = 0; ig < n_dim_; ++ig)
                    tmp[idx(ig, j, ld_psi_)] += spsi[idx(ig, i, ld_psi_)] * uinv;
            }
        std::copy(tmp.begin(), tmp.end(), spsi);

        cholesky_ok = is_s_orthonormal(psi, spsi, ncol);
    }
    catch (const std::runtime_error&) { cholesky_ok = false; }

    if (!cholesky_ok)
    {
        std::copy(psi_orig.begin(), psi_orig.end(), psi);
        std::copy(hpsi_orig.begin(), hpsi_orig.end(), hpsi);
        std::copy(spsi_orig.begin(), spsi_orig.end(), spsi);
        s_gram_schmidt(psi, hpsi, spsi, ncol);
    }
}

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

    validate_input(psi_in, eigenvalue_in, ethr_band, prec);
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

    std::vector<int> all_cols(ncol);
    std::iota(all_cols.begin(), all_cols.end(), 0);

    force_g0_real(psi_in, ncol);
    apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
    apply_s_current(psi_in, spsi_.data(), ncol);

    double avg_iter = 1.0;
    int iter = 1;
    std::vector<int> active_cols;

    std::ofstream residual_trace;
    if (const char* path = std::getenv("ABACUS_PPCG_RESIDUAL_TRACE"))
    {
        // Optional debug trace for plotting PPCG convergence curves.
        residual_trace.open(path);
        if (residual_trace)
            residual_trace << "iteration,stage,max_residual\n";
    }
    auto record_residual = [&](int iteration, const char* stage) {
        if (!residual_trace)
            return;
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

    // ---------------------------------------------------------------------------
    // Strategy dispatch
    // ---------------------------------------------------------------------------
    if (strategy_ == PpcgStrategy::BLOCK_SUBSPACE)
    {
        // Initialize with Rayleigh-Ritz.
        rayleigh_ritz(psi_in, eigenvalue_in, active_cols, ethr_band);
        // Recompute to keep hpsi/spi consistent with rotated psi.
        apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
        apply_s_current(psi_in, spsi_.data(), ncol);
        record_residual(0, "initial_rr");

        Real trG = trace_of_active_projected(psi_in, active_cols);
        Real trdif = static_cast<Real>(-1);

        while (!active_cols.empty() && iter <= maxiter_)
        {
            const int nact = static_cast<int>(active_cols.size());
            const int nsb = std::max(1, (nact + sbsize_ - 1) / sbsize_);
            const Real trtol = diag_thr_ * std::sqrt(static_cast<Real>(nact));

            // Precondition the residual.
            divide_by_preconditioner(active_cols, prec, w_);
            apply_s_current(w_.data(), sw_.data(), ncol);
            project_against(psi_in, spsi_.data(), all_cols, w_, sw_, active_cols);

            // Apply H to the search direction.
            std::vector<T> w_active;
            copy_cols(w_.data(), active_cols, w_active);
            force_g0_real(w_active.data(), nact);
            std::vector<T> hw_active(ld_psi_ * nact, T(0));
            scatter_cols(w_.data(), active_cols, w_active);
            apply_h(hpsi_func, w_active.data(), hw_active.data(), nact);
            scatter_cols(hw_.data(), active_cols, hw_active);
            apply_s_current(w_.data(), sw_.data(), ncol);

            avg_iter += static_cast<double>(nact) / static_cast<double>(ncol);

            // Use the 3-block [psi, w, p] subspace.
            // w and p are normalized to unit S-norm before building the
            // Gram matrix (see build_small_subspace), which keeps M
            // well-conditioned even when residuals are small.
            // When p is zero (first iteration) or nearly collinear with w,
            // we fall back to the 2-block subspace for this iteration;
            // update_one_block will still produce a valid p for the next
            // iteration from the w contribution.
            const bool use_p = true;
            bool use_p_now = use_p;
            if (use_p)
            {
                apply_s_current(p_.data(), sp_.data(), ncol);
                project_against(psi_in, spsi_.data(), all_cols, p_, sp_, active_cols);

                // Detect when p makes the subspace nearly rank-deficient:
                // p near-zero (first iteration, not yet built) or p nearly
                // collinear with w.  Either way the [w,p] block of the
                // Gram matrix becomes nearly singular.  We do NOT replace p
                // with H*w because H*w is close to lambda*w when w is
                // approximately an eigenvector. It does not fix the
                // collinearity. Instead
                // we simply skip p for this iteration.
                for (const int c : active_cols)
                {
                    Real p_nrm2 = 0, w_nrm2 = 0, pw_re = 0;
                    for (int ig = 0; ig < n_dim_; ++ig)
                    {
                        p_nrm2 += static_cast<Real>(std::norm(p_[idx(ig, c, ld_psi_)]));
                        w_nrm2 += static_cast<Real>(std::norm(w_[idx(ig, c, ld_psi_)]));
                        pw_re  += static_cast<Real>(
                            std::real(std::conj(p_[idx(ig, c, ld_psi_)])
                                      * w_[idx(ig, c, ld_psi_)]));
                    }
                    const Real denom = p_nrm2 * w_nrm2;
                    Real cos2 = -1;
                    if (denom > Real(1e-60))
                        cos2 = (pw_re * pw_re) / denom;
                    if (p_nrm2 <= Real(1e-30) ||
                        (denom > Real(1e-60) && cos2 > Real(0.99)))
                    {
                        use_p_now = false;
                        break;
                    }
                }
            }

            // Block subspace solve.
            for (int isb = 0; isb < nsb; ++isb)
            {
                const int i0 = isb * sbsize_;
                const int l = std::min(sbsize_, nact - i0);
                std::vector<int> cols(active_cols.begin() + i0,
                                      active_cols.begin() + i0 + l);

                SmallSubspace subspace;
                build_small_subspace(psi_in, cols, use_p_now, subspace);
                solve_small_generalized((use_p_now ? 3 : 2) * l, subspace);
                update_one_block(psi_in, cols, l, use_p_now, subspace);
            }

            // Periodic Rayleigh-Ritz.
            if (iter % rr_step_ == 0)
            {
                rayleigh_ritz(psi_in, eigenvalue_in, active_cols, ethr_band);
                apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
                apply_s_current(psi_in, spsi_.data(), ncol);
                trdif = static_cast<Real>(-1);
                trG = 0;
                for (const int c : active_cols)
                    trG += eigenvalue_in[c];
                record_residual(iter, "rayleigh_ritz");
            }
            else
            {
                chol_qr_active(psi_in, active_cols);

                // Compute updated eigenvalues and residuals.
                std::vector<T> psi_a, hpsi_a;
                copy_cols(psi_in, active_cols, psi_a);
                copy_cols(hpsi_.data(), active_cols, hpsi_a);

                const int na = static_cast<int>(active_cols.size());
                std::vector<T> ga(ncol * na, T(0));
                gram(psi_in, hpsi_a.data(), ncol, na, ga, ncol);

                set_zero(w_);
                for (int ja = 0; ja < na; ++ja)
                {
                    for (int ig = 0; ig < n_dim_; ++ig)
                    {
                        T sum = T(0);
                        for (int ia = 0; ia < ncol; ++ia)
                            sum += spsi_[idx(ig, ia, ld_psi_)] * ga[ia + ja * ncol];
                        w_[idx(ig, active_cols[ja], ld_psi_)] =
                            hpsi_a[idx(ig, ja, ld_psi_)] - sum;
                    }
                    eigenvalue_in[active_cols[ja]] =
                        static_cast<Real>(std::real(
                            ga[active_cols[ja] + ja * ncol]));
                }

                Real trG1 = 0;
                for (int ja = 0; ja < na; ++ja)
                    trG1 += static_cast<Real>(std::real(
                        ga[active_cols[ja] + ja * ncol]));

                trdif = std::abs(trG1 - trG);
                trG = trG1;

                lock_epairs(w_, ethr_band, active_cols);
                if (trdif >= 0 && trdif <= trtol)
                {
                    rayleigh_ritz(psi_in, eigenvalue_in, active_cols, ethr_band);
                    apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
                    apply_s_current(psi_in, spsi_.data(), ncol);
                    trdif = static_cast<Real>(-1);
                    record_residual(iter, "trace_rr");
                }
                else
                {
                    record_residual(iter, "block_update");
                }
            }

            ++iter;
        }

        if ((iter - 1) % rr_step_ != 0)
            rayleigh_ritz(psi_in, eigenvalue_in, active_cols, ethr_band);
        // Final consistency: ensure hpsi/spi match the converged psi.
        apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
        apply_s_current(psi_in, spsi_.data(), ncol);
        record_residual(iter - 1, "final");
    }
    else // CONJUGATE_GRADIENT
    {
        // Initialize with Rayleigh-Ritz, same as BLOCK_SUBSPACE.
        // Diagonal Rayleigh quotients are poor approximations for random
        // initial guesses; starting the CG loop with them produces wrong
        // gradients that drive the search toward high-energy bands.
        rayleigh_ritz(psi_in, eigenvalue_in, active_cols, ethr_band);
        apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
        apply_s_current(psi_in, spsi_.data(), ncol);
        record_residual(0, "initial_rr");

        std::vector<T> grad;
        calc_gradient(prec, hpsi_.data(), spsi_.data(), psi_in,
                      eigenvalue_in, grad);
        orth_gradient(psi_in, spsi_.data(), grad);

        std::vector<T> p;
        grad_old_.clear();
        z_old_.clear();
        beta_denom_.clear();
        update_polak_ribiere(grad, p, grad_old_, z_old_, beta_denom_, prec);

        // CG iteration loop.
        while (iter <= maxiter_)
        {
            // Apply H and S to search direction.
            std::vector<T> hp(ld_psi_ * ncol, T(0));
            std::vector<T> sp(ld_psi_ * ncol, T(0));
            apply_h(hpsi_func, p.data(), hp.data(), ncol);
            apply_s_current(p.data(), sp.data(), ncol);

            // Line minimization.
            line_minimize(psi_in, hpsi_.data(), spsi_.data(),
                          p.data(), hp.data(), sp.data(), ncol);

            const bool do_rr = (iter % rr_step_ == 0);
            if (do_rr)
            {
                // Rayleigh-Ritz: full subspace diagonalization.
                // We recompute H|psi> and S|psi> first because line_minimize
                // modified psi. We do NOT call orth_cholesky here; Cholesky
                // mixes bands through the upper-triangular U^{-1} factor,
                // contaminating low-energy bands with high-energy components
                // and driving the eigenvalues upward.
                apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
                apply_s_current(psi_in, spsi_.data(), ncol);

                std::vector<int> dummy_active;
                rayleigh_ritz(psi_in, eigenvalue_in, dummy_active, ethr_band);

                // Sync hpsi/spi to the rotated wavefunctions.
                apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
                apply_s_current(psi_in, spsi_.data(), ncol);

                // Reset PR state: the rotation changes the basis,
                // so old gradients / search directions are invalid.
                p.clear();
                grad_old_.clear();
                z_old_.clear();
                beta_denom_.clear();
                record_residual(iter, "rayleigh_ritz");
            }
            else
            {
                // Cholesky orthonormalization.
                orth_cholesky(psi_in, hpsi_.data(), spsi_.data(), ncol);

                // After Cholesky the bands are S-orthonormal, but the
                // upper-triangular U^{-1} transformation mixes high-energy
                // components into the low-energy bands.  Diagonal Rayleigh
                // quotients then overestimate the low eigenvalues and
                // produce wrong gradients that drive the CG search toward
                // high-energy states.
                //
                // Solve the subspace generalized eigenvalue problem to get
                // correct Ritz values. We do NOT rotate the states; that
                // would invalidate the Polak-Ribiere conjugate-direction
                // accumulators.  The Cholesky basis spans the same subspace,
                // so the Ritz values are exact for this subspace.
                std::vector<T> h_sub(ncol * ncol, T(0));
                std::vector<T> s_sub(ncol * ncol, T(0));
                for (int jj = 0; jj < ncol; ++jj)
                {
                    for (int ii = 0; ii < ncol; ++ii)
                    {
                        h_sub[ii + jj * ncol]
                            = complex_dot(psi_in + ii * ld_psi_,
                                          hpsi_.data() + jj * ld_psi_);
                        s_sub[ii + jj * ncol]
                            = complex_dot(psi_in + ii * ld_psi_,
                                          spsi_.data() + jj * ld_psi_);
                    }
                }

                std::vector<Real> eval_cg(ncol, static_cast<Real>(0));
                try
                {
                    PpcgLapack<T>::hegvd(ncol, h_sub.data(),
                                              s_sub.data(),
                                              eval_cg.data());
                }
                catch (const std::runtime_error&)
                {
                    // Fallback: diagonal Rayleigh quotients.
                    // h_sub and s_sub may be corrupted by sygvd; re-form them.
                    for (int jj = 0; jj < ncol; ++jj)
                    {
                        for (int ii = 0; ii < ncol; ++ii)
                        {
                            h_sub[ii + jj * ncol]
                                = complex_dot(psi_in + ii * ld_psi_,
                                              hpsi_.data() + jj * ld_psi_);
                            s_sub[ii + jj * ncol]
                                = complex_dot(psi_in + ii * ld_psi_,
                                              spsi_.data() + jj * ld_psi_);
                        }
                    }
                    for (int ii = 0; ii < ncol; ++ii)
                        eval_cg[ii] =
                            static_cast<Real>(std::real(h_sub[ii + ii * ncol]))
                            / std::max(static_cast<Real>(
                                           std::real(s_sub[ii + ii * ncol])),
                                       static_cast<Real>(1e-30));
                }
                for (int ii = 0; ii < ncol; ++ii)
                    eigenvalue_in[ii] = eval_cg[ii];
                record_residual(iter, "cg_step");
            }

            // Compute new gradient.
            calc_gradient(prec, hpsi_.data(), spsi_.data(), psi_in,
                          eigenvalue_in, grad);
            orth_gradient(psi_in, spsi_.data(), grad);

            // Polak-Ribiere update.
            update_polak_ribiere(grad, p, grad_old_, z_old_, beta_denom_, prec);

            // Convergence check.
            bool all_converged = true;
            for (int i = 0; i < ncol; ++i)
            {
                Real nrm2 = 0;
                for (int ig = 0; ig < n_dim_; ++ig)
                    nrm2 += static_cast<Real>(
                        std::norm(grad[idx(ig, i, ld_psi_)]));
                if (std::sqrt(nrm2) > std::max(static_cast<Real>(ethr_band[i]),
                                               diag_thr_))
                {
                    all_converged = false;
                    break;
                }
            }
            if (all_converged)
                break;

            ++iter;
        }

        avg_iter = static_cast<double>(iter);
    }

    return avg_iter;
}

// =============================================================================
// Explicit template instantiation (CPU only; extend for GPU as needed)
// =============================================================================
template class DiagoPPCG<std::complex<float>,  base_device::DEVICE_CPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_CPU>;

} // namespace hsolver
