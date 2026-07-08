#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace hsolver {

// ---------------------------------------------------------------------------
// Check S-orthonormality of a column block.
// ---------------------------------------------------------------------------
template <typename T, typename Device>
bool DiagoPPCG<T, Device>::is_s_orthonormal(
    const T* psi, const T* spsi, int ncol) const
{
    const Real orth_tol = static_cast<Real>(10)
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ > 4096)
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
            static_cast<Real>(1e-30)));
        Real inv_nrm = static_cast<Real>(1) / nrm;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ > 4096)
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
            rr_eval_[ii] = static_cast<Real>(std::real(rr_hsub_[ii + ii * n_band_]))
                     / std::max(static_cast<Real>(
                                    std::real(rr_ssub_[ii + ii * n_band_])),
                                static_cast<Real>(1e-30));
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
            eigenvalue[j] = rr_eval_[j];
    }

    // Compute residual: w_i = H|psi_i> - eps_i * S|psi_i>
    set_zero(w_);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) if (n_dim_ * n_band_ > 4096)
#endif
    for (int j = 0; j < n_band_; ++j)
        for (int ig = 0; ig < n_dim_; ++ig)
            w_[idx(ig, j, ld_psi_)] = hpsi_[idx(ig, j, ld_psi_)]
                                    - spsi_[idx(ig, j, ld_psi_)] * eigenvalue[j];

    lock_epairs(w_, ethr_band, active_cols);
}

} // namespace hsolver
