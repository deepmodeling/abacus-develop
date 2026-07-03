namespace hsolver {

// ---------------------------------------------------------------------------
// Back-substitute with upper triangular Cholesky factor: X *= R^{-1}
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::right_solve_upper(
    const std::vector<T>& r, int n, std::vector<T>& x) const
{
    std::vector<T> b = x;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ * n > 4096)
#endif
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
        HermitianLapack<T>::potrf(nact, s.data());
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
        HermitianLapack<T>::sygvd(n_band_, hsub.data(), ssub.data(),
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

        const T one = T(1);
        const T zero = T(0);
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         n_dim_,
                                         n_band_,
                                         n_band_,
                                         &one,
                                         psi_old.data(),
                                         ld_psi_,
                                         hsub.data(),
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
                                         spsi_old.data(),
                                         ld_psi_,
                                         hsub.data(),
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
                                         hpsi_old.data(),
                                         ld_psi_,
                                         hsub.data(),
                                         n_band_,
                                         &zero,
                                         hpsi_.data(),
                                         ld_psi_);

        for (int j = 0; j < n_band_; ++j)
        {
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
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) if (n_dim_ * n_band_ > 4096)
#endif
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

} // namespace hsolver
