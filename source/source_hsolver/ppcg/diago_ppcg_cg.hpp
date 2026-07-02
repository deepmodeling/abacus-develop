namespace hsolver {

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
            const T* pi = psi + i * ld_psi_;
            const T* gj = grad.data() + j * ld_psi_;
            const T coeff = complex_dot(pi, gj);
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
        reduce_pool_if_mpi_ready(beta_num_zr);
        reduce_pool_if_mpi_ready(beta_num_zo);

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
//   For each band j: find optimal step α by minimizing the Rayleigh quotient
//   in the 2D subspace spanned by |psi_j> and |p_j>.
//
//   The Rayleigh quotient:
//     R(α) = (h_ii + 2α h_ip + α² h_pp) / (s_ii + 2α s_ip + α² s_pp)
//
//   Setting dR/dα = 0 gives a QUADRATIC equation A α² + B α + C = 0 with:
//     A = s_ip * h_pp - h_ip * s_pp
//     B = s_ii * h_pp - h_ii * s_pp
//     C = s_ii * h_ip - h_ii * s_ip
//
//   The linear approximation α = -C / B (dropping the α² term) picks one of
//   the two stationary points more-or-less arbitrarily.  For bands far from
//   convergence this can select the MAXIMUM, driving ψ toward high-energy
//   states.  We solve the full quadratic and explicitly pick the root with
//   the lower Rayleigh quotient.
//
//   Update: |psi>  += α |p>
//           H|psi> += α H|p>
//           S|psi> += α S|p>
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
        HermitianLapack<T>::potrf(ncol, gram_s.data());
        HermitianLapack<T>::trtri(ncol, gram_s.data());

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

} // namespace hsolver
