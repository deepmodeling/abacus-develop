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
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ * n_band_ > 4096)
#endif
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
    std::vector<T> coeff(n_band_ * n_band_, T(0));
    gram(psi, grad.data(), n_band_, n_band_, coeff, n_band_);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ * n_band_ > 4096)
#endif
    for (int j = 0; j < n_band_; ++j)
    {
        for (int i = 0; i < n_band_; ++i)
        {
            const T cproj = coeff[i + j * n_band_];
            if (std::abs(cproj) <= std::numeric_limits<Real>::epsilon())
                continue;
            // grad_j -= S|psi_i> * coeff
            const T* si = spsi + i * ld_psi_;
            T* gj_out = grad.data() + j * ld_psi_;
            for (int ig = 0; ig < n_dim_; ++ig)
                gj_out[ig] -= si[ig] * cproj;
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
    std::vector<double> beta_nums(2 * n_band_, 0.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ * n_band_ > 4096)
#endif
    for (int j = 0; j < n_band_; ++j)
    {
        const T* g  = grad.data() + j * ld_psi_;
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
        beta_nums[j] = static_cast<double>(beta_num_zr);
        beta_nums[n_band_ + j] = static_cast<double>(beta_num_zo);
    }
    reduce_pool_if_mpi_ready(beta_nums.data(), static_cast<int>(beta_nums.size()));

    for (int j = 0; j < n_band_; ++j)
    {
        T* pj  = p.data() + j * ld_psi_;
        T* zn  = z_new.data() + j * ld_psi_;
        const Real beta_num_zr = static_cast<Real>(beta_nums[j]);
        const Real beta_num_zo = static_cast<Real>(beta_nums[n_band_ + j]);
        Real beta = 0;
        const Real denom = beta_denom[j];
        if (denom > static_cast<Real>(1.0e-30))
        {
            beta = (beta_num_zr - beta_num_zo) / denom;
            if (beta < 0)
                beta = 0;
        }

        // d_new = z_new + beta * d_old
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ > 4096)
#endif
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
    std::vector<double> real_coeffs(4 * ncol, 0.0);
    std::vector<T> mixed_coeffs(2 * ncol, T(0));

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ * ncol > 4096)
#endif
    for (int j = 0; j < ncol; ++j)
    {
        const int off = j * ld_psi_;
        const T* pj = psi + off;
        const T* hj = hpsi + off;
        const T* sj = spsi + off;
        const T* pp = p + off;
        const T* hpp = hp + off;
        const T* spp = sp + off;

        Real h_ii = 0;
        Real s_ii = 0;
        Real h_pp = 0;
        Real s_pp = 0;
        T h_ip = T(0);
        T s_ip = T(0);

        for (int ig = 0; ig < n_dim_; ++ig)
        {
            h_ii += static_cast<Real>(std::real(std::conj(pj[ig]) * hj[ig]));
            s_ii += static_cast<Real>(std::real(std::conj(pj[ig]) * sj[ig]));
            h_ip += std::conj(pj[ig]) * hpp[ig];
            s_ip += std::conj(pj[ig]) * spp[ig];
            h_pp += static_cast<Real>(std::real(std::conj(pp[ig]) * hpp[ig]));
            s_pp += static_cast<Real>(std::real(std::conj(pp[ig]) * spp[ig]));
        }

        real_coeffs[j] = static_cast<double>(h_ii);
        real_coeffs[ncol + j] = static_cast<double>(s_ii);
        real_coeffs[2 * ncol + j] = static_cast<double>(h_pp);
        real_coeffs[3 * ncol + j] = static_cast<double>(s_pp);
        mixed_coeffs[j] = h_ip;
        mixed_coeffs[ncol + j] = s_ip;
    }

    reduce_pool_if_mpi_ready(real_coeffs.data(), static_cast<int>(real_coeffs.size()));
    reduce_pool_if_mpi_ready(mixed_coeffs.data(), static_cast<int>(mixed_coeffs.size()));

    std::vector<T> steps(ncol, T(0));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ncol > 16)
#endif
    for (int j = 0; j < ncol; ++j)
    {
        Real h_ii = static_cast<Real>(real_coeffs[j]);
        Real s_ii = static_cast<Real>(real_coeffs[ncol + j]);
        const T h_ip_c = mixed_coeffs[j];
        const T s_ip_c = mixed_coeffs[ncol + j];
        Real h_pp = static_cast<Real>(real_coeffs[2 * ncol + j]);
        Real s_pp = static_cast<Real>(real_coeffs[3 * ncol + j]);

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

        steps[j] = T(alpha) * phase;
    }

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) if (n_dim_ * ncol > 4096)
#endif
    for (int j = 0; j < ncol; ++j)
    {
        for (int ig = 0; ig < n_dim_; ++ig)
        {
            const int off = idx(ig, j, ld_psi_);
            psi[off] += steps[j] * p[off];
            hpsi[off] += steps[j] * hp[off];
            spsi[off] += steps[j] * sp[off];
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
    std::vector<T> gram_s;
    gram(psi, spsi, ncol, ncol, gram_s, ncol);

    bool cholesky_ok = false;
    try
    {
        HermitianLapack<T>::potrf(ncol, gram_s.data());
        HermitianLapack<T>::trtri(ncol, gram_s.data());

        const T one = T(1);
        const T zero = T(0);
        std::vector<T> tmp(ld_psi_ * ncol, T(0));
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         n_dim_,
                                         ncol,
                                         ncol,
                                         &one,
                                         psi,
                                         ld_psi_,
                                         gram_s.data(),
                                         ncol,
                                         &zero,
                                         tmp.data(),
                                         ld_psi_);
        std::copy(tmp.begin(), tmp.end(), psi);

        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         n_dim_,
                                         ncol,
                                         ncol,
                                         &one,
                                         hpsi,
                                         ld_psi_,
                                         gram_s.data(),
                                         ncol,
                                         &zero,
                                         tmp.data(),
                                         ld_psi_);
        std::copy(tmp.begin(), tmp.end(), hpsi);

        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         n_dim_,
                                         ncol,
                                         ncol,
                                         &one,
                                         spsi,
                                         ld_psi_,
                                         gram_s.data(),
                                         ncol,
                                         &zero,
                                         tmp.data(),
                                         ld_psi_);
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
