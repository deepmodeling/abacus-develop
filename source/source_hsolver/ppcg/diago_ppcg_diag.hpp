#include <cstdlib>
#include <fstream>

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
            const int nact = static_cast<int>(active_cols.size());
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

            avg_iter += static_cast<double>(nact) / static_cast<double>(ncol);

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
            apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
            apply_s_current(psi_in, spsi_.data(), ncol);
            record_residual(iter, "rayleigh_ritz");

            ++iter;
        }

        // Final consistency: ensure hpsi/spi match the converged psi.
        apply_h(hpsi_func, psi_in, hpsi_.data(), ncol);
        apply_s_current(psi_in, spsi_.data(), ncol);
        record_residual(iter - 1, "final");
    }
    else // CONJUGATE_GRADIENT
    {
        // Initialize with Rayleigh-Ritz — same as BLOCK_SUBSPACE.
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
        z_old_.clear();
        beta_denom_.clear();
        update_polak_ribiere(grad, p, z_old_, beta_denom_, prec);

        // CG iteration loop.
        std::vector<T> hp;
        std::vector<T> sp;
        hp.reserve(sz);
        sp.reserve(sz);
        while (iter <= maxiter_)
        {
            // Apply H and S to search direction.
            hp.assign(ld_psi_ * ncol, T(0));
            sp.assign(ld_psi_ * ncol, T(0));
            apply_h(hpsi_func, p.data(), hp.data(), ncol);
            apply_s_current(p.data(), sp.data(), ncol);

            // Line minimization.
            line_minimize(psi_in, hpsi_.data(), spsi_.data(),
                          p.data(), hp.data(), sp.data(), ncol);

            const bool do_rr = (iter % rr_step_) == 0;
            if (do_rr)
            {
                // Rayleigh-Ritz: full subspace diagonalization.
                // We recompute H|psi> and S|psi> first because line_minimize
                // modified psi.  We do NOT call orth_cholesky here — Cholesky
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
                // correct Ritz values.  We do NOT rotate the states — that
                // would invalidate the Polak-Ribiere conjugate-direction
                // accumulators.  The Cholesky basis spans the same subspace,
                // so the Ritz values are exact for this subspace.
                std::vector<T> h_sub(ncol * ncol, T(0));
                std::vector<T> s_sub(ncol * ncol, T(0));
                gram(psi_in, hpsi_.data(), ncol, ncol, h_sub, ncol);
                gram(psi_in, spsi_.data(), ncol, ncol, s_sub, ncol);

                std::vector<Real> eval_cg(ncol, static_cast<Real>(0));
                try
                {
                    HermitianLapack<T>::sygvd(ncol, h_sub.data(),
                                              s_sub.data(),
                                              eval_cg.data());
                }
                catch (const std::runtime_error&)
                {
                    // Fallback: diagonal Rayleigh quotients.
                    // h_sub and s_sub may be corrupted by sygvd; re-form them.
                    gram(psi_in, hpsi_.data(), ncol, ncol, h_sub, ncol);
                    gram(psi_in, spsi_.data(), ncol, ncol, s_sub, ncol);
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
            update_polak_ribiere(grad, p, z_old_, beta_denom_, prec);

            // Convergence check.
            bool all_converged = true;
            std::vector<double> grad_nrm2(ncol, 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ * ncol > 4096)
#endif
            for (int i = 0; i < ncol; ++i)
            {
                double nrm2 = 0.0;
                for (int ig = 0; ig < n_dim_; ++ig)
                    nrm2 += static_cast<double>(
                        std::norm(grad[idx(ig, i, ld_psi_)]));
                grad_nrm2[i] = nrm2;
            }
            reduce_pool_if_mpi_ready(grad_nrm2.data(), ncol);
            for (int i = 0; i < ncol; ++i)
            {
                if (std::sqrt(static_cast<Real>(grad_nrm2[i]))
                    > std::max(static_cast<Real>(ethr_band[i]), diag_thr_))
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

} // namespace hsolver
