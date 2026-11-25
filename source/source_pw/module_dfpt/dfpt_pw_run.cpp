// ============================================================
// DFPT_PW::run driver dispatcher implementation with helper
// extraction, moved from dfpt_pw.cpp so the driver TU stays below
// the coding-rule 500-line budget. The helpers mirror Gonze-Lee
// (1992, 1997) order: q0 response -> occ projector -> position legs
// -> e-field SCF -> displacement solves -> 2n+1 assembly -> Born
// charges -> phonon matrix diag + LO-TO.
// ============================================================

#include "dfpt_pw_impl.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"

#include <algorithm>

namespace ModuleDFPT
{

void DFPT_PW::Impl::run_q0_pre(int q_idx)
{
    // C7 note: the uniform E-field / position-operator responses of the
    // periodic crystal are ill-defined as bare matrix elements, so they
    // are obtained instead as the well-defined periodic commutator
    // [H_SCf, r] (QE dfpt_kernel / dfpt_tetra / dvpsi_e layout). The
    // compute_q0_response kernel stashes the irrep info consumed by the
    // per-direction solves below.
    if (data_.get_compute_q0())
    {
        q0_.compute_q0_response(data_);
    }
    if (wired())
    {
        build_occ_kq(q_idx);
    }
    if (q_idx != 0 || !data_.get_compute_q0() || !wired())
    {
        return;
    }
    // position legs of the screened Born charges: the q = 0 Y solves
    // need the projector just built and must land before the two-pass
    // displacement solves below reuse the shifted-operator context
    solve_pos_resp(q_idx);
    // SCF E-field responses of the dielectric tensor: after the
    // bare Y legs they consume, before the displacement solves
    // reuse the slots; the epsilon contraction runs straight after
    // (QE solve_e -> dielec.f90 order)
    solve_efield_resp(q_idx);
    q0_.compute_eps(wg_, data_);
}

double DFPT_PW::Impl::run_displacement_irrep_pass(int q_idx, int irrep)
{
    // Per-irrep self-consistent outer pass. Ledger semantics (B4): one
    // outer pass solves every displacement to its own convergence
    // (solve_displacement restarts each from a zero input density), and
    // the pass residual is the worst final displacement residual. An
    // unconverged pass therefore re-runs the full solve, bounded by
    // max_iter_ outer passes.
    if (!wired())
    {
        // design-phase skeleton: no bases wired, converge at once
        data_.add_residual(q_idx, irrep, 0.0);
        data_.set_converged(q_idx, irrep, true);
        return 0.0;
    }
    const int nat = ucell_->nat;
    // two passes over the 3N displacement basis: first solve every
    // displacement to convergence (the 2n+1 accumulation of
    // displacement b needs the converged dpsi AND screened potential
    // of every column displacement a), then run the 2n+1 accumulation
    // for each
    double worst = 0.0;
    for (int iat = 0; iat < nat; ++iat)
    {
        for (int idir = 0; idir < 3; ++idir)
        {
            const double residual = solve_displacement(q_idx, iat, idir);
            worst = std::max(worst, residual);
        }
    }
    // nlcc: the 2n+1 kernel term needs the first-order XC potential of
    // each displacement's TOTAL response (valence + core). Precompute the
    // nlcc flag once per pass; the per-displacement core response is
    // recomputed here because the working drhoc slot holds the LAST
    // solved displacement, and the valence response comes from the
    // per-displacement stash.
    bool nlcc_any = false;
    for (int it = 0; it < ucell_->ntype; ++it)
    {
        if (ucell_->atoms[it].ncpp.nlcc)
        {
            nlcc_any = true;
            break;
        }
    }
    const int npw = pw_rho_->npw;
    const int nrxx = pw_rho_->nrxx;
    const ModuleBase::Vector3<double> q_cart = data_.get_qvec(q_idx) * ucell_->G;
    // per-q table of f_xc.(drho_val + drhoc) for every displacement slot:
    // the total-response first-order XC potential that the Hermitized
    // kernel term contracts against the core responses (QE dynmat_nlcc
    // semantics; built once, before the accumulation pass)
    std::vector<std::vector<std::complex<double>>> dvxc_core;
    if (nlcc_any && xc_ != nullptr)
    {
        // fill EVERY displacement slot: the kernel term X(i,j) contracts
        // the row's dvaux with the COLUMN's core response, so rows of
        // non-nlcc atoms need their dvaux too (f_xc applied to their
        // valence response alone -- they have no core piece). Leaving
        // their slots empty silently dropped the entire X(N,B) cross
        // block (the dominant mixed piece for a B-N cell).
        dvxc_core.resize(3 * nat);
        std::vector<std::complex<double>> drhoc_g;
        for (int iat = 0; iat < nat; ++iat)
        {
            const bool iat_nlcc = ucell_->atoms[ucell_->iat2it[iat]].ncpp.nlcc;
            for (int idir = 0; idir < 3; ++idir)
            {
                std::vector<std::complex<double>> a_r(nrxx, std::complex<double>(0.0, 0.0));
                if (iat_nlcc)
                {
                    rho_.compute_drhoc_g(*ucell_, iat, idir, q_cart, drhoc_g);
                    if (static_cast<int>(drhoc_g.size()) == npw)
                    {
                        pw_rho_->recip2real(drhoc_g.data(), a_r.data());
                    }
                }
                // add the stashed converged valence response of this
                // displacement (the working slots hold the last one)
                const std::vector<std::complex<double>> drho_a_g = data_.get_drho_disp(iat, idir);
                if (static_cast<int>(drho_a_g.size()) == npw)
                {
                    std::vector<std::complex<double>> dval_r(nrxx);
                    pw_rho_->recip2real(drho_a_g.data(), dval_r.data());
                    for (int ir = 0; ir < nrxx; ++ir)
                    {
                        a_r[ir] += dval_r[ir];
                    }
                }
                xc_->set_context(q_cart);
                xc_->apply(a_r, dvxc_core[3 * iat + idir]);
            }
        }
        phon_.set_dvxc_core_table(&dvxc_core);
    }
    for (int iat = 0; iat < nat; ++iat)
    {
        for (int idir = 0; idir < 3; ++idir)
        {
            std::vector<std::complex<double>> dvxc_r;
            if (nlcc_any && xc_ != nullptr)
            {
                const int slot = 3 * iat + idir;
                if (slot < static_cast<int>(dvxc_core.size())
                    && static_cast<int>(dvxc_core[slot].size()) == nrxx)
                {
                    dvxc_r = dvxc_core[slot];
                }
            }
            // 2n+1 accumulation of this converged displacement
            phon_.accumulate_electron(q_idx, iat, idir, gs_psi_, wg_, data_, dvxc_r);
        }
    }
    data_.add_residual(q_idx, irrep, worst);
    data_.set_converged(q_idx, irrep, worst < data_.get_conv_thr());
    return worst;
}

void DFPT_PW::Impl::run_q0_post(int q_idx)
{
    // Screened Born charges: the Gonze-Lee 2n+1 form consumes the
    // converged (screened) dpsi of every q = 0 displacement stashed by
    // solve_displacement, so it must run after the two-pass solves
    // above and before the LO-TO term below consumes it.
    if (q_idx != 0 || !data_.get_compute_q0() || !wired())
    {
        return;
    }
    q0_.compute_born(gs_psi_, wg_, eig_, data_);
}

void DFPT_PW::Impl::run_assemble(int q_idx)
{
    phon_.assemble(q_idx, data_);
    // translational sum-rule audit of the raw Gamma matrix (acoustic-branch
    // diagnostics; the loto-corrected matrix is audited separately below).
    // The raw audit does not consume frequencies, so the ordinary
    // eigensolve runs once after any applicable ASR correction
    if (asr_type_ != "no")
    {
        phon_.audit_asr("raw", q_idx, data_);
        // enforce the sum rule on the raw matrix (with the Z* sum rule),
        // diagonalize once and audit; the LO-TO term below then consumes
        // the corrected Born charges (QE dynmat set_asr order)
        phon_.apply_asr(asr_type_, q_idx, data_);
        phon_.diagonalize(q_idx, data_);
        phon_.audit_asr("asr", q_idx, data_);
    }
    else
    {
        phon_.diagonalize(q_idx, data_);
        phon_.audit_asr("raw", q_idx, data_);
    }
    if (q_idx != 0 || !data_.get_loto())
    {
        return;
    }
    // non-analytic LO-TO correction along the data-layer direction
    // (default isotropic (1,1,1)/sqrt(3) for cubic crystals;
    // set_loto_dir overrides, e.g. per irrep direction in stage A)
    phon_.add_loto(data_.get_loto_dir(), data_);
    phon_.diagonalize_loto(data_);
    phon_.audit_asr("loto", q_idx, data_);
}

void DFPT_PW::run()
{
    ModuleBase::TITLE("DFPT_PW", "run");
    ModuleBase::timer::start("DFPT_PW", "run");
    const int nq = pimpl_->qlist_.get_nq();
    for (int q_idx = 0; q_idx < nq; ++q_idx)
    {
        pimpl_->run_q0_pre(q_idx);

        const int nirr = pimpl_->data_.get_nirr(q_idx);
        for (int irrep = 0; irrep < nirr; ++irrep)
        {
            pimpl_->data_.set_converged(q_idx, irrep, false);
            pimpl_->data_.set_current_iter(q_idx, irrep, 0);
            while (!pimpl_->data_.get_converged(q_idx, irrep)
                   && pimpl_->data_.get_current_iter(q_idx, irrep) < pimpl_->max_iter_)
            {
                pimpl_->run_displacement_irrep_pass(q_idx, irrep);
                pimpl_->data_.set_current_iter(q_idx, irrep, pimpl_->data_.get_current_iter(q_idx, irrep) + 1);
            }
        }

        pimpl_->run_q0_post(q_idx);
        pimpl_->run_assemble(q_idx);
    }
    ModuleBase::timer::end("DFPT_PW", "run");
}

} // namespace ModuleDFPT
