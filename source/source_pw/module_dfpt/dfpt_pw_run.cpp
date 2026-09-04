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
    for (int iat = 0; iat < nat; ++iat)
    {
        for (int idir = 0; idir < 3; ++idir)
        {
            // 2n+1 accumulation of this converged displacement
            phon_.accumulate_electron(q_idx, iat, idir, gs_psi_, wg_, data_);
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
    phon_.diagonalize(q_idx, data_);
    if (q_idx != 0 || !data_.get_loto())
    {
        return;
    }
    // non-analytic LO-TO correction along the data-layer direction
    // (default isotropic (1,1,1)/sqrt(3) for cubic crystals;
    // set_loto_dir overrides, e.g. per irrep direction in stage A)
    phon_.add_loto(data_.get_loto_dir(), data_);
    phon_.diagonalize_loto(data_);
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
