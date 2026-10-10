#include "source_esolver/esolver_lr_lcao_tddft.h"
#include "source_lcao/module_lr/zeq_solver.h"
#include "source_lcao/module_lr/cal_edm.h"
#include "source_lcao/module_lr/lr_force.h"
#include "source_lcao/module_lr/gradient_inputs.h"
#include "source_lcao/module_lr/gradient_output.h"
#include "source_lcao/module_lr/lr_amp.h"
#include "source_lcao/module_lr/grad_degen.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include <algorithm>
#include <complex>
#include <numeric>
#include "source_estate/module_dm/dm_from_psi.h"
#include "source_io/module_output/output_log.h"

using namespace LR;


template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::setup_relax_target_()
{
    this->excited_relax_ = (this->inp_->calculation == "relax");
    if (!this->excited_relax_) { return; }

    // The Z-vector equation has no complex solver (see zeq_solver.hpp), so an
    // excited-state gradient only exists at gamma. Fail here rather than after the SCF.
    if (!std::is_same<T, double>::value)
    {
        ModuleBase::WARNING_QUIT("ESolver_LR",
            "excited-state relaxation currently requires gamma-only sampling: the complex "
            "Z-vector solver is not implemented.");
    }
    if (this->inp_->lr_target_state >= this->nstates)
    {
        ModuleBase::WARNING_QUIT("ESolver_LR",
            "lr_target_state is beyond the states actually solved (lr_nstates <= 0 expands to "
            "all particle-hole pairs, which may be fewer than requested).");
    }

    // `openshell` is only settled after the ground-state occupations have been read, which is why
    // this cannot live in the input-file checks
    const std::string& spin = this->inp_->lr_target_spin;
    if (this->openshell)
    {
        // An open-shell calculation solves one spin-conserving channel, so there is nothing to
        // choose: whatever lr_target_spin says, this is the state that gets relaxed. Only an
        // explicit `triplet` is worth mentioning -- `singlet` is the default and expresses no
        // intent, and `updown` is already the right name for this channel.
        if (spin == "triplet")
        {
            this->ofs_running_ << " WARNING: lr_target_spin=triplet is ignored. This is an"
                " open-shell calculation with a single spin-conserving channel (updown), which is"
                " what the relaxation will follow." << std::endl;
        }
        this->target_is_ = 0;
    }
    else
    {
        // Closed shell is the opposite case: singlet and triplet are genuinely different states
        // with different gradients, so `updown` here is ambiguous rather than redundant -- it
        // usually means lr_unrestricted was meant to be set.
        if (spin == "updown")
        {
            ModuleBase::WARNING_QUIT("ESolver_LR",
                "lr_target_spin=updown, but this is a closed-shell calculation, where singlet and "
                "triplet are separate states with separate gradients. Pick one of them, or set "
                "lr_unrestricted to run spin-unrestricted.");
        }
        this->target_is_ = (spin == "triplet") ? 1 : 0;
        if (this->target_is_ >= this->nspin)
        {
            ModuleBase::WARNING_QUIT("ESolver_LR",
                "lr_target_spin=triplet requires nspin=2: the triplet channel is not built here.");
        }
    }
    this->force_gs_.create(this->ucell_->nat, 3);
    this->lr_force_.create(this->ucell_->nat, 3);
    this->target_state_ = this->inp_->lr_target_state;   // seed; overlap takes over from step 2
    this->ofs_running_ << " Excited-state relaxation follows state " << this->inp_->lr_target_state
        << " of the " << (this->openshell ? "updown" : (this->target_is_ == 1 ? "triplet" : "singlet"))
        << " channel, tracked by cross-geometry orbital and amplitude overlap." << std::endl;
}

/// Compare roots in a common electron-hole basis using exact cross-geometry AO overlap.
template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::follow_target_state_(std::ofstream& ofs)
{
    ModuleBase::timer::start("ESolver_LR", "follow_target_state_");
    if (!this->excited_relax_)
    {
        ModuleBase::timer::end("ESolver_LR", "follow_target_state_");
        return;
    }
    const int channel = this->openshell ? 0 : this->target_is_;
    const T* const amplitudes = this->X[channel].template data<T>();
    const TwoCenterIntegrator& overlap_integrator = *this->tcb().overlap_orb;
    const bool reference_first = LR_Util::tolower(this->inp_->lr_degen_mode) == "jt";
    const LR::RootInputs<T> inputs{*this->ucell_, this->orb_cutoff_, overlap_integrator,
        *this->psi_ks, this->paraC_, this->paraMat_, this->paraX_, this->nocc, this->nvirt, this->openshell, this->target_is_, reference_first};
    const int energy_offset = channel * this->nstates;
    std::vector<double> energies(this->nstates);
    std::copy_n(this->pelec->ekb.c + energy_offset, this->nstates, energies.begin());
    const bool single_state = LR_Util::tolower(this->inp_->lr_degen_mode) == "state";
    const double group_threshold = single_state ? 0.0 : this->inp_->lr_degen_thr;
    LR::follow_cross_root(inputs, amplitudes, this->nloc_per_state, this->nstates,
        this->inp_->lr_target_state, energies, group_threshold, this->target_state_, this->target_X_prev_, this->target_basis_prev_, ofs);
    ModuleBase::timer::end("ESolver_LR", "follow_target_state_");
}

template<typename T, typename TR>
double ModuleESolver::ESolver_LR<T, TR>::cal_energy()
{
    // Outside a relaxation nothing consumes this, and returning a non-zero value would change
    // what the existing single-point outputs report.
    if (!this->excited_relax_) { return 0.0; }
    // `target_omega_()` is the multiplet average under `lr_degen_mode = average` and the
    // single state otherwise, matching whatever `cal_lr_force_relax_` produced the gradient of.
    // The energy-based optimisers line-search on this, so the two must not describe different
    // surfaces.
    return this->etot_gs_ + this->target_omega_();
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::cal_force(BaseCell& basecell, ModuleBase::matrix& force)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    const UnitCell& ucell = static_cast<const UnitCell&>(basecell);
    if (!this->excited_relax_)
    {   // single-point runs print the gradients of every state from `after_all_runners` instead
        return;
    }
    if (this->lr_force_.nr != ucell.nat)
    {
        ModuleBase::WARNING_QUIT("ESolver_LR::cal_force",
            "the excited-state gradient has not been computed for this geometry.");
    }
    // Both halves already follow the ABACUS force convention F = -dE/dR (Ry/Bohr), so they add.
    // The "Gradients of each excited state" heading that `cal_force(int)` prints under is a
    // misnomer: the finite-difference reference it was validated against (`abacus-fd lr-custom`)
    // computes (E(-h) - E(+h))/h, which is -d(Omega)/dR, and the two agree in sign.
    force.create(ucell.nat, 3);
    force = this->force_gs_ + this->lr_force_;
    ModuleIO::print_force(this->ofs_running_, ucell, "EXCITED-STATE TOTAL-FORCE (eV/Angstrom)", force, false);
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::cal_stress(BaseCell& basecell, ModuleBase::matrix& stress)
{
    static_cast<void>(stress);
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    ModuleBase::WARNING_QUIT("ESolver_LR::cal_stress",
        "the excited-state stress is not implemented (every LR gradient stress term is a "
        "dummy passed with isstress=false).");
}

template<typename T, typename TR>
ct::Tensor ModuleESolver::ESolver_LR<T, TR>::pad_group_to_z_(const int ispin,
    const std::vector<int>& group) const
{
    const int d = static_cast<int>(group.size());
    const int nloc_g = this->nloc_per_state_z_;
    ct::Tensor Xz = LR_Util::newTensor<T>({ d, nloc_g });
    Xz.zero();
    // One at a time rather than as a range: `group` is sorted, but nothing guarantees its members
    // are contiguous in the state list.
    for (int k = 0; k < d; ++k)
    {
        const ct::Tensor one = this->pad_X_to_z_(ispin, group[k], 1);
        std::copy(one.template data<T>(), one.template data<T>() + nloc_g,
            Xz.template data<T>() + static_cast<size_t>(k) * nloc_g);
    }
    return Xz;
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::resolve_target_multiplet_()
{
    // Keyed on `target_state_`, not on `lr_target_state`: the overlap tracking re-chooses the
    // followed root at every ionic step, and the multiplet has to be the one containing the root
    // actually being followed. They differ as soon as two surfaces have crossed.
    this->target_group_.clear();
    if (LR_Util::tolower(this->inp_->lr_degen_mode) == "state") { return; }
    const int ekb_off = this->openshell ? 0 : this->target_is_ * this->nstates;
    std::vector<double> omega(this->nstates);
    for (int ist = 0; ist < this->nstates; ++ist) { omega[ist] = this->pelec->ekb.c[ekb_off + ist]; }
    const std::vector<std::vector<int>> groups
        = LR::group_degenerate_states(omega, this->inp_->lr_degen_thr);
    for (const std::vector<int>& g : groups)
    {
        if (std::find(g.begin(), g.end(), this->target_state_) == g.end()) { continue; }
        // a one-member group means the target is not degenerate here, and averaging over it would
        // be the single-state path with extra steps
        const bool average = LR_Util::tolower(this->inp_->lr_degen_mode) == "average";
        if (average && g.size() > 3)
        {
            this->ofs_running_ << " WARNING: tracked group has " << g.size()
                << " roots; average supports at most three. Following the selected single root." << std::endl;
        }
        else if (g.size() > 1) { this->target_group_ = g; }
        break;
    }
}

template<typename T, typename TR>
double ModuleESolver::ESolver_LR<T, TR>::target_omega_() const
{
    if (this->target_group_.empty()) { return this->pelec->ekb.c[this->target_ekb_offset_()]; }
    const int ekb_off = this->openshell ? 0 : this->target_is_ * this->nstates;
    double sum = 0.0;
    for (const int ist : this->target_group_) { sum += this->pelec->ekb.c[ekb_off + ist]; }
    return sum / static_cast<double>(this->target_group_.size());
}

template<typename T, typename TR>
ModuleBase::matrix ModuleESolver::ESolver_LR<T, TR>::cal_lr_force_relax_(std::ofstream& ofs)
{
    this->resolve_target_multiplet_();
    if (this->target_group_.empty())
    {
        return this->cal_force(this->target_is_, this->target_state_)[0];
    }
    const bool jt_mode = (LR_Util::tolower(this->inp_->lr_degen_mode) == "jt");
    // `average` needs only the DIAGONAL of the gradient matrix: the average is basis-independent by
    // construction, so the off-diagonal part (and the extra d(d-1)/2 solves it costs) is not
    // involved. The Jahn-Teller direction is orthogonal to the average and does need them.
    const int d = static_cast<int>(this->target_group_.size());
    const int ekb_off = this->openshell ? 0 : this->target_is_ * this->nstates;
    const ct::Tensor Xz = this->pad_group_to_z_(this->target_is_, this->target_group_);
    std::vector<double> omega(d);
    for (int k = 0; k < d; ++k) { omega[k] = this->pelec->ekb.c[ekb_off + this->target_group_[k]]; }
    const std::vector<ModuleBase::matrix> forces = this->openshell
        ? this->cal_force_openshell_Xz(Xz, omega, this->target_group_.front())
        : this->cal_force_Xz(this->target_is_, Xz, omega, this->target_group_.front());
    ofs << " Followed state " << this->target_state_ << " is degenerate with";
    for (const int ist : this->target_group_)
    {
        if (ist != this->target_state_) { ofs << " " << ist; }
    }
    ofs << " (Omega_bar = " << this->target_omega_() << " Ry)." << std::endl;
    if (!jt_mode)
    {
        ofs << " lr_degen_mode=average: following the multiplet average, which keeps the"
            " geometry on the current average surface; split groups are not recombined." << std::endl;
        return average_forces(forces);
    }
    return this->cal_jt_force_(forces, ofs);
}

template<typename T, typename TR>
ModuleBase::matrix ModuleESolver::ESolver_LR<T, TR>::cal_jt_force_(
    const std::vector<ModuleBase::matrix>& diag, std::ofstream& ofs)
{
    // Regime (c): descend the Jahn-Teller branch. This needs the whole gradient matrix, so the
    // off-diagonal elements are assembled here -- d(d-1)/2 further Z-vector solves on top of the
    // d diagonal ones already in `diag`.
    const int d = static_cast<int>(this->target_group_.size());
    const std::vector<std::vector<ModuleBase::matrix>> g
        = this->cal_grad_matrix_degenerate(this->target_is_, this->target_group_, diag, ofs);
    const int nat = diag[0].nr;
    const int ncoord = nat * 3;
    // flatten to the layout `find_jt_direction` takes: one d x d block per nuclear coordinate.
    // FORCES go in, not gradients, so the direction that comes back points downhill.
    std::vector<double> gflat(static_cast<size_t>(ncoord) * d * d, 0.0);
    for (int iat = 0; iat < nat; ++iat)
    {
        for (int ixyz = 0; ixyz < 3; ++ixyz)
        {
            const int a = iat * 3 + ixyz;
            for (int k = 0; k < d; ++k)
            {
                for (int l = 0; l < d; ++l)
                {
                    gflat[static_cast<size_t>(a) * d * d + k * d + l] = g[k][l](iat, ixyz);
                }
            }
        }
    }
    const LR::JTDirection jt = LR::find_jt_direction(gflat, ncoord, d);
    const int channel = this->openshell ? 0 : this->target_is_;
    const T* const amplitudes = this->X[channel].template data<T>();
    LR::save_mixed_root(amplitudes, this->nloc_per_state, this->target_group_, jt.mixing, this->target_X_prev_);
    std::vector<double> jt_part;
    const std::vector<double> sym = LR::split_symmetric_part(gflat, ncoord, d, jt.mixing, jt_part);

    const double fac = ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
    ofs << " lr_degen_mode=jt: descending the steepest branch of the multiplet." << std::endl
        << "   |F| of that branch  = " << jt.slope * fac << " eV/Angstrom" << std::endl
        << "   mixing v            =";
    for (int k = 0; k < d; ++k) { ofs << " " << jt.mixing[k]; }
    ofs << std::endl
        << "   starts agreeing     = " << jt.restarts_agreeing << " of "
        << (d + (1 << (d - 1))) << " (iterations " << jt.iterations << ")" << std::endl;
    // The two halves matter separately: the symmetric part is common to the whole multiplet and
    // only relaxes the geometry, while the remainder is what actually breaks the degeneracy. At a
    // stationary point of the average surface the first is zero and the whole force is Jahn-Teller.
    double nsym = 0.0;
    double njt = 0.0;
    for (int a = 0; a < ncoord; ++a)
    {
        nsym += sym[a] * sym[a];
        njt += jt_part[a] * jt_part[a];
    }
    ofs << "   |symmetric part|    = " << std::sqrt(nsym) * fac << " eV/Angstrom (common to the"
        " multiplet; relaxes the geometry without splitting it)" << std::endl
        << "   |Jahn-Teller part|  = " << std::sqrt(njt) * fac << " eV/Angstrom (the symmetry-"
        "breaking remainder)" << std::endl;
    if (std::sqrt(njt) < 1e-8)
    {
        ofs << " WARNING: the symmetry-breaking part vanishes -- every branch of this multiplet has"
            " the same gradient, so there is no Jahn-Teller direction to descend here. This is the"
            " expected outcome for a linear molecule, where the effect is second order"
            " (Renner-Teller) and no first-order term exists." << std::endl;
    }
    ofs << " NOTE: this is the first-order DIRECTION only. The distortion amplitude also needs the"
        " harmonic term, and the step norm is Cartesian rather than mass-weighted." << std::endl;

    // The force handed back is that of the descending branch, q(v) with v the optimal mixing --
    // i.e. the branch's own gradient, which is what a relaxation must follow. Once a step has
    // split the multiplet, `resolve_target_multiplet_` finds no group and the ordinary
    // single-state path takes over, so this mode is self-limiting by construction.
    ModuleBase::matrix f(nat, 3);
    for (int iat = 0; iat < nat; ++iat)
    {
        for (int ixyz = 0; ixyz < 3; ++ixyz)
        {
            const int a = iat * 3 + ixyz;
            f(iat, ixyz) = sym[a] + jt_part[a];
        }
    }
    return f;
}

template<typename T, typename TR>
ModuleBase::matrix ModuleESolver::ESolver_LR<T, TR>::MultipletLVC::average_force() const
{
    assert(!this->g.empty());
    std::vector<ModuleBase::matrix> diag;
    for (int k = 0; k < this->dim(); ++k) { diag.push_back(this->g[k][k]); }
    return average_forces(diag);
}

template class ModuleESolver::ESolver_LR<double, double>;
template class ModuleESolver::ESolver_LR<std::complex<double>, double>;
