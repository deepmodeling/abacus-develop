#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_esolver/esolver_lr_lcao_tddft.h"
#include "source_io/module_output/output_log.h"
#include "source_lcao/module_lr/cal_edm.h"
#include "source_lcao/module_lr/grad_degen.h"
#include "source_lcao/module_lr/gradient_inputs.h"
#include "source_lcao/module_lr/gradient_output.h"
#include "source_lcao/module_lr/lr_amp.h"
#include "source_lcao/module_lr/lr_force.h"
#include "source_lcao/module_lr/zeq_solver.h"
#include <algorithm>
#include <complex>
#include <numeric>

using namespace LR;

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::setup_relax_target_()
{
    this->target_.setup(*this->inp_, this->openshell, this->nspin, this->nstates, this->ofs_running_);
    if (this->target_.active())
    {
        this->force_gs_.create(this->ucell_->nat, 3);
        this->lr_force_.create(this->ucell_->nat, 3);
    }
}

/// Compare roots in a common electron-hole basis using exact cross-geometry AO overlap.
template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::follow_target_state_(std::ofstream &ofs)
{
    ModuleBase::timer::start("ESolver_LR", "follow_target_state_");
    if (!this->target_.active())
    {
        ModuleBase::timer::end("ESolver_LR", "follow_target_state_");
        return;
    }
    const int channel = this->openshell ? 0 : this->target_.ispin();
    const T *const amplitudes = this->X[channel].template data<T>();
    const TwoCenterIntegrator &overlap_integrator = *this->tcb().overlap_orb;
    std::vector<double> energies(this->nstates);
    std::copy_n(this->pelec->ekb.c + this->target_.ekb_off(this->nstates), this->nstates, energies.begin());
    this->target_.follow(*this->inp_, *this->ucell_, this->orb_cutoff_, overlap_integrator, *this->psi_ks, this->paraC_,
                         this->paraMat_, this->paraX_, this->nocc, this->nvirt, amplitudes, this->nloc_per_state,
                         this->nstates, energies, ofs);
    ModuleBase::timer::end("ESolver_LR", "follow_target_state_");
}

template <typename T, typename TR>
double ModuleESolver::ESolver_LR<T, TR>::cal_energy()
{
    // Outside a relaxation nothing consumes this, and returning a non-zero value would change
    // what the existing single-point outputs report.
    if (!this->target_.active())
    {
        return 0.0;
    }
    // `omega` is the multiplet average under `lr_degen_mode = average` and the single state
    // otherwise, matching whatever `cal_lr_force_relax_` produced the gradient of. The
    // energy-based optimisers line-search on this, so the two must not describe different
    // surfaces.
    return this->etot_gs_ + this->target_.omega(this->pelec->ekb.c, this->nstates);
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::cal_force(BaseCell &basecell, ModuleBase::matrix &force)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    const UnitCell &ucell = static_cast<const UnitCell &>(basecell);
    if (!this->target_.active())
    { // single-point runs print the gradients of every state from `after_all_runners` instead
        return;
    }
    if (this->lr_force_.nr != ucell.nat)
    {
        ModuleBase::WARNING_QUIT(
            "ESolver_LR::cal_force", "the excited-state gradient has not been computed for this geometry.");
    }
    // Both halves already follow the ABACUS force convention F = -dE/dR (Ry/Bohr), so they add.
    // The "Gradients of each excited state" heading that `cal_force(int)` prints under is a
    // misnomer: the finite-difference reference it was validated against (`abacus-fd lr-custom`)
    // computes (E(-h) - E(+h))/h, which is -d(Omega)/dR, and the two agree in sign.
    force.create(ucell.nat, 3);
    force = this->force_gs_ + this->lr_force_;
    ModuleIO::print_force(this->ofs_running_, ucell, "EXCITED-STATE TOTAL-FORCE (eV/Angstrom)", force, false);
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::cal_stress(BaseCell &basecell, ModuleBase::matrix &stress)
{
    static_cast<void>(stress);
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    ModuleBase::WARNING_QUIT(
        "ESolver_LR::cal_stress", "the excited-state stress is not implemented (every LR gradient stress term is a "
                                  "dummy passed with isstress=false).");
}

template <typename T, typename TR>
ct::Tensor ModuleESolver::ESolver_LR<T, TR>::pad_group_to_z_(const int ispin, const std::vector<int> &group) const
{
    const int d = static_cast<int>(group.size());
    const int nloc_g = this->nloc_per_state_z_;
    ct::Tensor Xz = LR_Util::newTensor<T>({d, nloc_g});
    Xz.zero();
    // One at a time rather than as a range: `group` is sorted, but nothing guarantees its members
    // are contiguous in the state list.
    for (int k = 0; k < d; ++k)
    {
        const ct::Tensor one = this->pad_X_to_z_(ispin, group[k], 1);
        std::copy(
            one.template data<T>(), one.template data<T>() + nloc_g,
            Xz.template data<T>() + static_cast<size_t>(k) * nloc_g);
    }
    return Xz;
}

template <typename T, typename TR>
ModuleBase::matrix ModuleESolver::ESolver_LR<T, TR>::cal_lr_force_relax_(std::ofstream &ofs)
{
    this->target_.resolve_group(*this->inp_, this->pelec->ekb.c, this->nstates, this->ofs_running_);
    const std::vector<int> &group = this->target_.group();
    if (group.empty())
    {
        return this->cal_force(this->target_.ispin(), this->target_.state())[0];
    }
    const bool jt_mode = (LR_Util::tolower(this->inp_->lr_degen_mode) == "jt");
    // `average` needs only the DIAGONAL of the gradient matrix: the average is basis-independent by
    // construction, so the off-diagonal part (and the extra d(d-1)/2 solves it costs) is not
    // involved. The Jahn-Teller direction is orthogonal to the average and does need them.
    const int d = static_cast<int>(group.size());
    const int ekb_off = this->target_.ekb_off(this->nstates);
    const ct::Tensor Xz = this->pad_group_to_z_(this->target_.ispin(), group);
    std::vector<double> omega(d);
    for (int k = 0; k < d; ++k)
    {
        omega[k] = this->pelec->ekb.c[ekb_off + group[k]];
    }
    const std::vector<ModuleBase::matrix> forces =
        this->openshell ? this->cal_force_openshell_Xz(Xz, omega, group.front())
                        : this->cal_force_Xz(this->target_.ispin(), Xz, omega, group.front());
    ofs << " Followed state " << this->target_.state() << " is degenerate with";
    for (const int ist : group)
    {
        if (ist != this->target_.state())
        {
            ofs << " " << ist;
        }
    }
    ofs << " (Omega_bar = " << this->target_.omega(this->pelec->ekb.c, this->nstates) << " Ry)." << std::endl;
    if (!jt_mode)
    {
        ofs << " lr_degen_mode=average: following the multiplet average, which keeps the"
               " geometry on the current average surface; split groups are not recombined."
            << std::endl;
        return average_forces(forces);
    }
    return this->cal_jt_force_(forces, ofs);
}

template <typename T, typename TR>
ModuleBase::matrix ModuleESolver::ESolver_LR<T, TR>::cal_jt_force_(
    const std::vector<ModuleBase::matrix> &diag, std::ofstream &ofs)
{
    // Regime (c): descend the Jahn-Teller branch. This needs the whole gradient matrix, so the
    // off-diagonal elements are assembled here -- d(d-1)/2 further Z-vector solves on top of the
    // d diagonal ones already in `diag`.
    const std::vector<int> &group = this->target_.group();
    const int d = static_cast<int>(group.size());
    const std::vector<std::vector<ModuleBase::matrix>> g =
        this->cal_grad_matrix_degenerate(this->target_.ispin(), group, diag, ofs);
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
    const int channel = this->openshell ? 0 : this->target_.ispin();
    const T *const amplitudes = this->X[channel].template data<T>();
    this->target_.save_mixed(amplitudes, this->nloc_per_state, jt.mixing);
    std::vector<double> jt_part;
    const std::vector<double> sym = LR::split_symmetric_part(gflat, ncoord, d, jt.mixing, jt_part);

    const double fac = ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
    ofs << " lr_degen_mode=jt: descending the steepest branch of the multiplet." << std::endl
        << "   |F| of that branch  = " << jt.slope * fac << " eV/Angstrom" << std::endl
        << "   mixing v            =";
    for (int k = 0; k < d; ++k)
    {
        ofs << " " << jt.mixing[k];
    }
    ofs << std::endl
        << "   starts agreeing     = " << jt.restarts_agreeing << " of " << (d + (1 << (d - 1))) << " (iterations "
        << jt.iterations << ")" << std::endl;
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
    ofs << "   |symmetric part|    = " << std::sqrt(nsym) * fac
        << " eV/Angstrom (common to the"
           " multiplet; relaxes the geometry without splitting it)"
        << std::endl
        << "   |Jahn-Teller part|  = " << std::sqrt(njt) * fac
        << " eV/Angstrom (the symmetry-"
           "breaking remainder)"
        << std::endl;
    if (std::sqrt(njt) < 1e-8)
    {
        ofs << " WARNING: the symmetry-breaking part vanishes -- every branch of this multiplet has"
               " the same gradient, so there is no Jahn-Teller direction to descend here. This is the"
               " expected outcome for a linear molecule, where the effect is second order"
               " (Renner-Teller) and no first-order term exists."
            << std::endl;
    }
    ofs << " NOTE: this is the first-order DIRECTION only. The distortion amplitude also needs the"
           " harmonic term, and the step norm is Cartesian rather than mass-weighted."
        << std::endl;

    // The force handed back is that of the descending branch, q(v) with v the optimal mixing --
    // i.e. the branch's own gradient, which is what a relaxation must follow. Once a step has
    // split the multiplet, `resolve_group` finds no group and the ordinary single-state path
    // takes over, so this mode is self-limiting by construction.
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

template <typename T, typename TR>
ModuleBase::matrix ModuleESolver::ESolver_LR<T, TR>::MultipletLVC::average_force() const
{
    assert(!this->g.empty());
    std::vector<ModuleBase::matrix> diag;
    for (int k = 0; k < this->dim(); ++k)
    {
        diag.push_back(this->g[k][k]);
    }
    return average_forces(diag);
}

template class ModuleESolver::ESolver_LR<double, double>;
template class ModuleESolver::ESolver_LR<std::complex<double>, double>;
