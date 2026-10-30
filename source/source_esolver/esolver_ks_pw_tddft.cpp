#include "source_esolver/esolver_ks_pw_tddft.h"

#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/module_parallel/para_bridge.h"
#include "source_base/timer.h"
#include "source_estate/elecstate_pw.h"
#include "source_estate/elecstate_tools.h"
#include "source_estate/module_charge/chg_rhog_io.h"
#include "source_estate/module_charge/chg_symm.h"
#include "source_hamilt/hamilt_hs_adapter.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_ctrl/ctrl_output_pw.h"
#include "source_io/module_ctrl/ctrl_output_td_pw.h"
#include "source_io/module_efield/td_efield_io.h"
#include "source_io/module_efield/td_vecpot_io.h"
#include "source_io/module_energy/td_eij_pw.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_output/td_conserv_io.h"
#include "source_io/module_parameter/parameter.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
#include "source_pw/module_rt/td_endpoint_pw.h"
#include "source_pw/module_rt/td_pw.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

namespace ModuleESolver
{

template <typename T, typename Device>
ESolver_KS_PW_TDDFT<T, Device>::ESolver_KS_PW_TDDFT() : log_(GlobalV::ofs_running)
{
    this->classname = "ESolver_KS_PW_TDDFT";
    this->basisname = "PW";
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::before_all_runners(BaseCell& basecell, const Input_para& inp)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);
    const std::string& output_dir = PARAM.globalv.global_out_dir;

    this->td_field_manager_ = elecstate::create_td_field_manager(inp);
    if (inp.td_max_halvings > 0)
    {
        const int parts = 1 << inp.td_max_halvings;
        for (const elecstate::TDField& field: this->td_field_manager_->fields())
        {
            if (field.subdivisions() > std::numeric_limits<int>::max() / parts)
            {
                ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT", "td_max_halvings exceeds the field sampling integer range.");
            }
        }
    }
    if (inp.init_vecpot)
    {
        // An empty directory selects the calculation working directory, not OUT.<suffix>.
        const std::vector<ModuleBase::Vector3<double>> samples = ModuleIO::read_td_vecpot("");
        this->td_field_manager_->set_A_samples(samples);
    }

    ESolver_KS_PW<T, Device>::before_all_runners(ucell, inp);

    // The parent resolves the XC functional, including pseudopotential defaults.
    const int xc_type = XC_Functional::get_func_type();
    pw::check_td_input(inp, ucell, xc_type, this->general_exx_info_.cal_exx);

    current_req_.charge_total = inp.out_current == 1;
    current_req_.charge_per_k = current_req_.charge_total && inp.out_current_k;
    current_req_.spin_current = inp.nspin == 4 && inp.td_out_spincur;
    vel_.initialize(ucell, *this->pw_wfc, this->kv, this->ppcell);

    if (inp.nspin == 4)
    {
        spin_.initialize(ucell,
                         *this->pw_rhod,
                         *this->pw_wfc,
                         this->kv,
                         this->sf,
                         this->ppcell,
                         inp.td_out_mag,
                         inp.orbital_dir,
                         inp.onsite_radius,
                         inp.td_out_torque);
    }

    this->pelec->pot->set_td_field_manager(this->td_field_manager_);
#ifdef __MPI
    const hsolver::diag_comm_info comm(this->pw_wfc->pool_world, this->pw_wfc->poolrank, this->pw_wfc->poolnproc);
#else
    const hsolver::diag_comm_info comm(0, 1);
#endif
    const int npol = ucell.get_npol();
    this->td_prop_ = pw::make_td_prop<T, Device>(*this->pw_wfc, inp, npol, this->kv.ik2iktot, comm, log_);
    this->history_.prepare(*this->pelec->pot);

    // Preserve existing field history until input validation and initialization succeed.
    if (GlobalV::MY_RANK == 0)
    {
        if (inp.out_efield)
        {
            ModuleIO::prep_td_field_out(output_dir, *this->td_field_manager_, false);
        }
        if (inp.out_vecpot)
        {
            if (inp.init_vecpot)
            {
                // Reject output that aliases the input file before truncating it.
                ModuleIO::check_td_vecpot_paths("", output_dir);
            }
            ModuleIO::prep_td_vecpot_out(output_dir, false);
        }
    }
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::runner(BaseCell& basecell, const int istep)
{
    const int steps_per_md = this->inp_->estep_per_md;
    int first_step = 0;
    int last_step = 0;
    if (this->inp_->mdp.md_nstep == 1)
    {
        // Include the initial state before the fixed-ion propagation sequence.
        last_step = steps_per_md;
    }
    else if (istep > 0)
    {
        first_step = (istep - 1) * steps_per_md + 1;
        last_step = istep * steps_per_md;
    }

    for (int electronic_step = first_step; electronic_step <= last_step; ++electronic_step)
    {
        if (electronic_step > 0)
        {
            basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
            UnitCell& cell = static_cast<UnitCell&>(basecell);
            this->run_prop(cell, electronic_step);
        }
        else
        {
            ESolver_KS::runner(basecell, electronic_step);
        }
    }
}

template <typename T, typename Device>
bool ESolver_KS_PW_TDDFT<T, Device>::defer_iter_output(const int istep) const
{
    return istep > 0;
}

template <typename T, typename Device>
int ESolver_KS_PW_TDDFT<T, Device>::scf_iter_limit(const int istep) const
{
    return istep > 0 ? this->inp_->td_scf_nmax : this->inp_->scf_nmax;
}

template <typename T, typename Device>
bool ESolver_KS_PW_TDDFT<T, Device>::advance_interval(UnitCell& cell,
                                                      int step,
                                                      int index,
                                                      int depth,
                                                      const elecstate::TDFieldState& full,
                                                      pw::TDElecSnap& accepted,
                                                      hsolver::TDOrthStats& stats)
{
    const int parts = 1 << depth;
    if (depth > 0)
    {
        *this->stp.template get_psi_t<T, Device>() = this->history_.previous();
        pw::restore_td(accepted, this->pelec->f_en, this->pelec->ekb, this->chr);
        this->pelec->vnew_exist = false;
        this->td_prop_->invalidate_basis();
        this->td_prop_->orbitals().reset_orth_stats();

        const elecstate::TDFieldState trial = this->td_field_manager_->subinterval(full, index, parts);
        this->td_field_manager_->bind_field(trial);

        // Length-gauge fields are fixed potential components and must be rebuilt
        // at the trial endpoint, just as before_scf does for an original step.
        this->pelec->pot->init_pot(&this->chr);
        this->prepare_init_hamilt(cell);
    }

    const double left_time = (step - 1 + static_cast<double>(index) / parts) * full.dt * ModuleBase::AU_to_FS;
    const double right_time = (step - 1 + static_cast<double>(index + 1) / parts) * full.dt * ModuleBase::AU_to_FS;
    const std::streamsize precision = log_.precision();
    log_ << " TD trial: [" << std::setprecision(12) << left_time << ", " << right_time << "] fs, depth=" << depth << '\n';
    log_.precision(precision);

    std::ostringstream interval;
    interval << std::setprecision(12) << "[" << left_time << ", " << right_time << "] fs";
    const std::string interval_text = interval.str();
    if (depth > 0)
    {
        std::cout << " TD substep: " << interval_text << std::endl;
    }

    const int iter_limit = this->scf_iter_limit(step);
    this->td_prop_->begin_trial();
    const bool converged = this->scf_loop(cell, step, iter_limit);

    const bool accept_trial = converged || depth == this->inp_->td_max_halvings;
    log_ << " TD trial result: iterations=" << this->niter << ", residual=" << this->drho << ", converged=" << converged
         << ", accepted=" << accept_trial << '\n';
    if (!converged && depth < this->inp_->td_max_halvings)
    {
        std::cout << " TD retry: SCF not converged; splitting " << interval_text << " into two half steps." << std::endl;
        this->td_prop_->discard_trial();

        const int next_depth = depth + 1;
        const int left_index = 2 * index;
        const int right_index = left_index + 1;

        const bool left_converged = advance_interval(cell, step, left_index, next_depth, full, accepted, stats);
        const bool right_converged = advance_interval(cell, step, right_index, next_depth, full, accepted, stats);
        return left_converged && right_converged;
    }

    if (!converged)
    {
        psi::Psi<T, Device>& current = *this->stp.template get_psi_t<T, Device>();
        elecstate::ElecStatePW<T, Device>& estate = *static_cast<elecstate::ElecStatePW<T, Device>*>(this->pelec);
        hamilt::HamiltPW<T, Device>& hamilt = *static_cast<hamilt::HamiltPW<T, Device>*>(this->p_hamilt);
        const ModuleBase::Vector3<double>& vecpot = this->td_field_manager_->A_right_ha();
        pw::rebuild_td_end(current,
                           estate,
                           this->chr,
                           cell,
                           hamilt,
                           this->td_prop_->orbitals(),
                           *this->pw_wfc,
                           *this->pw_rhod,
                           vecpot,
                           this->inp_->nspin,
                           this->inp_->imp_sol,
                           this->inp_->sc_mag_switch,
                           this->inp_->dft_plus_u,
                           this->inp_->assume_isolated);
    }

    const hsolver::TDOrthStats& trial_stats = this->td_prop_->orbitals().orth_stats();
    stats.before = std::max(stats.before, trial_stats.before);
    stats.after = std::max(stats.after, trial_stats.after);

    const psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
    this->history_.save(*current, *this->pelec->pot);
    this->td_prop_->finish_trial();
    if (this->inp_->td_max_halvings > 0)
    {
        accepted = pw::capture_td(this->pelec->f_en, this->pelec->ekb, this->chr);
    }

    return converged;
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::run_prop(UnitCell& cell, int step)
{
    ModuleBase::timer::start("ESolver_KS_PW_TDDFT", "run_prop");
    const ModuleBase::TimePoint started = ModuleBase::get_time();

    pw::TDElecSnap accepted;
    if (this->inp_->td_max_halvings > 0)
    {
        accepted = pw::capture_td(this->pelec->f_en, this->pelec->ekb, this->chr);
    }

    this->before_scf(cell, step);
    ModuleBase::GlobalFunc::DONE(log_, "INIT SCF");

    const elecstate::TDFieldState full = this->td_field_manager_->field_state();
    hsolver::TDOrthStats stats;
    const bool converged = advance_interval(cell, step, 0, 0, full, accepted, stats);

    // Restore original-grid sampling metadata without altering the accepted endpoint.
    this->td_field_manager_->bind_field(full);
    this->td_prop_->orbitals().set_orth_stats(stats);
    this->chr.save_rho_before_sum_band();

    const System_para& system = PARAM.globalv;
    const bool write_restart = this->inp_->out_chg[0] != -1;
    const std::string prefix = system.global_out_dir + this->inp_->suffix;
    const bool write_tau = XC_Functional::get_ked_flag();
    const bool writer_pool = GlobalV::MY_POOL == 0;
    const Parallel::ParaWorld pw_world = Parallel::make_pw_world();
    module_charge::write_density_restart(this->chr,
                                         *this->pw_rhod,
                                         cell.GT,
                                         prefix,
                                         this->inp_->nspin,
                                         write_restart,
                                         system.gamma_only_pw,
                                         write_tau,
                                         writer_pool,
                                         pw_world,
                                         &GlobalV::ofs_warning);

    this->write_iter_bands();
    ModuleIO::ctrl_iter_pw(step, this->niter, converged, this->stp, this->kv, this->pw_wfc, *this->inp_, true, log_);

    const double duration = ModuleBase::get_duration(started, ModuleBase::get_time());
    this->record_iteration(cell, duration);

    this->output_field(step);

    this->after_scf(cell, step, converged);
    this->conv_esolver = converged;

    ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "run_prop");
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::print_step(const int istep) const
{
    const int evolution_step = istep + 1;
    std::cout << " -------------------------------------------" << std::endl;
    std::cout << " STEP OF ELECTRON EVOLVE : " << evolution_step << std::endl;
    std::cout << " -------------------------------------------" << std::endl;
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::before_scf(UnitCell& ucell, const int istep)
{
    this->print_step(istep);
    this->prep_td_step(istep);
    this->td_prop_->orbitals().reset_orth_stats();

    const bool basis_updated = ucell.cell_parameter_updated;
    ESolver_KS_PW<T, Device>::before_scf(ucell, istep);
    if (basis_updated)
    {
        this->td_prop_->invalidate_basis();
    }

    this->history_.prepare(*this->pelec->pot);
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::prepare_init_hamilt(const UnitCell& ucell)
{
    if (this->td_field_manager_->gauge() == 1)
    {
        // Both endpoints also bound any interpolated sample inside the interval.
        const double q_unshifted = pw::td_momentum_bound(*this->pw_wfc, ucell.tpiba);
        const ModuleBase::Vector3<double>& left = this->td_field_manager_->A_left_ha();
        const ModuleBase::Vector3<double>& right = this->td_field_manager_->A_right_ha();
        pw::ensure_td_vnl(ucell, q_unshifted, right, left, &this->ppcell);

        elecstate::Potential& potential = *this->pelec->pot;
        const Real* veff = potential.template get_veff_smooth_data<Real>();
        auto* hamiltonian = static_cast<hamilt::HamiltPW<T, Device>*>(this->p_hamilt);
        // A file can supply nonzero A(0), which must enter H before the initial wavefunctions are prepared.
        hamiltonian->bind_td_state(veff, this->td_field_manager_->A_right_ha());
    }
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::prep_td_step(const int istep)
{
    ModuleBase::timer::start("ESolver_KS_PW_TDDFT", "prep_td_step");

    if (this->td_field_manager_->current_step() == istep)
    {
        ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "prep_td_step");
        return;
    }

    if (this->td_field_manager_->gauge() == 0 || istep == 0)
    {
        this->td_field_manager_->prepare_sample(istep);
    }
    else
    {
        const int interval_left = istep - 1;
        this->td_field_manager_->prepare_interval(interval_left, istep);
    }

    if (this->td_field_manager_->gauge() == 1)
    {
        const ModuleBase::Vector3<double>& A_right_ha = this->td_field_manager_->A_right_ha();
        ModuleBase::GlobalFunc::OUT(log_, "Cartesian vector potential Ax(t)", A_right_ha[0]);
        ModuleBase::GlobalFunc::OUT(log_, "Cartesian vector potential Ay(t)", A_right_ha[1]);
        ModuleBase::GlobalFunc::OUT(log_, "Cartesian vector potential Az(t)", A_right_ha[2]);
    }

    if (!this->defer_iter_output(istep))
    {
        this->output_field(istep);
    }

    ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "prep_td_step");
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::output_field(int istep) const
{
    if (GlobalV::MY_RANK == 0)
    {
        if (this->inp_->out_vecpot)
        {
            ModuleIO::write_td_vecpot(PARAM.globalv.global_out_dir, istep, this->td_field_manager_->A_right_ha());
        }
        if (this->inp_->out_efield)
        {
            ModuleIO::write_td_field(*this->td_field_manager_, PARAM.globalv.global_out_dir);
        }
    }
}

template <typename T, typename Device>
std::string ESolver_KS_PW_TDDFT<T, Device>::diag_policy(const int istep) const
{
    if (istep == 0)
    {
        return "ksdft";
    }
    return ESolver_KS::diag_policy(istep);
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::iter_init(UnitCell& ucell, const int istep, const int iter)
{
    int ionic_step = 0;
    if (this->inp_->mdp.md_nstep > 1 && istep > 0)
    {
        ionic_step = (istep - 1) / this->inp_->estep_per_md + 1;
    }
    ModuleIO::write_head_td(GlobalV::ofs_running, ionic_step, istep, iter, this->basisname);
    ESolver_KS_PW<T, Device>::iter_init(ucell, istep, iter);
}

template <typename T, typename Device>
ESolver_KS::DensityStage ESolver_KS_PW_TDDFT<T, Device>::density_stage(const int istep, const int iter) const
{
    if (istep > 0 && this->td_prop_->stage(iter) == pw::PropStage::predictor)
    {
        return ESolver_KS::DensityStage::predictor;
    }
    return ESolver_KS::DensityStage::standard;
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr)
{
    ModuleBase::timer::start("ESolver_KS_PW_TDDFT", "hamilt2rho_single");

    if (istep == 0)
    {
        ESolver_KS_PW<T, Device>::hamilt2rho_single(ucell, istep, iter, ethr);

        psi::Psi<T, Device>& current = *this->stp.template get_psi_t<T, Device>();
        elecstate::ElecStatePW<T, Device>& estate = *static_cast<elecstate::ElecStatePW<T, Device>*>(this->pelec);
        hamilt::Hamilt<T, Device>* hamilt = static_cast<hamilt::Hamilt<T, Device>*>(this->p_hamilt);
        hamilt::HamiltHSOperator<T, Device> op(hamilt, this->pw_wfc);
        const bool do_orth = this->inp_->td_orthonormal != "none";
        pw::orth_init_pw(this->td_prop_->orbitals(),
                         op,
                         current,
                         estate,
                         this->chr,
                         ucell,
                         *this->pw_rhod,
                         this->inp_->nspin,
                         this->inp_->nbands,
                         iter,
                         do_orth);

        ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "hamilt2rho_single");
        return;
    }

    ModuleBase::TITLE("ESolver_KS_PW_TDDFT", "hamilt2rho_single (RT-TDDFT)");
    psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
    elecstate::Potential& potential = *this->pelec->pot;
    hamilt::HamiltPW<T, Device>* hamiltonian = static_cast<hamilt::HamiltPW<T, Device>*>(this->p_hamilt);
    hamilt::HamiltHSOperator<T, Device> op(hamiltonian, this->pw_wfc);

    const bool detailed = this->inp_->out_level == "ie";
    this->td_prop_->evolve(*hamiltonian, potential, this->history_, *this->td_field_manager_, current, istep, iter, detailed);

    elecstate::ElecStatePW<T, Device>* estate = static_cast<elecstate::ElecStatePW<T, Device>*>(this->pelec);
    // Propagation keeps the occupations of the initial state fixed.
    estate->psiToRho(*current);
    module_charge::symmetrize_rho(this->inp_->nspin, this->chr, this->pw_rhod, ucell.symm);

    // Pair the band sum with KS double-counting terms using this same potential.
    // The later KS potential update belongs to the next SCF iteration.
    this->td_prop_->orbitals().cal_band_energy(op, *current, &this->pelec->ekb);
    elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

    ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "hamilt2rho_single");
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    if (istep == 0 && this->inp_->td_orthonormal == "none")
    {
        const psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
        this->td_prop_->orbitals().check_initial(*current, this->niter);
    }
    if (!conv_esolver && istep == 0)
    {
        psi::Psi<T, Device>& current = *this->stp.template get_psi_t<T, Device>();
        elecstate::ElecStatePW<T, Device>& estate = *static_cast<elecstate::ElecStatePW<T, Device>*>(this->pelec);
        hamilt::HamiltPW<T, Device>& hamilt = *static_cast<hamilt::HamiltPW<T, Device>*>(this->p_hamilt);
        const ModuleBase::Vector3<double>& vecpot = this->td_field_manager_->A_right_ha();
        pw::rebuild_td_end(current,
                           estate,
                           this->chr,
                           ucell,
                           hamilt,
                           this->td_prop_->orbitals(),
                           *this->pw_wfc,
                           *this->pw_rhod,
                           vecpot,
                           this->inp_->nspin,
                           this->inp_->imp_sol,
                           this->inp_->sc_mag_switch,
                           this->inp_->dft_plus_u,
                           this->inp_->assume_isolated);
    }

    // Endpoint diagnostics use the accepted endpoint, not the propagation midpoint or a mixed trial density.
    hamilt::HamiltPW<T, Device>* endpoint_hamilt = static_cast<hamilt::HamiltPW<T, Device>*>(this->p_hamilt);
    const Real* endpoint_veff = this->pelec->pot->template get_veff_smooth_data<Real>();
    endpoint_hamilt->bind_td_state(endpoint_veff, this->td_field_manager_->A_right_ha());
    if (this->inp_->td_print_eij >= 0.0)
    {
        psi::Psi<T, Device>& current = *this->stp.template get_psi_t<T, Device>();
        hamilt::HamiltHSOperator<T, Device> op(endpoint_hamilt, this->pw_wfc);
        ModuleIO::print_eij_pw(op,
                               this->td_prop_->orbitals(),
                               current,
                               this->kv,
                               this->inp_->nspin,
                               istep,
                               this->niter,
                               this->inp_->td_print_eij,
                               log_);
    }

    ESolver_KS_PW<T, Device>::after_scf(ucell, istep, conv_esolver);

    const std::streamsize original_precision = std::cout.precision();
    std::cout << " Potential (Ry): " << std::setprecision(15) << this->pelec->f_en.etot << std::endl;
    std::cout.precision(original_precision);

    if (istep >= 0)
    {
        // Propagation intervals already saved their accepted endpoint before reaching this hook.
        if (istep == 0)
        {
            psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
            this->history_.save(*current, *this->pelec->pot);
        }

        const psi::Psi<T, Device>& current = *this->stp.template get_psi_t<T, Device>();
        const ModuleBase::Vector3<double>& vecpot = this->td_field_manager_->A_right_ha();
        const pw::TDVelResultPW currents = vel_.calculate(current, *this->pelec, vecpot, current_req_);

        pw::TDSpinResultPW moments;
        if (this->inp_->nspin == 4 && (this->inp_->td_out_mag != "none" || this->inp_->td_out_torque))
        {
            moments = spin_.calculate(current, *this->pelec, vecpot);
        }

        const double dt_ha = this->td_field_manager_->dt_ha();
        const int nkstot = this->kv.get_nkstot();
        ModuleIO::ctrl_output_td_pw(currents,
                                    moments,
                                    current_req_,
                                    this->inp_->td_out_mag,
                                    this->inp_->orbital_dir,
                                    this->inp_->onsite_radius,
                                    this->inp_->td_out_torque,
                                    ucell.orbital_fn,
                                    istep,
                                    dt_ha,
                                    this->inp_->nspin,
                                    nkstot,
                                    GlobalV::MY_RANK,
                                    PARAM.globalv.global_out_dir,
                                    log_);
    }

    if (this->inp_->td_out_stat)
    {
        const double local_npsi = this->td_prop_->orbitals().wave_electrons(this->pelec->wg);
        const hsolver::TDOrthStats& stats = this->td_prop_->orbitals().orth_stats();
        const bool full_gram = this->inp_->td_orthonormal != "none";
        const pw::TDConserv values
            = pw::cal_td_conserv(local_npsi, this->chr, ucell.omega, this->pw_wfc->poolrank, stats.before, stats.after, full_gram);
        if (istep == 0)
        {
            initial_wave_electrons_ = values.npsi;
            initial_rho_electrons_ = values.nrho;
        }
        ModuleIO::print_conserv(values,
                                initial_wave_electrons_,
                                initial_rho_electrons_,
                                istep,
                                full_gram,
                                this->inp_->td_orthonormal,
                                log_);
    }
}

template class ESolver_KS_PW_TDDFT<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_KS_PW_TDDFT<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class ESolver_KS_PW_TDDFT<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_KS_PW_TDDFT<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace ModuleESolver
