#include "source_esolver/esolver_ks_pw_tddft.h"

#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/module_parallel/para_bridge.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_device.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_estate/elecstate_pw.h"
#include "source_estate/elecstate_tools.h"
#include "source_estate/module_charge/chg_rhog_io.h"
#include "source_estate/module_charge/chg_symm.h"
#include "source_hamilt/hamilt_hs_adapter.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_ctrl/ctrl_output_pw.h"
#include "source_io/module_current/td_current_pw.h"
#include "source_io/module_current/td_spin_pw_io.h"
#include "source_io/module_efield/td_efield_io.h"
#include "source_io/module_efield/td_vecpot_io.h"
#include "source_io/module_energy/td_eij_io.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/parameter.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
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

    spin_options_.method = inp.td_out_mag;
    spin_options_.orbital_dir = inp.orbital_dir;
    spin_options_.radius = inp.onsite_radius;
    spin_options_.torque = inp.td_out_torque;

    vel_request_.charge_total = inp.out_current == 1;
    vel_request_.charge_per_k = vel_request_.charge_total && inp.out_current_k;
    vel_request_.spin_current = inp.nspin == 4 && inp.td_out_spincur;
    vel_.initialize(ucell, *this->pw_wfc, this->kv, this->ppcell);

    if (inp.nspin == 4)
    {
        spin_.initialize(ucell, *this->pw_rhod, *this->pw_wfc, this->kv, this->sf, this->ppcell, spin_options_);
    }

    this->pelec->pot->set_td_field_manager(this->td_field_manager_);
#ifdef __MPI
    const hsolver::diag_comm_info comm(this->pw_wfc->pool_world, this->pw_wfc->poolrank, this->pw_wfc->poolnproc);
#else
    const hsolver::diag_comm_info comm(0, 1);
#endif
    const int npol = ucell.get_npol();
    const pw::TDPropOpts options = pw::td_prop_options(inp, npol, this->kv.ik2iktot);
    this->td_prop_.reset(new pw::TDPropPW<T, Device>(*this->pw_wfc, options, comm, log_));
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
            this->run_propagation(cell, electronic_step);
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
typename ESolver_KS_PW_TDDFT<T, Device>::AcceptedState ESolver_KS_PW_TDDFT<T, Device>::capture_state() const
{
    AcceptedState state;
    state.energy = this->pelec->f_en;
    state.eigenvalues = this->pelec->ekb;

    state.density.resize(this->chr.nspin * this->chr.nrxx);
    for (int spin = 0; spin < this->chr.nspin; ++spin)
    {
        std::copy(this->chr.rho[spin], this->chr.rho[spin] + this->chr.nrxx, state.density.begin() + spin * this->chr.nrxx);
    }

    return state;
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::restore_state(const AcceptedState& state)
{
    ModuleBase::timer::start("ESolver_KS_PW_TDDFT", "restore_state");

    *this->stp.template get_psi_t<T, Device>() = this->history_.previous();
    this->pelec->f_en = state.energy;
    this->pelec->ekb = state.eigenvalues;
    this->pelec->vnew_exist = false;

    for (int spin = 0; spin < this->chr.nspin; ++spin)
    {
        const double* density = state.density.data() + spin * this->chr.nrxx;
        std::copy(density, density + this->chr.nrxx, this->chr.rho[spin]);
    }
    this->chr.save_rho_before_sum_band();

    // The caller binds the trial field and rebuilds the endpoint potential. The saved
    // left-endpoint potential used by the corrector is already held in history_.
    this->td_prop_->invalidate_basis();
    this->td_prop_->orbitals().reset_orth_stats();

    ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "restore_state");
}

template <typename T, typename Device>
bool ESolver_KS_PW_TDDFT<T, Device>::advance_interval(UnitCell& cell,
                                                      int step,
                                                      int index,
                                                      int depth,
                                                      const elecstate::TDFieldState& full,
                                                      AcceptedState& accepted,
                                                      hsolver::TDOrthStats& stats)
{
    const int parts = 1 << depth;
    if (depth > 0)
    {
        this->restore_state(accepted);

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
        this->finalize_endpoint(cell);
    }

    const hsolver::TDOrthStats& trial_stats = this->td_prop_->orbitals().orth_stats();
    stats.before = std::max(stats.before, trial_stats.before);
    stats.after = std::max(stats.after, trial_stats.after);

    const psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
    this->history_.save(*current, *this->pelec->pot);
    this->td_prop_->finish_trial();
    if (this->inp_->td_max_halvings > 0)
    {
        accepted = this->capture_state();
    }

    return converged;
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::run_propagation(UnitCell& cell, int step)
{
    ModuleBase::timer::start("ESolver_KS_PW_TDDFT", "run_propagation");
    const ModuleBase::TimePoint started = ModuleBase::get_time();

    AcceptedState accepted;
    if (this->inp_->td_max_halvings > 0)
    {
        accepted = this->capture_state();
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

    const auto& system = PARAM.globalv;
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

    ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "run_propagation");
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
        const auto& left = this->td_field_manager_->A_left_ha();
        const auto& right = this->td_field_manager_->A_right_ha();
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

        psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
        // Padding is a storage invariant even with orthogonalization disabled; clearing it does not change the density.
        this->td_prop_->orbitals().zero_padding(current);

        if (this->inp_->td_orthonormal != "none")
        {
            int changed = this->td_prop_->orbitals().orth_initial(current, iter);
            // Occupations and density construction contain collectives across k-point pools.
            Parallel_Reduce::reduce_all(changed);
            if (changed > 0)
            {
                hamilt::Hamilt<T, Device>* hamiltonian = static_cast<hamilt::Hamilt<T, Device>*>(this->p_hamilt);
                hamilt::HamiltHSOperator<T, Device> op(hamiltonian, this->pw_wfc);
                this->td_prop_->orbitals().cal_band_energy(op, *current, &this->pelec->ekb);

                elecstate::calculate_weights(this->pelec->ekb,
                                             this->pelec->wg,
                                             this->pelec->klist,
                                             this->pelec->eferm,
                                             this->pelec->f_en,
                                             this->pelec->nelec_spin,
                                             this->inp_->nbands,
                                             this->pelec->skip_weights);
                elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

                elecstate::ElecStatePW<T, Device>* estate = static_cast<elecstate::ElecStatePW<T, Device>*>(this->pelec);
                estate->psiToRho(*current);
                module_charge::symmetrize_rho(this->inp_->nspin, this->chr, this->pw_rhod, ucell.symm);
            }
        }

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
        finalize_endpoint(ucell);
    }

    // Endpoint diagnostics use the accepted endpoint, not the propagation midpoint or a mixed trial density.
    auto* endpoint_hamilt = static_cast<hamilt::HamiltPW<T, Device>*>(this->p_hamilt);
    const Real* endpoint_veff = this->pelec->pot->template get_veff_smooth_data<Real>();
    endpoint_hamilt->bind_td_state(endpoint_veff, this->td_field_manager_->A_right_ha());
    if (this->inp_->td_print_eij >= 0.0)
    {
        print_eij(istep);
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

        sample_endpoint(ucell, istep);
    }

    if (this->inp_->td_out_stat)
    {
        report_orth(ucell, istep);
    }
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::finalize_endpoint(UnitCell& cell)
{
    ModuleBase::timer::start("ESolver_KS_PW_TDDFT", "finalize_endpoint");

    psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
    auto* estate = static_cast<elecstate::ElecStatePW<T, Device>*>(this->pelec);
    // The last SCF iteration may have mixed the density. Keep the actual orbitals and their occupations.
    estate->psiToRho(*current);
    module_charge::symmetrize_rho(this->inp_->nspin, this->chr, this->pw_rhod, cell.symm);
    this->chr.save_rho_before_sum_band();
    this->pelec->pot->update_from_charge(&this->chr, &cell);

    auto* hamiltonian = static_cast<hamilt::HamiltPW<T, Device>*>(this->p_hamilt);
    const Real* veff = this->pelec->pot->template get_veff_smooth_data<Real>();
    hamiltonian->bind_td_state(veff, this->td_field_manager_->A_right_ha());

    hamilt::HamiltHSOperator<T, Device> op(hamiltonian, this->pw_wfc);
    this->td_prop_->orbitals().cal_band_energy(op, *current, &this->pelec->ekb);
    elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);
    this->pelec->f_en.deband = this->pelec->cal_delta_eband(cell);
    this->pelec->f_en.descf = 0.0;

    // Do not reuse a self-consistency force correction from an earlier endpoint.
    this->pelec->vnew_exist = false;
    this->pelec->cal_energies(2, this->inp_->imp_sol, this->inp_->sc_mag_switch, this->inp_->dft_plus_u, this->inp_->assume_isolated);

    ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "finalize_endpoint");
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::sample_endpoint(UnitCell& cell, int step)
{
    psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
    const ModuleBase::Vector3<double>& a = this->td_field_manager_->A_right_ha();

    // All endpoint currents share one velocity action; disabled requests perform no collectives.
    const pw::TDVelResultPW currents = vel_.calculate(*current, *this->pelec, a, vel_request_);

    pw::TDSpinResultPW moments;
    if (this->inp_->nspin == 4 && (spin_options_.method != "none" || spin_options_.torque))
    {
        moments = spin_.calculate(*current, *this->pelec, a);
    }

    const std::string& directory = PARAM.globalv.global_out_dir;
    const int rank = GlobalV::MY_RANK;
    if (rank == 0)
    {
        if (this->inp_->out_current == 1)
        {
            const int spin_mult = this->inp_->nspin == 2 ? 2 : 1;
            const int nk = this->kv.get_nkstot() / spin_mult;
            ModuleIO::write_pw_current(currents, step, nk, this->inp_->out_current_k, directory);
        }
        if (this->inp_->nspin == 4)
        {
            const double time = step * this->td_field_manager_->dt_ha() * ModuleBase::AU_to_FS;
            ModuleIO::write_td_spin_pw(directory,
                                       step,
                                       time,
                                       cell.orbital_fn,
                                       spin_options_,
                                       vel_request_.spin_current,
                                       moments,
                                       currents,
                                       log_);
        }
    }
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::print_eij(const int istep)
{
    ModuleBase::timer::start("ESolver_KS_PW_TDDFT", "print_eij");

    psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
    const int original_k = current->get_current_k();
    const int nband = current->get_nbands();
    const size_t elements = static_cast<size_t>(nband) * nband;
    const int nspin = this->inp_->nspin;
    const int nk = this->kv.get_nkstot() / (nspin == 2 ? 2 : 1);
    const Parallel_Kpoints& pools = this->kv.para_k;
    const bool main_rank = pools.my_pool == 0 && pools.rank_in_pool == 0;

    auto* hamiltonian = static_cast<hamilt::HamiltPW<T, Device>*>(this->p_hamilt);
    const Real* veff = this->pelec->pot->template get_veff_smooth_data<Real>();
    // cal_converged() updates the potential from the final density, after the last band-energy evaluation.
    hamiltonian->bind_td_state(veff, this->td_field_manager_->A_right_ha());
    hamilt::HamiltHSOperator<T, Device> op(hamiltonian, this->pw_wfc);

    for (int global_k = 0; global_k < this->kv.get_nkstot(); ++global_k)
    {
        const int physical_k = global_k % nk;
        const int spin = global_k / nk;
        int owner = 0;
        int local_k = global_k;
#ifdef __MPI
        owner = pools.whichpool[physical_k];
        local_k = physical_k - pools.startk_pool[owner] + spin * pools.nks_pool[owner];
#endif

        std::vector<std::complex<double>> eij;
        if (pools.my_pool == owner)
        {
            eij = this->td_prop_->orbitals().cal_eij(op, *current, local_k);
        }

#ifdef __MPI
        // Each pool already has the full matrix; transfer only one copy to the main log.
        if (owner != 0 && pools.rank_in_pool == 0 && (main_rank || pools.my_pool == owner))
        {
            if (main_rank)
            {
                eij.resize(elements);
            }
            MPI_Comm world = MPI_COMM_WORLD;
            const int source = pools.get_startpro_pool(owner);
            const size_t max_count = std::numeric_limits<int>::max();
            for (size_t offset = 0; offset < elements; offset += max_count)
            {
                const int count = static_cast<int>(std::min(max_count, elements - offset));
                if (main_rank)
                {
                    Parallel_Common::recv_data(eij.data() + offset, count, source, 0, world, MPI_STATUS_IGNORE);
                }
                else
                {
                    Parallel_Common::send_data(eij.data() + offset, count, 0, 0, world);
                }
            }
        }
#endif

        if (main_rank)
        {
            const ModuleIO::EijOutputInfo info{this->inp_->td_print_eij, istep + 1, this->niter, physical_k + 1, nspin, spin + 1};
            ModuleIO::print_eij(eij.data(), nband, nband, info, log_);
        }
    }

    current->fix_k(original_k);
    op.update_k(original_k);

    ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "print_eij");
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::report_orth(const UnitCell& ucell, int istep)
{
    ModuleBase::timer::start("ESolver_KS_PW_TDDFT", "report_orth");

    double electrons = this->td_prop_->orbitals().wave_electrons(this->pelec->wg);
    // Norms are replicated within a pool; only its root contributes to the total.
    if (this->pw_wfc->poolrank != 0)
    {
        electrons = 0.0;
    }
    Parallel_Reduce::reduce_all(electrons);

    double rho_electrons = 0.0;
    const int density_spins = this->chr.nspin == 2 ? 2 : 1;
    for (int spin = 0; spin < density_spins; ++spin)
    {
        for (int ir = 0; ir < this->chr.nrxx; ++ir)
        {
            rho_electrons += this->chr.rho[spin][ir];
        }
    }
    rho_electrons *= ucell.omega / this->chr.rhopw->nxyz;
    Parallel_Reduce::reduce_pool(rho_electrons);

    if (istep == 0)
    {
        initial_wave_electrons_ = electrons;
        initial_rho_electrons_ = rho_electrons;
    }

    const hsolver::TDOrthStats& stats = this->td_prop_->orbitals().orth_stats();
    const bool full_gram = this->inp_->td_orthonormal != "none";
    double maxima[] = {stats.before, stats.after};
    if (full_gram)
    {
        Parallel_Reduce::reduce_max(maxima, 2);
    }

    if (!log_.good())
    {
        ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "report_orth");
        return;
    }

    const char* stage = istep == 0 ? "initial state" : "propagation";
    const int evolution_step = istep + 1;
    const double delta_wave = electrons - initial_wave_electrons_;
    const double delta_rho = rho_electrons - initial_rho_electrons_;

    const std::ios::fmtflags saved_flags = log_.flags();
    const std::streamsize saved_precision = log_.precision();
    log_ << " PW RT-TDDFT conservation: evolution step " << evolution_step << ", " << stage << ", " << this->inp_->td_orthonormal << '\n';
    log_ << std::right << std::fixed << std::setprecision(10);
    ModuleBase::GlobalFunc::OUT(log_, "Npsi", electrons);
    ModuleBase::GlobalFunc::OUT(log_, "Nrho", rho_electrons);

    log_ << std::scientific << std::setprecision(6);
    ModuleBase::GlobalFunc::OUT(log_, "dNpsi", delta_wave);
    ModuleBase::GlobalFunc::OUT(log_, "dNrho", delta_rho);

    if (full_gram)
    {
        if (istep == 0)
        {
            ModuleBase::GlobalFunc::OUT(log_, "initial_orth_error", maxima[1]);
        }
        else
        {
            ModuleBase::GlobalFunc::OUT(log_, "orth_before_step_max", maxima[0]);
            ModuleBase::GlobalFunc::OUT(log_, "orth_after_step_max", maxima[1]);
        }
    }

    log_.flags(saved_flags);
    log_.precision(saved_precision);

    ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "report_orth");
}

template class ESolver_KS_PW_TDDFT<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_KS_PW_TDDFT<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class ESolver_KS_PW_TDDFT<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_KS_PW_TDDFT<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace ModuleESolver
