#include "source_esolver/esolver_ks_pw_tddft.h"

#include "source_base/global_variable.h"
#include "source_base/parallel_comm.h"
#include "source_base/timer.h"
#include "source_estate/elecstate_pw.h"
#include "source_estate/elecstate_tools.h"
#include "source_estate/module_charge/chg_symm.h"
#include "source_estate/module_pot/h_tddft_pw.h"
#include "source_hamilt/hamilt_hs_adapter.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_efield/td_efield_io.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/parameter.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
#include "source_pw/module_pwdft/td_pw.h"

namespace ModuleESolver
{

template <typename T, typename Device>
ESolver_KS_PW_TDDFT<T, Device>::ESolver_KS_PW_TDDFT()
{
    this->classname = "ESolver_KS_PW_TDDFT";
    this->basisname = "PW";
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::before_all_runners(BaseCell& basecell, const Input_para& inp)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    this->td_field_manager_ = elecstate::create_td_field_manager(inp);
    elecstate::H_TDDFT_pw::set_field_state(*this->td_field_manager_);

    ESolver_KS_PW<T, Device>::before_all_runners(ucell, inp);
    // The parent resolves the XC functional, including pseudopotential defaults.
    pw::check_td_input(inp, ucell, XC_Functional::get_ked_flag(), this->general_exx_info_.cal_exx);
    this->pelec->pot->set_td_field_manager(this->td_field_manager_);
#ifdef __MPI
    const hsolver::diag_comm_info comm(this->pw_wfc->pool_world, this->pw_wfc->poolrank, this->pw_wfc->poolnproc);
#else
    const hsolver::diag_comm_info comm(0, 1);
#endif
    this->td_solver_.reset(new hsolver::HSolverPWTDDFT<T, Device>(*this->pw_wfc,
                                                                  inp.lin_solver,
                                                                  inp.lin_precond,
                                                                  inp.lin_thr,
                                                                  inp.lin_maxiter,
                                                                  inp.t_in_h,
                                                                  comm,
                                                                  GlobalV::ofs_running));
    this->history_.prepare(*this->pelec->pot, XC_Functional::get_ked_flag());
    // Preserve existing field history until input validation and initialization succeed.
    if (inp.out_efield && GlobalV::MY_RANK == 0)
    {
        ModuleIO::prepare_td_field_output(PARAM.globalv.global_out_dir, this->td_field_manager_->fields().size(), false);
    }
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::before_scf(UnitCell& ucell, const int istep)
{
    this->prepare_td_step(istep);
    ESolver_KS_PW<T, Device>::before_scf(ucell, istep);
    this->history_.prepare(*this->pelec->pot, XC_Functional::get_ked_flag());
    if (this->td_field_manager_->gauge() == 1)
    {
        // Refresh after the parent updates the cell and distributed basis.
        this->q_unshifted_ = pw::td_momentum_bound(*this->pw_wfc, ucell.tpiba);
        pw::ensure_td_vnl(ucell,
                          this->q_unshifted_,
                          this->td_field_manager_->A_right_ha(),
                          this->td_field_manager_->A_prop_ha(),
                          &this->ppcell);
    }
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::prepare_td_step(const int istep)
{
    ModuleBase::timer::start("ESolver_KS_PW_TDDFT", "prepare_td_step");
    if (this->td_field_manager_->current_step() == istep)
    {
        ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "prepare_td_step");
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
    elecstate::H_TDDFT_pw::set_field_state(*this->td_field_manager_);
    if (this->inp_->out_efield && GlobalV::MY_RANK == 0)
    {
        ModuleIO::write_td_field_values(*this->td_field_manager_, PARAM.globalv.global_out_dir);
    }
    ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "prepare_td_step");
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
    // With estep_per_md=1, the ionic and electronic step indices coincide.
    ModuleIO::write_head_td(GlobalV::ofs_running, istep, istep, iter, this->basisname);
    ESolver_KS_PW<T, Device>::iter_init(ucell, istep, iter);
}

template <typename T, typename Device>
ESolver_KS::DensityStage ESolver_KS_PW_TDDFT<T, Device>::density_stage(const int istep, const int iter) const
{
    if (istep > 0 && iter == 1)
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
        ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "hamilt2rho_single");
        return;
    }

    ModuleBase::TITLE("ESolver_KS_PW_TDDFT", "hamilt2rho_single (RT-TDDFT)");
    psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
    elecstate::Potential& potential = *this->pelec->pot;
    hamilt::HamiltPW<T, Device>* hamiltonian = static_cast<hamilt::HamiltPW<T, Device>*>(this->p_hamilt);
    hamilt::HamiltHSOperator<T, Device> op(hamiltonian, this->pw_wfc);

    // Predict with the current potential, then correct with the midpoint.
    const Real* endpoint_veff = potential.template get_veff_smooth_data<Real>();
    const Real* endpoint_vofk = potential.template get_vofk_smooth_data<Real>();
    const ESolver_KS::DensityStage stage = this->density_stage(istep, iter);
    const bool corrector = stage != ESolver_KS::DensityStage::predictor;
    const pw::TDPotentialView<Real> propagation = this->history_.prop_pot(potential, XC_Functional::get_ked_flag(), corrector);
    ModuleBase::Vector3<double> momentum_shift(0.0, 0.0, 0.0);
    if (this->td_field_manager_->gauge() == 1)
    {
        momentum_shift = this->td_field_manager_->A_prop_ha();
    }
    hamiltonian->bind_td_state(propagation.veff, propagation.vofk, momentum_shift);
    this->td_solver_->solve(op,
                            this->history_.previous(),
                            current,
                            this->td_field_manager_->dt_ha(),
                            momentum_shift,
                            istep,
                            iter,
                            this->inp_->out_level == "ie",
                            GlobalV::ofs_running);

    // Restore the endpoint Hamiltonian before evaluating density and energy.
    elecstate::H_TDDFT_pw::set_field_state(*this->td_field_manager_);
    hamiltonian->bind_td_state(endpoint_veff, endpoint_vofk, this->td_field_manager_->A_right_ha());
    elecstate::ElecStatePW<T, Device>* estate = static_cast<elecstate::ElecStatePW<T, Device>*>(this->pelec);
    // Propagation keeps the occupations of the initial state fixed.
    estate->psiToRho(*current);
    module_charge::symmetrize_rho(this->inp_->nspin, this->chr, this->pw_rhod, ucell.symm);

    // Pair the band sum with KS double-counting terms using this same potential.
    // The later KS potential update belongs to the next SCF iteration.
    this->td_solver_->cal_band_energy(op, *current, &this->pelec->ekb);
    elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);
    ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "hamilt2rho_single");
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver)
{
    ESolver_KS_PW<T, Device>::iter_finish(ucell, istep, iter, conv_esolver);
    if (istep > 0 && !conv_esolver && iter >= this->inp_->scf_nmax)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT", "The time-dependent SCF step did not converge.");
    }
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    if (!conv_esolver)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT", "Cannot propagate an unconverged electronic state.");
    }
    ESolver_KS_PW<T, Device>::after_scf(ucell, istep, conv_esolver);
    if (istep >= 0)
    {
        psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
        this->history_.save(*current, *this->pelec->pot, XC_Functional::get_ked_flag());
        if (istep == 0)
        {
            std::cout << "[RT-TDDFT] Ground state SCF finished. Historical wavefunction and V_eff initialized." << std::endl;
        }
        if (this->inp_->out_current == 1)
        {
            const ModuleBase::Vector3<double>& A_right_ha = this->td_field_manager_->A_right_ha();
            this->current_output_.write(istep,
                                        ucell,
                                        this->pw_wfc,
                                        current,
                                        this->pelec,
                                        this->kv,
                                        &this->ppcell,
                                        this->td_field_manager_->gauge(),
                                        A_right_ha,
                                        this->inp_->out_current_k,
                                        PARAM.globalv.global_out_dir,
                                        GlobalV::MY_RANK);
        }
    }
}

template class ESolver_KS_PW_TDDFT<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_KS_PW_TDDFT<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class ESolver_KS_PW_TDDFT<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_KS_PW_TDDFT<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace ModuleESolver
