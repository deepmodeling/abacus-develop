#include "esolver_ks_pw.h"

#include "source_estate/cal_ux.h"
#include "source_estate/elecstate_pw.h"
#include "source_estate/module_charge/symmetry_rho.h"

#include "source_hsolver/diago_iter_assist.h"
#include "source_hsolver/hsolver_pw.h"
#include "source_hsolver/diago_params.h"

#include "source_hsolver/kernels/hegvd_op.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_pw/module_pwdft/onsite_projector.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_pw/module_pwdft/vsep_pw.h"
#include "source_pw/module_pwdft/hamilt_pw.h"

#include "source_pw/module_pwdft/forces.h"
#include "source_pw/module_pwdft/stress_pw.h"
#include "source_hamilt/module_xc/xc_functional.h" // use XC_Functional

#ifdef __DSP
#include "source_base/kernels/dsp/dsp_connector.h"
#endif

#include "source_pw/module_pwdft/setup_pot.h" // mohan add 20250929
#include "source_estate/setup_estate_pw.h" // mohan add 20251005
#include "source_io/module_ctrl/ctrl_output_pw.h" // mohan add 20250927
#include "source_estate/module_charge/chgmixing.h" // use charge mixing, mohan add 20251006 
#include "source_estate/update_pot.h" // mohan add 20251016
#include "source_pw/module_pwdft/update_cell_pw.h" // mohan add 20250309
#include "source_pw/module_pwdft/dftu_pw.h" // mohan add 20250309
#include "source_pw/module_pwdft/deltaspin_pw.h" // mohan add 20250309

#include "source_hamilt/module_xc/exx_info.h" // use GlobalC::exx_info

namespace ModuleESolver
{

template <typename T, typename Device>
ESolver_KS_PW<T, Device>::ESolver_KS_PW()
{
    this->classname = "ESolver_KS_PW";
    this->basisname = "PW";
}

template <typename T, typename Device>
ESolver_KS_PW<T, Device>::~ESolver_KS_PW()
{
    //****************************************************
    // do not add any codes in this deconstructor funcion
    //****************************************************
    // delete Hamilt
    this->deallocate_hamilt();

    // mohan add 2025-10-12
    this->stp.clean();
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::allocate_hamilt(const UnitCell& ucell)
{
	this->p_hamilt = new hamilt::HamiltPW<T, Device>(
			this->pelec->pot, 
			this->pw_wfc, 
			&this->kv, 
			&this->ppcell, 
			&this->dftu,
			&ucell);
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::deallocate_hamilt()
{
    if (this->p_hamilt != nullptr)
    {
        delete reinterpret_cast<hamilt::HamiltPW<T, Device>*>(this->p_hamilt);
        this->p_hamilt = nullptr;
    }
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    ESolver_KS::before_all_runners(ucell, inp);

    //! setup and allocation for pelec, potentials, etc. 
    elecstate::setup_estate_pw<T, Device>(ucell, this->kv, this->sf, this->pelec, this->chr,
      this->locpp, this->ppcell, this->vsep_cell, this->pw_wfc, this->pw_rho,
      this->pw_rhod, this->pw_big, this->solvent, inp);

    this->stp.before_runner(ucell, this->kv, this->sf, *this->pw_wfc, this->ppcell, PARAM.inp);

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT BASIS");

    //! Initialize exx pw
    this->exx_helper.init(ucell, inp, this->pelec->wg);
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::before_scf(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_KS_PW", "before_scf");
    ModuleBase::timer::tick("ESolver_KS_PW", "before_scf");

    ESolver_KS::before_scf(ucell, istep);

    //! Init variables (once the cell has changed)
    pw::update_cell_pw(ucell, this->ppcell, this->kv, this->pw_wfc, PARAM.inp);

    if (ucell.cell_parameter_updated)
    {
        auto* p_psi_init = static_cast<psi::PSIPrepare<T, Device>*>(this->stp.p_psi_init);
        p_psi_init->prepare_init(PARAM.inp.pw_seed);
    }

    //! Init Hamiltonian (cell changed)
    //! Operators in HamiltPW should be reallocated once cell changed
    //! delete Hamilt if not first scf
    this->deallocate_hamilt();

    //! Allocate HamiltPW
    this->allocate_hamilt(ucell);

    //! Setup potentials (local, non-local, sc, +U, DFT-1/2)
    // note: init DFT+U is done here for pw basis for every scf iteration, however, 
    // init DFT+U is done in "before_all_runners" in LCAO basis. This should be refactored, mohan note 2025-11-06
    pw::setup_pot(istep, ucell, this->kv, this->sf, this->pelec, this->Pgrid,
              this->chr, this->locpp, this->ppcell, this->dftu, this->vsep_cell,
              this->stp.get_psi_t(), static_cast<hamilt::Hamilt<T, Device>*>(this->p_hamilt), this->pw_wfc, this->pw_rhod, PARAM.inp);

    // setup psi (electronic wave functions)
    this->stp.init(this->p_hamilt);

    //! Setup EXX helper for Hamiltonian and psi
    exx_helper.before_scf(this->p_hamilt, this->stp.get_psi_t(), PARAM.inp);

    ModuleBase::timer::tick("ESolver_KS_PW", "before_scf");
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::iter_init(UnitCell& ucell, const int istep, const int iter)
{
    ESolver_KS::iter_init(ucell, istep, iter);

    module_charge::chgmixing_ks_pw(iter, this->p_chgmix, this->dftu, PARAM.inp);

    // mohan move harris functional here, 2012-06-05
    // use 'rho(in)' and 'v_h and v_xc'(in)
    this->pelec->f_en.deband_harris = this->pelec->cal_delta_eband(ucell);

    // update local occupations for DFT+U
    // should before lambda loop in DeltaSpin
    pw::iter_init_dftu_pw(iter, istep, this->dftu, this->stp.get_psi_t(), this->pelec->wg, ucell, PARAM.inp);
}

// Temporary, it should be replaced by hsolver later.
template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr)
{
    ModuleBase::timer::tick("ESolver_KS_PW", "hamilt2rho_single");

    // reset energy
    this->pelec->f_en.eband = 0.0;
    this->pelec->f_en.demet = 0.0;

    // setup diagonalization parameters
    hsolver::setup_diago_params_pw<T, Device>(istep, iter, ethr, PARAM.inp);

    bool skip_charge = PARAM.inp.calculation == "nscf" ? true : false;

    // run the inner lambda loop to contrain atomic moments with the DeltaSpin method
    bool skip_solve = pw::run_deltaspin_lambda_loop(iter - 1, this->drho, PARAM.inp);

    if (!skip_solve)
    {
        hsolver::HSolverPW<T, Device> hsolver_pw_obj(this->pw_wfc,
                                                     PARAM.inp.calculation,
                                                     PARAM.inp.basis_type,
                                                     PARAM.inp.ks_solver,
                                                     false,
                                                     PARAM.globalv.use_uspp,
                                                     PARAM.inp.nspin,
                                                     hsolver::DiagoIterAssist<T, Device>::SCF_ITER,
                                                     hsolver::DiagoIterAssist<T, Device>::PW_DIAG_NMAX,
                                                     hsolver::DiagoIterAssist<T, Device>::PW_DIAG_THR,
                                                     hsolver::DiagoIterAssist<T, Device>::need_subspace,
                                                     PARAM.inp.use_k_continuity);

        hsolver_pw_obj.solve(static_cast<hamilt::Hamilt<T, Device>*>(this->p_hamilt), *this->stp.get_psi_t(), this->pelec, this->pelec->ekb.c,
          GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL, skip_charge, ucell.tpiba, ucell.nat);
    }

    // symmetrize the charge density
    Symmetry_rho::symmetrize_rho(PARAM.inp.nspin, this->chr, this->pw_rhod, ucell.symm);

    ModuleBase::timer::tick("ESolver_KS_PW", "hamilt2rho_single");
}


template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver)
{
    // Related to EXX
    if (GlobalC::exx_info.info_global.cal_exx && !exx_helper.op_exx->first_iter)
    {
        this->pelec->set_exx(exx_helper.cal_exx_energy(this->stp.get_psi_t()));
    }

    // deband is calculated from "output" charge density
    this->pelec->f_en.deband = this->pelec->cal_delta_eband(ucell);

    // Call iter_finish() of ESolver_KS
    ESolver_KS::iter_finish(ucell, istep, iter, conv_esolver);

    // D in USPP needs vloc, thus needs update when veff updated
    // calculate the effective coefficient matrix for non-local
    // pp projectors, liuyu 2023-10-24
    if (PARAM.globalv.use_uspp)
    {
        ModuleBase::matrix veff = this->pelec->pot->get_eff_v();
        this->ppcell.cal_effective_D(veff, this->pw_rhod, ucell);
    }

    // Handle EXX-related operations after SCF iteration
    exx_helper.iter_finish(this->pelec, &this->chr, this->stp.get_psi_t(), ucell, PARAM.inp, conv_esolver, iter);

    // check if oscillate for delta_spin method
    pw::check_deltaspin_oscillation(iter, this->drho, this->p_chgmix, PARAM.inp);

    // the output quantities
    ModuleIO::ctrl_iter_pw(istep, iter, conv_esolver, this->stp.psi_cpu, 
              this->kv, this->pw_wfc, PARAM.inp);
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS_PW", "after_scf");
    ModuleBase::timer::tick("ESolver_KS_PW", "after_scf");

    // Calculate kinetic energy density tau for ELF if needed
    if (PARAM.inp.out_elf[0] > 0)
    {
        this->pelec->cal_tau(*(this->stp.psi_cpu));
    }

    ESolver_KS::after_scf(ucell, istep, conv_esolver);

    // Output quantities
    ModuleIO::ctrl_scf_pw<T, Device>(istep, ucell, this->pelec, this->chr, this->kv, this->pw_wfc,
              this->pw_rho, this->pw_rhod, this->pw_big, this->stp,
              this->Pgrid, PARAM.inp);

    ModuleBase::timer::tick("ESolver_KS_PW", "after_scf");
}

template <typename T, typename Device>
double ESolver_KS_PW<T, Device>::cal_energy()
{
    return this->pelec->f_en.etot;
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::cal_force(UnitCell& ucell, ModuleBase::matrix& force)
{
    Forces<double, Device> ff(ucell.nat);

    // mohan add 2025-10-12
    this->stp.update_psi_d();

    // Calculate forces
    ff.cal_force(ucell, force, *this->pelec, this->pw_rhod, &ucell.symm,
                 &this->sf, this->solvent, &this->dftu, &this->locpp, &this->ppcell, 
                 &this->kv, this->pw_wfc, this->stp.get_psi_d());
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::cal_stress(UnitCell& ucell, ModuleBase::matrix& stress)
{
    Stress_PW<double, Device> ss(this->pelec);

    // mohan add 2025-10-12
    this->stp.update_psi_d();

    ss.cal_stress(stress, ucell, this->dftu, this->locpp, this->ppcell, this->pw_rhod,
                  &ucell.symm, &this->sf, &this->kv, this->pw_wfc, this->stp.get_psi_d());

    // external stress
    double unit_transform = 0.0;
    unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    double external_stress[3] = {PARAM.inp.press1, PARAM.inp.press2, PARAM.inp.press3};
    for (int i = 0; i < 3; i++)
    {
        stress(i, i) -= external_stress[i] / unit_transform;
    }
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::after_all_runners(UnitCell& ucell)
{
    ESolver_KS::after_all_runners(ucell);

    ModuleIO::ctrl_runner_pw<T, Device>(ucell, this->pelec, this->pw_wfc, 
            this->pw_rho, this->pw_rhod, this->chr, this->kv, this->stp, 
            this->sf, this->ppcell, this->solvent, this->Pgrid, PARAM.inp); 

    elecstate::teardown_estate_pw<T, Device>(this->pelec, this->vsep_cell);
    
}

template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace ModuleESolver
