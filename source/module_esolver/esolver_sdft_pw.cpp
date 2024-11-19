#include "esolver_sdft_pw.h"

#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_elecstate/elecstate_pw_sdft.h"
#include "module_hamilt_pw/hamilt_stodft/sto_dos.h"
#include "module_hamilt_pw/hamilt_stodft/sto_elecond.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/hsolver_pw_sdft.h"
#include "module_io/cube_io.h"
#include "module_io/output_log.h"
#include "module_io/write_istate_info.h"
#include "module_parameter/parameter.h"

#include <algorithm>
#include <fstream>

//-------------------Temporary------------------
#include "module_base/global_variable.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
//----------------------------------------------
//-----force-------------------
#include "module_hamilt_pw/hamilt_stodft/sto_forces.h"
//-----stress------------------
#include "module_hamilt_pw/hamilt_stodft/sto_stress_pw.h"
//---------------------------------------------------

namespace ModuleESolver
{

template <typename T, typename Device>
ESolver_SDFT_PW<T, Device>::ESolver_SDFT_PW()
    : stoche(PARAM.inp.nche_sto, PARAM.inp.method_sto, PARAM.inp.emax_sto, PARAM.inp.emin_sto)
{
    this->classname = "ESolver_SDFT_PW";
    this->basisname = "PW";
}

template <typename T, typename Device>
ESolver_SDFT_PW<T, Device>::~ESolver_SDFT_PW()
{
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::before_all_runners(const Input_para& inp, UnitCell& ucell)
{
    // 1) initialize parameters from int Input class
    this->nche_sto = inp.nche_sto;
    this->method_sto = inp.method_sto;

    // 2) run "before_all_runners" in ESolver_KS
    ESolver_KS<T, Device>::before_all_runners(inp, ucell);

    // 3) initialize the pointer for electronic states of SDFT
    this->pelec = new elecstate::ElecStatePW_SDFT<T, Device>(this->pw_wfc,
                                                             &(this->chr),
                                                             &(this->kv),
                                                             &ucell,
                                                             &(GlobalC::ppcell),
                                                             this->pw_rhod,
                                                             this->pw_rho,
                                                             this->pw_big);

    // 4) inititlize the charge density.
    this->pelec->charge->allocate(PARAM.inp.nspin);
    this->pelec->omega = ucell.omega;

    // 5) initialize the potential.
    if (this->pelec->pot == nullptr)
    {
        this->pelec->pot = new elecstate::Potential(this->pw_rhod,
                                                    this->pw_rho,
                                                    &ucell,
                                                    &(GlobalC::ppcell.vloc),
                                                    &(this->sf),
                                                    &(this->pelec->f_en.etxc),
                                                    &(this->pelec->f_en.vtxc));
    }

    // 6) prepare some parameters for electronic wave functions initilization
    this->p_wf_init = new psi::WFInit<T, Device>(PARAM.inp.init_wfc,
                                                 PARAM.inp.ks_solver,
                                                 PARAM.inp.basis_type,
                                                 PARAM.inp.psi_initializer,
                                                 &this->wf,
                                                 this->pw_wfc);
    // 7) set occupatio, redundant?
    if (PARAM.inp.ocp)
    {
        this->pelec->fixed_weights(PARAM.inp.ocp_kb, PARAM.inp.nbands, PARAM.inp.nelec);
    }

    // 8) initialize the global classes
    this->Init_GlobalC(inp, ucell, GlobalC::ppcell); // temporary

    // 9) initialize the stochastic wave functions
    this->stowf.init(&this->kv, this->pw_wfc->npwk_max);
    if (inp.nbands_sto != 0)
    {
        if (inp.initsto_ecut < inp.ecutwfc)
        {
            this->stowf.init_sto_orbitals(inp.seed_sto);
        }
        else
        {
            this->stowf.init_sto_orbitals_Ecut(inp.seed_sto, this->kv, *this->pw_wfc, inp.initsto_ecut);
        }
    }
    else
    {
        this->stowf.init_com_orbitals();
    }
    if (this->method_sto == 2)
    {
        this->stowf.allocate_chiallorder(this->nche_sto);
    }
    this->stowf.sync_chi0();

    // 10) allocate spaces for \sqrt(f(H))|chi> and |\tilde{chi}>
    size_t size = stowf.chi0->size();
    this->stowf.shchi
        = new psi::Psi<T, Device>(this->kv.get_nks(), this->stowf.nchip_max, this->wf.npwx, this->kv.ngk.data());
    ModuleBase::Memory::record("SDFT::shchi", size * sizeof(T));

    if (PARAM.inp.nbands > 0)
    {
        this->stowf.chiortho
            = new psi::Psi<T, Device>(this->kv.get_nks(), this->stowf.nchip_max, this->wf.npwx, this->kv.ngk.data());
        ModuleBase::Memory::record("SDFT::chiortho", size * sizeof(T));
    }

    return;
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::before_scf(const int istep)
{
    ESolver_KS_PW<T, Device>::before_scf(istep);
    delete reinterpret_cast<hamilt::HamiltPW<double>*>(this->p_hamilt);
    this->p_hamilt = new hamilt::HamiltSdftPW<T, Device>(this->pelec->pot,
                                                         this->pw_wfc,
                                                         &this->kv,
                                                         PARAM.globalv.npol,
                                                         &this->stoche.emin_sto,
                                                         &this->stoche.emax_sto);
    this->p_hamilt_sto = static_cast<hamilt::HamiltSdftPW<T, Device>*>(this->p_hamilt);

    if (istep > 0 && PARAM.inp.nbands_sto != 0 && PARAM.inp.initsto_freq > 0 && istep % PARAM.inp.initsto_freq == 0)
    {
        this->stowf.update_sto_orbitals(PARAM.inp.seed_sto);
    }
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::iter_finish(const int istep, int& iter)
{
    // call iter_finish() of ESolver_KS
    ESolver_KS<T, Device>::iter_finish(istep, iter);
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::after_scf(const int istep)
{
    // 1) call after_scf() of ESolver_KS_PW
    ESolver_KS_PW<T, Device>::after_scf(istep);
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::hamilt2density_single(int istep, int iter, double ethr)
{
    ModuleBase::TITLE("ESolver_SDFT_PW", "hamilt2density");
    ModuleBase::timer::tick("ESolver_SDFT_PW", "hamilt2density");

    // reset energy
    this->pelec->f_en.eband = 0.0;
    this->pelec->f_en.demet = 0.0;
    // choose if psi should be diag in subspace
    // be careful that istep start from 0 and iter start from 1
    if (istep == 0 && iter == 1 || PARAM.inp.calculation == "nscf")
    {
        hsolver::DiagoIterAssist<T, Device>::need_subspace = false;
    }
    else
    {
        hsolver::DiagoIterAssist<T, Device>::need_subspace = true;
    }

    bool skip_charge = PARAM.inp.calculation == "nscf" ? true : false;
    hsolver::DiagoIterAssist<T, Device>::PW_DIAG_THR = ethr;
    hsolver::DiagoIterAssist<T, Device>::PW_DIAG_NMAX = PARAM.inp.pw_diag_nmax;

    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------for psi init guess!!!!--------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    if (!PARAM.inp.psi_initializer && PARAM.inp.basis_type == "pw" && this->init_psi == false)
    {
        for (int ik = 0; ik < this->pw_wfc->nks; ++ik)
        {
            //! Update Hamiltonian from other kpoint to the given one
            this->p_hamilt->updateHk(ik);

            if (this->kspw_psi->get_nbands() > 0 && GlobalV::MY_STOGROUP == 0)
            {
                //! Fix the wavefunction to initialize at given kpoint
                this->kspw_psi->fix_k(ik);

                /// for psi init guess!!!!
                hamilt::diago_PAO_in_pw_k2(this->ctx, ik, *(this->kspw_psi), this->pw_wfc, &this->wf, this->p_hamilt);
            }

        }
    }
    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------END: for psi init guess!!!!--------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------


    // hsolver only exists in this function
    hsolver::HSolverPW_SDFT<T, Device> hsolver_pw_sdft_obj(&this->kv,
                                                           this->pw_wfc,
                                                           this->stowf,
                                                           this->stoche,
                                                           this->p_hamilt_sto,
                                                           PARAM.inp.calculation,
                                                           PARAM.inp.basis_type,
                                                           PARAM.inp.ks_solver,
                                                           PARAM.inp.use_paw,
                                                           PARAM.globalv.use_uspp,
                                                           PARAM.inp.nspin,
                                                           hsolver::DiagoIterAssist<T, Device>::SCF_ITER,
                                                           hsolver::DiagoIterAssist<T, Device>::PW_DIAG_NMAX,
                                                           hsolver::DiagoIterAssist<T, Device>::PW_DIAG_THR,
                                                           hsolver::DiagoIterAssist<T, Device>::need_subspace);

    hsolver_pw_sdft_obj.solve(this->p_hamilt,
                              this->kspw_psi[0],
                              this->psi[0],
                              this->pelec,
                              this->pw_wfc,
                              this->stowf,
                              istep,
                              iter,
                              skip_charge);
    this->init_psi = true;

    // set_diagethr need it
    this->esolver_KS_ne = hsolver_pw_sdft_obj.stoiter.KS_ne;

    if (GlobalV::MY_STOGROUP == 0)
    {
        Symmetry_rho srho;
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            srho.begin(is, *(this->pelec->charge), this->pw_rho, GlobalC::ucell.symm);
        }
        this->pelec->f_en.deband = this->pelec->cal_delta_eband();
    }
    else
    {
#ifdef __MPI
        if (ModuleSymmetry::Symmetry::symm_flag == 1)
        {
            MPI_Barrier(MPI_COMM_WORLD);
        }
#endif
    }
#ifdef __MPI
    MPI_Bcast(&(this->pelec->f_en.deband), 1, MPI_DOUBLE, 0, PARAPW_WORLD);
#endif
    ModuleBase::timer::tick("ESolver_SDFT_PW", "hamilt2density");
}

template <typename T, typename Device>
double ESolver_SDFT_PW<T, Device>::cal_energy()
{
    return this->pelec->f_en.etot;
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::cal_force(ModuleBase::matrix& force)
{
    Sto_Forces<double, Device> ff(GlobalC::ucell.nat);

    ff.cal_stoforce(force,
                    *this->pelec,
                    this->pw_rho,
                    &GlobalC::ucell.symm,
                    &this->sf,
                    &this->kv,
                    this->pw_wfc,
                    GlobalC::ppcell,
                    GlobalC::ucell,
                    *this->kspw_psi,
                    this->stowf);
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::cal_stress(ModuleBase::matrix& stress)
{
    Sto_Stress_PW<double, Device> ss;
    ss.cal_stress(stress,
                  *this->pelec,
                  this->pw_rho,
                  &GlobalC::ucell.symm,
                  &this->sf,
                  &this->kv,
                  this->pw_wfc,
                  *this->kspw_psi,
                  this->stowf,
                  this->pelec->charge,
                  &GlobalC::ppcell,
                  GlobalC::ucell);
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::after_all_runners()
{
    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << this->pelec->f_en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;
    ModuleIO::write_istate_info(this->pelec->ekb, this->pelec->wg, this->kv, &(GlobalC::Pkpoints));
}

template <>
void ESolver_SDFT_PW<std::complex<double>, base_device::DEVICE_CPU>::after_all_runners()
{

    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << this->pelec->f_en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;
    ModuleIO::write_istate_info(this->pelec->ekb, this->pelec->wg, this->kv, &(GlobalC::Pkpoints));

    if (this->method_sto == 2)
    {
        stowf.clean_chiallorder(); // release lots of memories
    }
    if (PARAM.inp.out_dos)
    {
        Sto_DOS sto_dos(this->pw_wfc, &this->kv, this->pelec, this->psi, this->p_hamilt, this->stoche, &stowf);
        sto_dos.decide_param(PARAM.inp.dos_nche,
                             PARAM.inp.emin_sto,
                             PARAM.inp.emax_sto,
                             PARAM.globalv.dos_setemin,
                             PARAM.globalv.dos_setemax,
                             PARAM.inp.dos_emin_ev,
                             PARAM.inp.dos_emax_ev,
                             PARAM.inp.dos_scale);
        sto_dos.caldos(PARAM.inp.dos_sigma, PARAM.inp.dos_edelta_ev, PARAM.inp.npart_sto);
    }

    // sKG cost memory, and it should be placed at the end of the program
    if (PARAM.inp.cal_cond)
    {
        Sto_EleCond sto_elecond(&GlobalC::ucell,
                                &this->kv,
                                this->pelec,
                                this->pw_wfc,
                                this->psi,
                                &GlobalC::ppcell,
                                this->p_hamilt,
                                this->stoche,
                                &stowf);
        sto_elecond.decide_nche(PARAM.inp.cond_dt, 1e-8, this->nche_sto, PARAM.inp.emin_sto, PARAM.inp.emax_sto);
        sto_elecond.sKG(PARAM.inp.cond_smear,
                        PARAM.inp.cond_fwhm,
                        PARAM.inp.cond_wcut,
                        PARAM.inp.cond_dw,
                        PARAM.inp.cond_dt,
                        PARAM.inp.cond_nonlocal,
                        PARAM.inp.npart_sto);
    }
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::others(const int istep)
{
    ModuleBase::TITLE("ESolver_SDFT_PW", "others");

    ModuleBase::WARNING_QUIT("ESolver_SDFT_PW<T, Device>::others", "CALCULATION type not supported");

    return;
}

// template class ESolver_SDFT_PW<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_SDFT_PW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
// template class ESolver_SDFT_PW<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_SDFT_PW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace ModuleESolver
