#include "esolver_sdft_pw.h"

#include "module_parameter/parameter.h"
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

ESolver_SDFT_PW::ESolver_SDFT_PW()
    : stoche(PARAM.inp.nche_sto, PARAM.inp.method_sto, PARAM.inp.emax_sto, PARAM.inp.emin_sto)
{
    classname = "ESolver_SDFT_PW";
    basisname = "PW";
}

ESolver_SDFT_PW::~ESolver_SDFT_PW()
{
}

void ESolver_SDFT_PW::before_all_runners(const Input_para& inp, UnitCell& ucell)
{
    // 1) initialize parameters from int Input class
    this->nche_sto = inp.nche_sto;
    this->method_sto = inp.method_sto;

    // 2) run "before_all_runners" in ESolver_KS
    ESolver_KS::before_all_runners(inp, ucell);

    // 3) initialize the pointer for electronic states of SDFT
    this->pelec = new elecstate::ElecStatePW_SDFT(pw_wfc,
                                                  &(chr),
                                                  (K_Vectors*)(&(kv)),
                                                  &ucell,
                                                  &(GlobalC::ppcell),
                                                  this->pw_rhod,
                                                  this->pw_rho,
                                                  pw_big);

    // 4) inititlize the charge density.
    this->pelec->charge->allocate(PARAM.inp.nspin);
    this->pelec->omega = ucell.omega;

    // 5) initialize the potential.
    if (this->pelec->pot == nullptr)
    {
        this->pelec->pot = new elecstate::Potential(pw_rhod,
                                                    pw_rho,
                                                    &ucell,
                                                    &(GlobalC::ppcell.vloc),
                                                    &(sf),
                                                    &(this->pelec->f_en.etxc),
                                                    &(this->pelec->f_en.vtxc));
        GlobalTemp::veff = &(this->pelec->pot->get_effective_v());
    }

    // 6) prepare some parameters for electronic wave functions initilization
    this->p_wf_init = new psi::WFInit<std::complex<double>>(PARAM.inp.init_wfc,
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
    stowf.init(&kv, pw_wfc->npwk_max);
    if (inp.nbands_sto != 0)
    {
        if (inp.initsto_ecut < inp.ecutwfc)
        {
            Init_Sto_Orbitals(this->stowf, inp.seed_sto);
        }
        else
        {
            Init_Sto_Orbitals_Ecut(this->stowf, inp.seed_sto, kv, *pw_wfc, inp.initsto_ecut);
        }
    }
    else
    {
        Init_Com_Orbitals(this->stowf);
    }
    if (this->method_sto == 2)
    {
        stowf.allocate_chiallorder(this->nche_sto);
    }

    size_t size = stowf.chi0->size();

    this->stowf.shchi = new psi::Psi<std::complex<double>>(kv.get_nks(), stowf.nchip_max, wf.npwx, kv.ngk.data());

    ModuleBase::Memory::record("SDFT::shchi", size * sizeof(std::complex<double>));

    if (PARAM.inp.nbands > 0)
    {
        this->stowf.chiortho
            = new psi::Psi<std::complex<double>>(kv.get_nks(), stowf.nchip_max, wf.npwx, kv.ngk.data());
        ModuleBase::Memory::record("SDFT::chiortho", size * sizeof(std::complex<double>));
    }

    return;
}

void ESolver_SDFT_PW::before_scf(const int istep)
{
    ESolver_KS_PW::before_scf(istep);
    if (istep > 0 && PARAM.inp.nbands_sto != 0 && PARAM.inp.initsto_freq > 0 && istep % PARAM.inp.initsto_freq == 0)
    {
        Update_Sto_Orbitals(this->stowf, PARAM.inp.seed_sto);
    }
}

void ESolver_SDFT_PW::iter_finish(int& iter)
{
    // call iter_finish() of ESolver_KS
    ESolver_KS<std::complex<double>>::iter_finish(iter);
}

void ESolver_SDFT_PW::after_scf(const int istep)
{
    // 1) call after_scf() of ESolver_KS_PW
    ESolver_KS_PW<std::complex<double>>::after_scf(istep);
}

void ESolver_SDFT_PW::hamilt2density(int istep, int iter, double ethr)
{
    // reset energy
    this->pelec->f_en.eband = 0.0;
    this->pelec->f_en.demet = 0.0;
    // choose if psi should be diag in subspace
    // be careful that istep start from 0 and iter start from 1
    if (istep == 0 && iter == 1)
    {
        hsolver::DiagoIterAssist<std::complex<double>>::need_subspace = false;
    }
    else
    {
        hsolver::DiagoIterAssist<std::complex<double>>::need_subspace = true;
    }

    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = ethr;

    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = PARAM.inp.pw_diag_nmax;

    // hsolver only exists in this function
    hsolver::HSolverPW_SDFT hsolver_pw_sdft_obj(
                    &this->kv, 
                    this->pw_wfc, 
                    &this->wf, 
                    this->stowf, 
                    this->stoche, 
                    PARAM.inp.calculation,
                    PARAM.inp.basis_type,
                    PARAM.inp.ks_solver,
                    PARAM.inp.use_paw,
                    PARAM.globalv.use_uspp,
                    PARAM.inp.nspin,
                    hsolver::DiagoIterAssist<std::complex<double>>::SCF_ITER,
                    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX,
                    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR,
                    hsolver::DiagoIterAssist<std::complex<double>>::need_subspace,
                    this->init_psi);

    hsolver_pw_sdft_obj.solve(this->p_hamilt,
                              this->psi[0],
                              this->pelec,
                              pw_wfc,
                              this->stowf,
                              istep,
                              iter,
                              false);
    this->init_psi = true;

    // set_diagethr need it
    this->esolver_KS_ne = hsolver_pw_sdft_obj.stoiter.KS_ne;

    if (GlobalV::MY_STOGROUP == 0)
    {
        Symmetry_rho srho;
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            srho.begin(is, *(this->pelec->charge), pw_rho, GlobalC::ucell.symm);
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
}

double ESolver_SDFT_PW::cal_energy()
{
    return this->pelec->f_en.etot;
}

void ESolver_SDFT_PW::cal_force(ModuleBase::matrix& force)
{
    Sto_Forces ff(GlobalC::ucell.nat);

    ff.cal_stoforce(force, *this->pelec, pw_rho, &GlobalC::ucell.symm, &sf, &kv, pw_wfc, this->psi, this->stowf);
}

void ESolver_SDFT_PW::cal_stress(ModuleBase::matrix& stress)
{
    Sto_Stress_PW ss;
    ss.cal_stress(stress,
                  *this->pelec,
                  pw_rho,
                  &GlobalC::ucell.symm,
                  &sf,
                  &kv,
                  pw_wfc,
                  this->psi,
                  this->stowf,
                  pelec->charge,
                  &GlobalC::ppcell,
                  GlobalC::ucell);
}

void ESolver_SDFT_PW::after_all_runners()
{

    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << this->pelec->f_en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;
    ModuleIO::write_istate_info(this->pelec->ekb, this->pelec->wg, kv, &(GlobalC::Pkpoints));

    if (this->method_sto == 2)
    {
        stowf.clean_chiallorder(); // release lots of memories
    }
    if (PARAM.inp.out_dos)
    {
        Sto_DOS sto_dos(this->pw_wfc,
                        &this->kv,
                        this->pelec,
                        this->psi,
                        this->p_hamilt,
                        this->stoche,
                        &stowf);
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

void ESolver_SDFT_PW::others(const int istep)
{
    ModuleBase::TITLE("ESolver_SDFT_PW", "others");

    if (PARAM.inp.calculation == "nscf")
    {
        this->nscf();
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_SDFT_PW::others", "CALCULATION type not supported");
    }

    return;
}

void ESolver_SDFT_PW::nscf()
{
    ModuleBase::TITLE("ESolver_SDFT_PW", "nscf");
    ModuleBase::timer::tick("ESolver_SDFT_PW", "nscf");

    const int istep = 0;

    const int iter = 1;

    const double diag_thr = std::max(std::min(1e-5, 0.1 * PARAM.inp.scf_thr / std::max(1.0, PARAM.inp.nelec)), 1e-12);

    std::cout << " DIGA_THR          : " << diag_thr << std::endl;

    this->before_scf(istep);

    this->hamilt2density(istep, iter, diag_thr);

    this->pelec->cal_energies(2);

    ModuleBase::timer::tick("ESolver_SDFT_PW", "nscf");
    return;
}
} // namespace ModuleESolver
