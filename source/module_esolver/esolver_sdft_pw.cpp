#include "esolver_sdft_pw.h"

#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_elecstate/elecstate_pw_sdft.h"
#include "module_hamilt_pw/hamilt_stodft/sto_dos.h"
#include "module_hamilt_pw/hamilt_stodft/sto_elecond.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/hsolver_pw_sdft.h"
#include "module_io/output_log.h"
#include "module_io/rho_io.h"
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
{
    classname = "ESolver_SDFT_PW";
    basisname = "PW";
}

ESolver_SDFT_PW::~ESolver_SDFT_PW()
{
}

void ESolver_SDFT_PW::before_all_runners(Input& inp, UnitCell& ucell)
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
    this->pelec->charge->allocate(GlobalV::NSPIN);
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
    this->p_wf_init = new psi::WFInit<std::complex<double>>(GlobalV::init_wfc,
                                                            GlobalV::KS_SOLVER,
                                                            GlobalV::BASIS_TYPE,
                                                            GlobalV::psi_initializer,
                                                            &this->wf,
                                                            this->pw_wfc);
    // 7) set occupatio, redundant?
    if (GlobalV::ocp)
    {
        this->pelec->fixed_weights(GlobalV::ocp_kb, GlobalV::NBANDS, GlobalV::nelec);
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

    size_t size = stowf.chi0->size();

    this->stowf.shchi = new psi::Psi<std::complex<double>>(kv.get_nks(), stowf.nchip_max, wf.npwx, kv.ngk.data());

    ModuleBase::Memory::record("SDFT::shchi", size * sizeof(std::complex<double>));

    if (GlobalV::NBANDS > 0)
    {
        this->stowf.chiortho
            = new psi::Psi<std::complex<double>>(kv.get_nks(), stowf.nchip_max, wf.npwx, kv.ngk.data());
        ModuleBase::Memory::record("SDFT::chiortho", size * sizeof(std::complex<double>));
    }

    // 9) initialize the hsolver
    this->phsol = new hsolver::HSolverPW_SDFT(&kv, pw_wfc, &wf, this->stowf, inp.method_sto);

    return;
}

void ESolver_SDFT_PW::before_scf(const int istep)
{
    ESolver_KS_PW::before_scf(istep);
    if (istep > 0 && INPUT.nbands_sto != 0 && INPUT.initsto_freq > 0 && istep % INPUT.initsto_freq == 0)
    {
        Update_Sto_Orbitals(this->stowf, INPUT.seed_sto);
    }
}

void ESolver_SDFT_PW::iter_finish(int iter)
{
    // this->pelec->print_eigenvalue(GlobalV::ofs_running);
    this->pelec->cal_energies(2);
}

void ESolver_SDFT_PW::after_scf(const int istep)
{
    // save charge difference into files for charge extrapolation
    if (GlobalV::CALCULATION != "scf")
    {
        this->CE.save_files(istep,
                            GlobalC::ucell,
#ifdef __MPI
                            this->pw_big,
#endif
                            this->pelec->charge,
                            &this->sf);
    }

    if (GlobalV::out_chg > 0)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            std::stringstream ssc;
            ssc << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG.cube";
            const double ef_tmp = this->pelec->eferm.get_efval(is);
            ModuleIO::write_rho(
#ifdef __MPI
                pw_big->bz,
                pw_big->nbz,
                pw_rho->nplane,
                pw_rho->startz_current,
#endif
                pelec->charge->rho_save[is],
                is,
                GlobalV::NSPIN,
                0,
                ssc.str(),
                pw_rho->nx,
                pw_rho->ny,
                pw_rho->nz,
                ef_tmp,
                &(GlobalC::ucell));
        }
    }

    ModuleIO::output_convergence_after_scf(this->conv_elec, this->pelec->f_en.etot);
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

    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = GlobalV::PW_DIAG_NMAX;

    this->phsol->solve(this->p_hamilt, this->psi[0], this->pelec, pw_wfc, this->stowf, istep, iter, GlobalV::KS_SOLVER);

    if (GlobalV::MY_STOGROUP == 0)
    {
        Symmetry_rho srho;
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            srho.begin(is, *(this->pelec->charge), pw_rho, GlobalC::Pgrid, GlobalC::ucell.symm);
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
                  pelec->charge);
}

void ESolver_SDFT_PW::after_all_runners()
{

    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << this->pelec->f_en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;
    ModuleIO::write_istate_info(this->pelec->ekb, this->pelec->wg, kv, &(GlobalC::Pkpoints));

    ((hsolver::HSolverPW_SDFT*)phsol)->stoiter.cleanchiallorder(); // release lots of memories

    if (INPUT.out_dos)
    {
        Sto_DOS sto_dos(this->pw_wfc,
                        &this->kv,
                        this->pelec,
                        this->psi,
                        this->p_hamilt,
                        (hsolver::HSolverPW_SDFT*)phsol,
                        &stowf);
        sto_dos.decide_param(INPUT.dos_nche,
                             INPUT.emin_sto,
                             INPUT.emax_sto,
                             INPUT.dos_setemin,
                             INPUT.dos_setemax,
                             INPUT.dos_emin_ev,
                             INPUT.dos_emax_ev,
                             INPUT.dos_scale);
        sto_dos.caldos(INPUT.dos_sigma, INPUT.dos_edelta_ev, INPUT.npart_sto);
    }

    // sKG cost memory, and it should be placed at the end of the program
    if (INPUT.cal_cond)
    {
        Sto_EleCond sto_elecond(&GlobalC::ucell,
                                &this->kv,
                                this->pelec,
                                this->pw_wfc,
                                this->psi,
                                &GlobalC::ppcell,
                                this->p_hamilt,
                                (hsolver::HSolverPW_SDFT*)phsol,
                                &stowf);
        sto_elecond
            .decide_nche(INPUT.cond_dt, INPUT.cond_dtbatch, 1e-8, this->nche_sto, INPUT.emin_sto, INPUT.emax_sto);
        sto_elecond.sKG(INPUT.cond_smear,
                        INPUT.cond_fwhm,
                        INPUT.cond_wcut,
                        INPUT.cond_dw,
                        INPUT.cond_dt,
                        INPUT.cond_nonlocal,
                        INPUT.cond_dtbatch,
                        INPUT.npart_sto);
    }
}

void ESolver_SDFT_PW::others(const int istep)
{
    ModuleBase::TITLE("ESolver_SDFT_PW", "others");

    if (GlobalV::CALCULATION == "nscf")
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

    const double diag_thr = std::max(std::min(1e-5, 0.1 * GlobalV::SCF_THR / std::max(1.0, GlobalV::nelec)), 1e-12);

    std::cout << " DIGA_THR          : " << diag_thr << std::endl;

    this->before_scf(istep);

    this->hamilt2density(istep, iter, diag_thr);

    this->pelec->cal_energies(2);

    ModuleBase::timer::tick("ESolver_SDFT_PW", "nscf");
    return;
}
} // namespace ModuleESolver
