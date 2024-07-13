#include "esolver_ks.h"

#include <ctime>
#ifdef __MPI
#include <mpi.h>
#else
#include <chrono>
#endif
#include <iostream>

#include "module_base/timer.h"
#include "module_io/input.h"
#include "module_io/json_output/init_info.h"
#include "module_io/print_info.h"
#include "module_parameter/parameter.h"
//--------------Temporary----------------
#include "module_base/global_variable.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
//---------------------------------------
#ifdef USE_PAW
#include "module_base/parallel_common.h"
#include "module_cell/module_paw/paw_cell.h"
#endif
#include "module_io/json_output/output_info.h"

namespace ModuleESolver
{

//------------------------------------------------------------------------------
//! the 1st function of ESolver_KS: constructor
//! mohan add 2024-05-11
// in future, the initialize of ESolver_KS should not be based on the
// assumption that INPUT has been initialized, mohan 2024-05-12
//------------------------------------------------------------------------------
template <typename T, typename Device>
ESolver_KS<T, Device>::ESolver_KS()
{
    classname = "ESolver_KS";
    basisname = "PLEASE ADD BASISNAME FOR CURRENT ESOLVER.";

    // should not use GlobalV here, mohan 2024-05-12
    scf_thr = GlobalV::SCF_THR;
    drho = 0.0;

    // should not use GlobalV here, mohan 2024-05-12
    maxniter = GlobalV::SCF_NMAX;
    niter = maxniter;

    // should not use GlobalV here, mohan 2024-05-12
    out_freq_elec = GlobalV::OUT_FREQ_ELEC;

    // pw_rho = new ModuleBase::PW_Basis();
    // temporary, it will be removed
    pw_wfc = new ModulePW::PW_Basis_K_Big(GlobalV::device_flag, GlobalV::precision_flag);
    ModulePW::PW_Basis_K_Big* tmp = static_cast<ModulePW::PW_Basis_K_Big*>(pw_wfc);

    // should not use INPUT here, mohan 2024-05-12
    tmp->setbxyz(PARAM.inp.bx, PARAM.inp.by, PARAM.inp.bz);

    ///----------------------------------------------------------
    /// charge mixing
    ///----------------------------------------------------------
    p_chgmix = new Charge_Mixing();
    p_chgmix->set_rhopw(this->pw_rho, this->pw_rhod);

    ///----------------------------------------------------------
    /// wavefunc
    ///----------------------------------------------------------
    this->wf.init_wfc = PARAM.inp.init_wfc;
    this->wf.mem_saver = PARAM.inp.mem_saver;
    this->wf.out_wfc_pw = PARAM.inp.out_wfc_pw;
    this->wf.out_wfc_r = PARAM.inp.out_wfc_r;
}

//------------------------------------------------------------------------------
//! the 2nd function of ESolver_KS: deconstructor
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename T, typename Device>
ESolver_KS<T, Device>::~ESolver_KS()
{
    delete this->psi;
    delete this->pw_wfc;
    delete this->p_hamilt;
    delete this->phsol;
    delete this->p_chgmix;
}

//------------------------------------------------------------------------------
//! the 3rd function of ESolver_KS: before_all_runners
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename T, typename Device>
void ESolver_KS<T, Device>::before_all_runners(const Input_para& inp, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_KS", "before_all_runners");

    //! 1) initialize "before_all_runniers" in ESolver_FP
    ESolver_FP::before_all_runners(inp, ucell);

    //! 2) setup the charge mixing parameters
    p_chgmix->set_mixing(GlobalV::MIXING_MODE,
                         GlobalV::MIXING_BETA,
                         GlobalV::MIXING_NDIM,
                         GlobalV::MIXING_GG0,
                         GlobalV::MIXING_TAU,
                         GlobalV::MIXING_BETA_MAG,
                         GlobalV::MIXING_GG0_MAG,
                         GlobalV::MIXING_GG0_MIN,
                         GlobalV::MIXING_ANGLE,
                         GlobalV::MIXING_DMR);

    /// PAW Section
#ifdef USE_PAW
    if (GlobalV::use_paw)
    {
        int* atom_type = nullptr;
        double** atom_coord = nullptr;
        std::vector<std::string> filename_list;

        atom_type = new int[ucell.nat];
        atom_coord = new double*[ucell.nat];
        filename_list.resize(ucell.ntype);

        for (int ia = 0; ia < ucell.nat; ia++)
        {
            atom_coord[ia] = new double[3];
        }

        int iat = 0;
        for (int it = 0; it < ucell.ntype; it++)
        {
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                atom_type[iat] = it;
                atom_coord[iat][0] = ucell.atoms[it].taud[ia].x;
                atom_coord[iat][1] = ucell.atoms[it].taud[ia].y;
                atom_coord[iat][2] = ucell.atoms[it].taud[ia].z;
                iat++;
            }
        }

        if (GlobalV::MY_RANK == 0)
        {
            std::ifstream ifa(GlobalV::stru_file.c_str(), std::ios::in);
            if (!ifa)
            {
                ModuleBase::WARNING_QUIT("set_libpaw_files", "can not open stru file");
            }

            std::string line;
            while (!ifa.eof())
            {
                getline(ifa, line);
                if (line.find("PAW_FILES") != std::string::npos) {
                    break;
                }
            }

            for (int it = 0; it < ucell.ntype; it++)
            {
                ifa >> filename_list[it];
            }
        }
#ifdef __MPI
        for (int it = 0; it < ucell.ntype; it++)
        {
            Parallel_Common::bcast_string(filename_list[it]);
        }
#endif

        GlobalC::paw_cell.init_paw_cell(inp.ecutwfc,
                                        inp.cell_factor,
                                        ucell.omega,
                                        ucell.nat,
                                        ucell.ntype,
                                        atom_type,
                                        (const double**)atom_coord,
                                        filename_list);

        for (int iat = 0; iat < ucell.nat; iat++)
        {
            delete[] atom_coord[iat];
        }
        delete[] atom_coord;
        delete[] atom_type;
    }
#endif
    /// End PAW

    //! 3) calculate the electron number
    ucell.cal_nelec(GlobalV::nelec);

    //! 4) it has been established that
    // xc_func is same for all elements, therefore
    // only the first one if used
    if (GlobalV::use_paw)
    {
        XC_Functional::set_xc_type(GlobalV::DFT_FUNCTIONAL);
    }
    else
    {
        XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    //! 5) ESolver depends on the Symmetry module
    // symmetry analysis should be performed every time the cell is changed
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        ucell.symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    ucell.print_cell_cif("STRU.cif");

    //! 6) Setup the k points according to symmetry.
    this->kv.set(ucell.symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, ucell.G, ucell.latvec, GlobalV::ofs_running);

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

    //! 7) print information
    Print_Info::setup_parameters(ucell, this->kv);

    //! 8) new plane wave basis, fft grids, etc.
#ifdef __MPI
    this->pw_wfc->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif

    this->pw_wfc->initgrids(inp.ref_cell_factor * ucell.lat0,
                            ucell.latvec,
                            this->pw_rho->nx,
                            this->pw_rho->ny,
                            this->pw_rho->nz);

    this->pw_wfc->initparameters(false, inp.ecutwfc, this->kv.get_nks(), this->kv.kvec_d.data());

    // the MPI allreduce should not be here, mohan 2024-05-12
#ifdef __MPI
    if (inp.pw_seed > 0)
    {
        MPI_Allreduce(MPI_IN_PLACE, &this->pw_wfc->ggecut, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    // qianrui add 2021-8-13 to make different kpar parameters can get the same
    // results
#endif

    this->pw_wfc->ft.fft_mode = inp.fft_mode;

    this->pw_wfc->setuptransform();

    //! 9) initialize the number of plane waves for each k point
    for (int ik = 0; ik < this->kv.get_nks(); ++ik)
    {
        this->kv.ngk[ik] = this->pw_wfc->npwk[ik];
    }

    this->pw_wfc->collect_local_pw(inp.erf_ecut, inp.erf_height, inp.erf_sigma);

    this->print_wfcfft(inp, GlobalV::ofs_running);

    //! 10) initialize the real-space uniform grid for FFT and parallel
    //! distribution of plane waves
    GlobalC::Pgrid.init(this->pw_rhod->nx,
                        this->pw_rhod->ny,
                        this->pw_rhod->nz,
                        this->pw_rhod->nplane,
                        this->pw_rhod->nrxx,
                        pw_big->nbz,
                        pw_big->bz);

    //! 11) calculate the structure factor
    this->sf.setup_structure_factor(&ucell, this->pw_rhod);

    //! 12) initialize the charge extrapolation method if necessary
    CE.Init_CE(ucell.nat);

#ifdef USE_PAW
    if (GlobalV::use_paw)
    {
        GlobalC::paw_cell.set_libpaw_ecut(inp.ecutwfc / 2.0,
                                          inp.ecutwfc / 2.0); // in Hartree
        GlobalC::paw_cell.set_libpaw_fft(this->pw_wfc->nx,
                                         this->pw_wfc->ny,
                                         this->pw_wfc->nz,
                                         this->pw_wfc->nx,
                                         this->pw_wfc->ny,
                                         this->pw_wfc->nz,
                                         this->pw_wfc->startz,
                                         this->pw_wfc->numz);
#ifdef __MPI
        if (GlobalV::RANK_IN_POOL == 0)
        {
            GlobalC::paw_cell.prepare_paw();
        }
#else
        GlobalC::paw_cell.prepare_paw();
#endif
        GlobalC::paw_cell.set_sij();

        GlobalC::paw_cell.set_eigts(this->pw_wfc->nx,
                                    this->pw_wfc->ny,
                                    this->pw_wfc->nz,
                                    this->sf.eigts1.c,
                                    this->sf.eigts2.c,
                                    this->sf.eigts3.c);

        std::vector<std::vector<double>> rhoijp;
        std::vector<std::vector<int>> rhoijselect;
        std::vector<int> nrhoijsel;
#ifdef __MPI
        if (GlobalV::RANK_IN_POOL == 0)
        {
            GlobalC::paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

            for (int iat = 0; iat < ucell.nat; iat++)
            {
                GlobalC::paw_cell.set_rhoij(iat,
                                            nrhoijsel[iat],
                                            rhoijselect[iat].size(),
                                            rhoijselect[iat].data(),
                                            rhoijp[iat].data());
            }
        }
#else
        GlobalC::paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

        for (int iat = 0; iat < ucell.nat; iat++)
        {
            GlobalC::paw_cell.set_rhoij(iat,
                                        nrhoijsel[iat],
                                        rhoijselect[iat].size(),
                                        rhoijselect[iat].data(),
                                        rhoijp[iat].data());
        }
#endif
    }
#endif
}

//------------------------------------------------------------------------------
//! the 4th function of ESolver_KS: init_after_vc
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename T, typename Device>
void ESolver_KS<T, Device>::init_after_vc(const Input_para& inp, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_KS", "init_after_vc");

    ESolver_FP::init_after_vc(inp, ucell);

    if (GlobalV::md_prec_level == 2)
    {
        // initialize the real-space uniform grid for FFT and parallel
        // distribution of plane waves
        GlobalC::Pgrid.init(this->pw_rhod->nx,
                            this->pw_rhod->ny,
                            this->pw_rhod->nz,
                            this->pw_rhod->nplane,
                            this->pw_rhod->nrxx,
                            pw_big->nbz,
                            pw_big->bz); // mohan add 2010-07-22, update 2011-05-04

        // Calculate Structure factor
        this->sf.setup_structure_factor(&ucell, this->pw_rhod);
    }
}

//------------------------------------------------------------------------------
//! the 5th function of ESolver_KS: hamilt2density
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename T, typename Device>
void ESolver_KS<T, Device>::hamilt2density(const int istep, const int iter, const double ethr)
{
    ModuleBase::timer::tick(this->classname, "hamilt2density");
    // Temporarily, before HSolver is constructed, it should be overrided by
    // LCAO, PW, SDFT and TDDFT.
    // After HSolver is constructed, LCAO, PW, SDFT should delete their own
    // hamilt2density() and use:
    // this->phsol->solve(this->phamilt, this->pes, this->wf, ETHR);
    ModuleBase::timer::tick(this->classname, "hamilt2density");
}

//------------------------------------------------------------------------------
//! the 6th function of ESolver_KS: print_wfcfft
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename T, typename Device>
void ESolver_KS<T, Device>::print_wfcfft(const Input_para& inp, std::ofstream& ofs)
{
    ofs << "\n\n\n\n";
    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
           ">>>>"
        << std::endl;
    ofs << " |                                                                 "
           "   |"
        << std::endl;
    ofs << " | Setup plane waves of wave functions:                            "
           "   |"
        << std::endl;
    ofs << " | Use the energy cutoff and the lattice vectors to generate the   "
           "   |"
        << std::endl;
    ofs << " | dimensions of FFT grid. The number of FFT grid on each "
           "processor   |"
        << std::endl;
    ofs << " | is 'nrxx'. The number of plane wave basis in reciprocal space "
           "is   |"
        << std::endl;
    ofs << " | different for charege/potential and wave functions. We also set "
           "   |"
        << std::endl;
    ofs << " | the 'sticks' for the parallel of FFT. The number of plane wave "
           "of  |"
        << std::endl;
    ofs << " | each k-point is 'npwk[ik]' in each processor                    "
           "   |"
        << std::endl;
    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
           "<<<<"
        << std::endl;
    ofs << "\n\n\n\n";
    ofs << "\n SETUP PLANE WAVES FOR WAVE FUNCTIONS" << std::endl;

    double ecut = inp.ecutwfc;
    if (std::abs(ecut - this->pw_wfc->gk_ecut * this->pw_wfc->tpiba2) > 1e-6)
    {
        ecut = this->pw_wfc->gk_ecut * this->pw_wfc->tpiba2;
        ofs << "Energy cutoff for wavefunc is incompatible with nx, ny, nz and "
               "it will be reduced!"
            << std::endl;
    }
    ModuleBase::GlobalFunc::OUT(ofs, "energy cutoff for wavefunc (unit:Ry)", ecut);
    ModuleBase::GlobalFunc::OUT(ofs,
                                "fft grid for wave functions",
                                this->pw_wfc->nx,
                                this->pw_wfc->ny,
                                this->pw_wfc->nz);
    ModuleBase::GlobalFunc::OUT(ofs, "number of plane waves", this->pw_wfc->npwtot);
    ModuleBase::GlobalFunc::OUT(ofs, "number of sticks", this->pw_wfc->nstot);

    ofs << "\n PARALLEL PW FOR WAVE FUNCTIONS" << std::endl;
    ofs << " " << std::setw(8) << "PROC" << std::setw(15) << "COLUMNS(POT)" << std::setw(15) << "PW" << std::endl;

    for (int i = 0; i < GlobalV::NPROC_IN_POOL; ++i)
    {
        ofs << " " << std::setw(8) << i + 1 << std::setw(15) << this->pw_wfc->nst_per[i] << std::setw(15)
            << this->pw_wfc->npw_per[i] << std::endl;
    }

    ofs << " --------------- sum -------------------" << std::endl;
    ofs << " " << std::setw(8) << GlobalV::NPROC_IN_POOL << std::setw(15) << this->pw_wfc->nstot << std::setw(15)
        << this->pw_wfc->npwtot << std::endl;
    ModuleBase::GlobalFunc::DONE(ofs, "INIT PLANEWAVE");
}

//------------------------------------------------------------------------------
//! the 7th function of ESolver_KS: run
//! mohan add 2024-05-11
//! 2) before_scf (electronic iteration loops)
//! 3) run charge density
//! 4) SCF iterations
//! 5) write head
//! 6) initialization of SCF iterations
//! 7) use Hamiltonian to obtain charge density
//! 8) for MPI: STOGROUP? need to rewrite
//! 9) update potential
//! 10) finish scf iterations
//! 11) get mtaGGA related parameters
//! 12) Json, need to be moved to somewhere else
//! 13) check convergence
//! 14) add Json of efermi energy converge
//! 15) after scf
//! 16) Json again
//------------------------------------------------------------------------------
template <typename T, typename Device>
void ESolver_KS<T, Device>::runner(const int istep, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_KS", "runner");

    ModuleBase::timer::tick(this->classname, "runner");

    // 2) before_scf (electronic iteration loops)
    this->before_scf(istep);

    // 3) write charge density
    if (PARAM.inp.dm_to_rho)
    {
        ModuleBase::timer::tick(this->classname, "runner");
        return; // nothing further is needed
    }

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT SCF");

    bool firstscf = true;
    this->conv_elec = false;
    this->niter = this->maxniter;

    // 4) SCF iterations
    std::cout << " * * * * * *\n << Start SCF iteration." << std::endl;
    for (int iter = 1; iter <= this->maxniter; ++iter)
    {
        // 5) write head
        this->write_head(GlobalV::ofs_running, istep, iter);

#ifdef __MPI
        auto iterstart = MPI_Wtime();
#else
        auto iterstart = std::chrono::system_clock::now();
#endif
        double diag_ethr = this->phsol->set_diagethr(istep, iter, drho);

        // 6) initialization of SCF iterations
        this->iter_init(istep, iter);

        // 7) use Hamiltonian to obtain charge density
        this->hamilt2density(istep, iter, diag_ethr);

        // 8) for MPI: STOGROUP? need to rewrite
        //<Temporary> It may be changed when more clever parallel algorithm is
        // put forward.
        // When parallel algorithm for bands are adopted. Density will only be
        // treated in the first group.
        //(Different ranks should have abtained the same, but small differences
        // always exist in practice.)
        // Maybe in the future, density and wavefunctions should use different
        // parallel algorithms, in which they do not occupy all processors, for
        // example wavefunctions uses 20 processors while density uses 10.
        if (GlobalV::MY_STOGROUP == 0)
        {
            // double drho = this->estate.caldr2();
            // EState should be used after it is constructed.

            drho = p_chgmix->get_drho(pelec->charge, GlobalV::nelec);
            double hsolver_error = 0.0;
            if (firstscf)
            {
                firstscf = false;
                hsolver_error = this->phsol->cal_hsolerror();
                // The error of HSolver is larger than drho,
                // so a more precise HSolver should be excuconv_elected.
                if (hsolver_error > drho)
                {
                    diag_ethr = this->phsol->reset_diagethr(GlobalV::ofs_running, hsolver_error, drho);
                    this->hamilt2density(istep, iter, diag_ethr);
                    drho = p_chgmix->get_drho(pelec->charge, GlobalV::nelec);
                    hsolver_error = this->phsol->cal_hsolerror();
                }
            }
            // mixing will restart at this->p_chgmix->mixing_restart steps
            if (drho <= GlobalV::MIXING_RESTART && GlobalV::MIXING_RESTART > 0.0
                && this->p_chgmix->mixing_restart_step > iter)
            {
                this->p_chgmix->mixing_restart_step = iter + 1;
            }

            // drho will be 0 at this->p_chgmix->mixing_restart step, which is
            // not ground state
            bool not_restart_step = !(iter == this->p_chgmix->mixing_restart_step && GlobalV::MIXING_RESTART > 0.0);
            // SCF will continue if U is not converged for uramping calculation
            bool is_U_converged = true;
            // to avoid unnecessary dependence on dft+u, refactor is needed
#ifdef __LCAO
            if (GlobalV::dft_plus_u)
            {
                is_U_converged = GlobalC::dftu.u_converged();
            }
#endif

            this->conv_elec = (drho < this->scf_thr && not_restart_step && is_U_converged);

            // If drho < hsolver_error in the first iter or drho < scf_thr, we
            // do not change rho.
            if (drho < hsolver_error || this->conv_elec)
            {
                if (drho < hsolver_error)
                {
                    GlobalV::ofs_warning << " drho < hsolver_error, keep "
                                            "charge density unchanged."
                                         << std::endl;
                }
            }
            else
            {
                //----------charge mixing---------------
                // mixing will restart after this->p_chgmix->mixing_restart
                // steps
                if (GlobalV::MIXING_RESTART > 0 && iter == this->p_chgmix->mixing_restart_step - 1
                    && drho <= GlobalV::MIXING_RESTART)
                {
                    // do not mix charge density
                }
                else
                {
                    p_chgmix->mix_rho(pelec->charge); // update chr->rho by mixing
                }
                if (GlobalV::SCF_THR_TYPE == 2)
                {
                    pelec->charge->renormalize_rho(); // renormalize rho in R-space would
                                                      // induce a error in K-space
                }
                //----------charge mixing done-----------
            }
        }
#ifdef __MPI
        MPI_Bcast(&drho, 1, MPI_DOUBLE, 0, PARAPW_WORLD);
        MPI_Bcast(&this->conv_elec, 1, MPI_DOUBLE, 0, PARAPW_WORLD);
        MPI_Bcast(pelec->charge->rho[0], this->pw_rhod->nrxx, MPI_DOUBLE, 0, PARAPW_WORLD);
#endif

        // 9) update potential
        // Hamilt should be used after it is constructed.
        // this->phamilt->update(conv_elec);
        this->update_pot(istep, iter);

        // 10) finish scf iterations
        this->iter_finish(iter);
#ifdef __MPI
        double duration = (double)(MPI_Wtime() - iterstart);
#else
        double duration
            = (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - iterstart))
                  .count()
              / static_cast<double>(1e6);
#endif

        // 11) get mtaGGA related parameters
        double dkin = 0.0; // for meta-GGA
        if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
        {
            dkin = p_chgmix->get_dkin(pelec->charge, GlobalV::nelec);
        }
        this->print_iter(iter, drho, dkin, duration, diag_ethr);

        // 12) Json, need to be moved to somewhere else
#ifdef __RAPIDJSON
        // add Json of scf mag
        Json::add_output_scf_mag(GlobalC::ucell.magnet.tot_magnetization,
                                 GlobalC::ucell.magnet.abs_magnetization,
                                 this->pelec->f_en.etot * ModuleBase::Ry_to_eV,
                                 this->pelec->f_en.etot_delta * ModuleBase::Ry_to_eV,
                                 drho,
                                 duration);
#endif //__RAPIDJSON

        // 13) check convergence
        if (this->conv_elec)
        {
            this->niter = iter;
            bool stop = this->do_after_converge(iter);
            if (stop)
            {
                break;
            }
        }

        // notice for restart
        if (GlobalV::MIXING_RESTART > 0 && iter == this->p_chgmix->mixing_restart_step - 1 && iter != GlobalV::SCF_NMAX)
        {
            std::cout << " SCF restart after this step!" << std::endl;
        }
    } // end scf iterations
    std::cout << " >> Leave SCF iteration.\n * * * * * *" << std::endl;
#ifdef __RAPIDJSON
    // 14) add Json of efermi energy converge
    Json::add_output_efermi_converge(this->pelec->eferm.ef * ModuleBase::Ry_to_eV, this->conv_elec);
#endif //__RAPIDJSON
    // 15) after scf
    this->after_scf(istep);
    ModuleBase::timer::tick(this->classname, "runner");

    // 16) Json again
#ifdef __RAPIDJSON
    // add nkstot,nkstot_ibz to output json
    int Jnkstot = this->pelec->klist->get_nkstot();
    Json::add_nkstot(Jnkstot);
#endif //__RAPIDJSON
    return;
};

//------------------------------------------------------------------------------
//! the 8th function of ESolver_KS: print_head
//! mohan add 2024-05-12
//------------------------------------------------------------------------------
template <typename T, typename Device>
void ESolver_KS<T, Device>::print_head()
{
    std::cout << " " << std::setw(7) << "ITER";

    if (GlobalV::NSPIN == 2)
    {
        std::cout << std::setw(10) << "TMAG";
        std::cout << std::setw(10) << "AMAG";
    }

    std::cout << std::setw(15) << "ETOT(eV)";
    std::cout << std::setw(15) << "EDIFF(eV)";
    std::cout << std::setw(11) << "DRHO";

    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
        std::cout << std::setw(11) << "DKIN";
    }

    std::cout << std::setw(11) << "TIME(s)" << std::endl;
}

//------------------------------------------------------------------------------
//! the 8th function of ESolver_KS: print_iter
//! mohan add 2024-05-12
//------------------------------------------------------------------------------
template <typename T, typename Device>
void ESolver_KS<T, Device>::print_iter(const int iter,
                                       const double drho,
                                       const double dkin,
                                       const double duration,
                                       const double ethr)
{
    this->pelec->print_etot(this->conv_elec, iter, drho, dkin, duration, PARAM.inp.printe, ethr);
}

//------------------------------------------------------------------------------
//! the 9th function of ESolver_KS: write_head
//! mohan add 2024-05-12
//------------------------------------------------------------------------------
template <typename T, typename Device>
void ESolver_KS<T, Device>::write_head(std::ofstream& ofs_running, const int istep, const int iter)
{
    ofs_running << "\n " << this->basisname << " ALGORITHM --------------- ION=" << std::setw(4) << istep + 1
                << "  ELEC=" << std::setw(4) << iter << "--------------------------------\n";
}

//------------------------------------------------------------------------------
//! the 10th function of ESolver_KS: getnieter
//! mohan add 2024-05-12
//------------------------------------------------------------------------------
template <typename T, typename Device>
int ESolver_KS<T, Device>::get_niter()
{
    return this->niter;
}

//------------------------------------------------------------------------------
//! the 11th function of ESolver_KS: get_maxniter
//! tqzhao add 2024-05-15
//------------------------------------------------------------------------------
template <typename T, typename Device>
int ESolver_KS<T, Device>::get_maxniter()
{
    return this->maxniter;
}

//------------------------------------------------------------------------------
//! the 12th function of ESolver_KS: get_conv_elec
//! tqzhao add 2024-05-15
//------------------------------------------------------------------------------
template <typename T, typename Device>
bool ESolver_KS<T, Device>::get_conv_elec()
{
    return this->conv_elec;
}

//------------------------------------------------------------------------------
//! the 13th function of ESolver_KS: create_Output_Rho
//! mohan add 2024-05-12
//------------------------------------------------------------------------------
template <typename T, typename Device>
ModuleIO::Output_Rho ESolver_KS<T, Device>::create_Output_Rho(int is, int iter, const std::string& prefix)
{
    const int precision = 3;
    std::string tag = "CHG";
    if (PARAM.inp.dm_to_rho)
    {
        return ModuleIO::Output_Rho(this->pw_big,
                                    this->pw_rhod,
                                    is,
                                    GlobalV::NSPIN,
                                    pelec->charge->rho[is],
                                    iter,
                                    this->pelec->eferm.get_efval(is),
                                    &(GlobalC::ucell),
                                    GlobalV::global_out_dir,
                                    precision,
                                    tag,
                                    prefix);
    }
    return ModuleIO::Output_Rho(this->pw_big,
                                this->pw_rhod,
                                is,
                                GlobalV::NSPIN,
                                pelec->charge->rho_save[is],
                                iter,
                                this->pelec->eferm.get_efval(is),
                                &(GlobalC::ucell),
                                GlobalV::global_out_dir,
                                precision,
                                tag,
                                prefix);
}

//------------------------------------------------------------------------------
//! the 14th function of ESolver_KS: create_Output_Kin
//! mohan add 2024-05-12
//------------------------------------------------------------------------------
template <typename T, typename Device>
ModuleIO::Output_Rho ESolver_KS<T, Device>::create_Output_Kin(int is, int iter, const std::string& prefix)
{
    const int precision = 11;
    std::string tag = "TAU";
    return ModuleIO::Output_Rho(this->pw_big,
                                this->pw_rhod,
                                is,
                                GlobalV::NSPIN,
                                pelec->charge->kin_r_save[is],
                                iter,
                                this->pelec->eferm.get_efval(is),
                                &(GlobalC::ucell),
                                GlobalV::global_out_dir,
                                precision,
                                tag,
                                prefix);
}

//------------------------------------------------------------------------------
//! the 15th function of ESolver_KS: create_Output_Potential
//! mohan add 2024-05-12
//------------------------------------------------------------------------------
template <typename T, typename Device>
ModuleIO::Output_Potential ESolver_KS<T, Device>::create_Output_Potential(int iter, const std::string& prefix)
{
    const int precision = 3;
    std::string tag = "POT";
    return ModuleIO::Output_Potential(this->pw_big,
                                      this->pw_rhod,
                                      GlobalV::NSPIN,
                                      iter,
                                      GlobalV::out_pot,
                                      this->pelec->pot->get_effective_v(),
                                      this->pelec->pot->get_fixed_v(),
                                      &(GlobalC::ucell),
                                      pelec->charge,
                                      precision,
                                      GlobalV::global_out_dir,
                                      tag,
                                      prefix);
}

//------------------------------------------------------------------------------
//! the 16th-20th functions of ESolver_KS
//! mohan add 2024-05-12
//------------------------------------------------------------------------------
//! This is for mixed-precision pw/LCAO basis sets.
template class ESolver_KS<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_KS<std::complex<double>, base_device::DEVICE_CPU>;

//! This is for GPU codes.
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_KS<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_KS<std::complex<double>, base_device::DEVICE_GPU>;
#endif

//! This is for LCAO basis set.
#ifdef __LCAO
template class ESolver_KS<double, base_device::DEVICE_CPU>;
#endif
} // namespace ModuleESolver
