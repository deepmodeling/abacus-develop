#include "esolver_ks.h"

#include "module_base/timer.h"
#include "module_cell/cal_atoms_info.h"
#include "module_io/json_output/init_info.h"
#include "module_io/json_output/output_info.h"
#include "module_io/output_log.h"
#include "module_io/print_info.h"
#include "module_io/write_istate_info.h"
#include "module_parameter/parameter.h"

#include <ctime>
#include <iostream>
//--------------Temporary----------------
#include "module_base/global_variable.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
//---------------------------------------
#ifdef USE_PAW
#include "module_base/parallel_common.h"
#include "module_cell/module_paw/paw_cell.h"
#endif

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
    scf_thr = PARAM.inp.scf_thr;
    scf_ene_thr = PARAM.inp.scf_ene_thr;
    drho = 0.0;

    // should not use GlobalV here, mohan 2024-05-12
    maxniter = PARAM.inp.scf_nmax;
    niter = maxniter;

    // should not use GlobalV here, mohan 2024-05-12
    out_freq_elec = PARAM.inp.out_freq_elec;

    // pw_rho = new ModuleBase::PW_Basis();
    // temporary, it will be removed
    std::string fft_device = PARAM.inp.device;
    // LCAO basis doesn't support GPU acceleration on FFT currently
    if(PARAM.inp.basis_type == "lcao")
    {
        fft_device = "cpu";
    }
    pw_wfc = new ModulePW::PW_Basis_K_Big(fft_device, PARAM.inp.precision);
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
    p_chgmix->set_mixing(PARAM.inp.mixing_mode,
                         PARAM.inp.mixing_beta,
                         PARAM.inp.mixing_ndim,
                         PARAM.inp.mixing_gg0,
                         PARAM.inp.mixing_tau,
                         PARAM.inp.mixing_beta_mag,
                         PARAM.inp.mixing_gg0_mag,
                         PARAM.inp.mixing_gg0_min,
                         PARAM.inp.mixing_angle,
                         PARAM.inp.mixing_dmr);

    /// PAW Section
#ifdef USE_PAW
    if (PARAM.inp.use_paw)
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
            std::ifstream ifa(PARAM.globalv.global_in_stru.c_str(), std::ios::in);
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
        CalAtomsInfo ca;
        ca.cal_atoms_info(ucell.atoms, ucell.ntype, PARAM);
    }
#endif
    /// End PAW

    //! 4) it has been established that
    // xc_func is same for all elements, therefore
    // only the first one if used
    if (PARAM.inp.use_paw)
    {
        XC_Functional::set_xc_type(PARAM.inp.dft_functional);
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

    //! 6) Setup the k points according to symmetry.
    this->kv.set(ucell.symm, PARAM.inp.kpoint_file, PARAM.inp.nspin, ucell.G, ucell.latvec, GlobalV::ofs_running);

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

    //! 7) print information
    ModuleIO::setup_parameters(ucell, this->kv);

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

    this->pw_wfc->fft_bundle.initfftmode(inp.fft_mode);
    this->pw_wfc->setuptransform();

    //! 9) initialize the number of plane waves for each k point
    for (int ik = 0; ik < this->kv.get_nks(); ++ik)
    {
        this->kv.ngk[ik] = this->pw_wfc->npwk[ik];
    }

    this->pw_wfc->collect_local_pw(inp.erf_ecut, inp.erf_height, inp.erf_sigma);

    ModuleIO::print_wfcfft(inp, *this->pw_wfc, GlobalV::ofs_running);

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

#ifdef USE_PAW
    if (PARAM.inp.use_paw)
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
//! the 5th function of ESolver_KS: hamilt2density_single
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename T, typename Device>
void ESolver_KS<T, Device>::hamilt2density_single(const int istep, const int iter, const double ethr)
{
    ModuleBase::timer::tick(this->classname, "hamilt2density_single");
    // Temporarily, before HSolver is constructed, it should be overrided by
    // LCAO, PW, SDFT and TDDFT.
    // After HSolver is constructed, LCAO, PW, SDFT should delete their own
    // hamilt2density_single() and use:
    ModuleBase::timer::tick(this->classname, "hamilt2density_single");
}

template <typename T, typename Device>
void ESolver_KS<T, Device>::hamilt2density(const int istep, const int iter, const double ethr)
{
    // 7) use Hamiltonian to obtain charge density
    this->hamilt2density_single(istep, iter, diag_ethr);

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

        drho = p_chgmix->get_drho(pelec->charge, PARAM.inp.nelec);
        hsolver_error = 0.0;
        if (iter == 1 && PARAM.inp.calculation != "nscf")
        {
            hsolver_error
                = hsolver::cal_hsolve_error(PARAM.inp.basis_type, PARAM.inp.esolver_type, diag_ethr, PARAM.inp.nelec);

            // The error of HSolver is larger than drho,
            // so a more precise HSolver should be executed.
            if (hsolver_error > drho)
            {
                diag_ethr = hsolver::reset_diag_ethr(GlobalV::ofs_running,
                                                     PARAM.inp.basis_type,
                                                     PARAM.inp.esolver_type,
                                                     PARAM.inp.precision,
                                                     hsolver_error,
                                                     drho,
                                                     diag_ethr,
                                                     PARAM.inp.nelec);

                this->hamilt2density_single(istep, iter, diag_ethr);

                drho = p_chgmix->get_drho(pelec->charge, PARAM.inp.nelec);

                hsolver_error = hsolver::cal_hsolve_error(PARAM.inp.basis_type,
                                                          PARAM.inp.esolver_type,
                                                          diag_ethr,
                                                          PARAM.inp.nelec);
            }
        }
    }
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
    ModuleBase::timer::tick(this->classname, "before_scf");
    this->before_scf(istep);
    ModuleBase::timer::tick(this->classname, "before_scf");

    // 3) write charge density
    if (PARAM.inp.dm_to_rho)
    {
        ModuleBase::timer::tick(this->classname, "runner");
        return; // nothing further is needed
    }

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT SCF");

    // 4) SCF iterations
    this->conv_esolver = false;
    this->niter = this->maxniter;
    this->diag_ethr = PARAM.inp.pw_diag_thr;
    for (int iter = 1; iter <= this->maxniter; ++iter)
    {
        // 6) initialization of SCF iterations
        this->iter_init(istep, iter);

        this->hamilt2density(istep, iter, diag_ethr);

        // 10) finish scf iterations
        this->iter_finish(istep, iter);

        // 13) check convergence
        if (this->conv_esolver || this->oscillate_esolver)
        {
            this->niter = iter;
            if (this->oscillate_esolver)
            {
                std::cout << " !! Density oscillation is found, STOP HERE !!" << std::endl;
            }
            break;
        }
    } // end scf iterations

    // 15) after scf
    ModuleBase::timer::tick(this->classname, "after_scf");
    this->after_scf(istep);
    ModuleBase::timer::tick(this->classname, "after_scf");

    ModuleBase::timer::tick(this->classname, "runner");
    return;
};

template <typename T, typename Device>
void ESolver_KS<T, Device>::iter_init(const int istep, const int iter)
{
    ModuleIO::write_head(GlobalV::ofs_running, istep, iter, this->basisname);

#ifdef __MPI
    iter_time = MPI_Wtime();
#else
    iter_time = std::chrono::system_clock::now();
#endif

    if (PARAM.inp.esolver_type == "ksdft")
    {
        diag_ethr = hsolver::set_diagethr_ks(PARAM.inp.basis_type,
                                             PARAM.inp.esolver_type,
                                             PARAM.inp.calculation,
                                             PARAM.inp.init_chg,
                                             PARAM.inp.precision,
                                             istep,
                                             iter,
                                             drho,
                                             PARAM.inp.pw_diag_thr,
                                             diag_ethr,
                                             PARAM.inp.nelec);
    }
    else if (PARAM.inp.esolver_type == "sdft")
    {
        diag_ethr = hsolver::set_diagethr_sdft(PARAM.inp.basis_type,
                                               PARAM.inp.esolver_type,
                                               PARAM.inp.calculation,
                                               PARAM.inp.init_chg,
                                               istep,
                                               iter,
                                               drho,
                                               PARAM.inp.pw_diag_thr,
                                               diag_ethr,
                                               PARAM.inp.nbands,
                                               esolver_KS_ne);
    }

    // 1) save input rho
    this->pelec->charge->save_rho_before_sum_band();
}

template <typename T, typename Device>
void ESolver_KS<T, Device>::iter_finish(const int istep, int& iter)
{
    if (PARAM.inp.out_bandgap)
    {
        if (!PARAM.globalv.two_fermi)
        {
            this->pelec->cal_bandgap();
        }
        else
        {
            this->pelec->cal_bandgap_updw();
        }
    }

    for (int ik = 0; ik < this->kv.get_nks(); ++ik)
    {
        this->pelec->print_band(ik, PARAM.inp.printe, iter);
    }

    // compute magnetization, only for LSDA(spin==2)
    GlobalC::ucell.magnet.compute_magnetization(this->pelec->charge->nrxx,
                                                this->pelec->charge->nxyz,
                                                this->pelec->charge->rho,
                                                this->pelec->nelec_spin.data());

    if (GlobalV::MY_STOGROUP == 0)
    {
        // mixing will restart at this->p_chgmix->mixing_restart steps
        if (drho <= PARAM.inp.mixing_restart && PARAM.inp.mixing_restart > 0.0
            && this->p_chgmix->mixing_restart_step > iter)
        {
            this->p_chgmix->mixing_restart_step = iter + 1;
        }

        if (PARAM.inp.scf_os_stop) // if oscillation is detected, SCF will stop
        {
            this->oscillate_esolver
                = this->p_chgmix->if_scf_oscillate(iter, drho, PARAM.inp.scf_os_ndim, PARAM.inp.scf_os_thr);
        }

        // drho will be 0 at this->p_chgmix->mixing_restart step, which is
        // not ground state
        bool not_restart_step = !(iter == this->p_chgmix->mixing_restart_step && PARAM.inp.mixing_restart > 0.0);
        // SCF will continue if U is not converged for uramping calculation
        bool is_U_converged = true;
        // to avoid unnecessary dependence on dft+u, refactor is needed
#ifdef __LCAO
        if (PARAM.inp.dft_plus_u)
        {
            is_U_converged = GlobalC::dftu.u_converged();
        }
#endif

        this->conv_esolver = (drho < this->scf_thr && not_restart_step && is_U_converged);

        // add energy threshold for SCF convergence
        if (this->scf_ene_thr > 0.0)
        {
            // calculate energy of output charge density
            this->update_pot(istep, iter);
            this->pelec->cal_energies(2); // 2 means Kohn-Sham functional
            // now, etot_old is the energy of input density, while etot is the energy of output density
            this->pelec->f_en.etot_delta = this->pelec->f_en.etot - this->pelec->f_en.etot_old;
            // output etot_delta
            GlobalV::ofs_running << " DeltaE_womix = " << this->pelec->f_en.etot_delta * ModuleBase::Ry_to_eV << " eV"
                                 << std::endl;
            if (iter > 1 && this->conv_esolver == 1) // only check when density is converged
            {
                // update the convergence flag
                this->conv_esolver
                    = (std::abs(this->pelec->f_en.etot_delta * ModuleBase::Ry_to_eV) < this->scf_ene_thr);
            }
        }

        // If drho < hsolver_error in the first iter or drho < scf_thr, we
        // do not change rho.
        if (drho < hsolver_error || this->conv_esolver || PARAM.inp.calculation == "nscf")
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
            if (PARAM.inp.mixing_restart > 0 && iter == this->p_chgmix->mixing_restart_step - 1
                && drho <= PARAM.inp.mixing_restart)
            {
                // do not mix charge density
            }
            else
            {
                p_chgmix->mix_rho(pelec->charge); // update chr->rho by mixing
            }
            if (PARAM.inp.scf_thr_type == 2)
            {
                pelec->charge->renormalize_rho(); // renormalize rho in R-space would
                                                  // induce a error in K-space
            }
            //----------charge mixing done-----------
        }
    }

#ifdef __MPI
    MPI_Bcast(&drho, 1, MPI_DOUBLE, 0, PARAPW_WORLD);
    MPI_Bcast(&this->conv_esolver, 1, MPI_DOUBLE, 0, PARAPW_WORLD);
    MPI_Bcast(pelec->charge->rho[0], this->pw_rhod->nrxx, MPI_DOUBLE, 0, PARAPW_WORLD);
#endif

    // update potential
    // Hamilt should be used after it is constructed.
    // this->phamilt->update(conv_esolver);
    this->update_pot(istep, iter);

    // 1 means Harris-Foulkes functional
    // 2 means Kohn-Sham functional
    this->pelec->cal_energies(1);
    this->pelec->cal_energies(2);

    if (iter == 1)
    {
        this->pelec->f_en.etot_old = this->pelec->f_en.etot;
    }
    this->pelec->f_en.etot_delta = this->pelec->f_en.etot - this->pelec->f_en.etot_old;
    this->pelec->f_en.etot_old = this->pelec->f_en.etot;

#ifdef __MPI
    double duration = (double)(MPI_Wtime() - iter_time);
#else
    double duration
        = (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - iter_time)).count()
          / static_cast<double>(1e6);
#endif

    // get mtaGGA related parameters
    double dkin = 0.0; // for meta-GGA
    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
        dkin = p_chgmix->get_dkin(pelec->charge, PARAM.inp.nelec);
    }
    this->pelec->print_etot(this->conv_esolver, iter, drho, dkin, duration, PARAM.inp.printe, diag_ethr);

    // Json, need to be moved to somewhere else
#ifdef __RAPIDJSON
    // add Json of scf mag
    Json::add_output_scf_mag(GlobalC::ucell.magnet.tot_magnetization,
                             GlobalC::ucell.magnet.abs_magnetization,
                             this->pelec->f_en.etot * ModuleBase::Ry_to_eV,
                             this->pelec->f_en.etot_delta * ModuleBase::Ry_to_eV,
                             drho,
                             duration);
#endif //__RAPIDJSON

    // notice for restart
    if (PARAM.inp.mixing_restart > 0 && iter == this->p_chgmix->mixing_restart_step - 1 && iter != PARAM.inp.scf_nmax)
    {
        this->p_chgmix->mixing_restart_last = iter;
        std::cout << " SCF restart after this step!" << std::endl;
    }
}

//! Something to do after SCF iterations when SCF is converged or comes to the max iter step.
template <typename T, typename Device>
void ESolver_KS<T, Device>::after_scf(const int istep)
{
    // 1) call after_scf() of ESolver_FP
    ESolver_FP::after_scf(istep);

    // 2) write eigenvalues
    if (istep % PARAM.inp.out_interval == 0)
    {
        this->pelec->print_eigenvalue(GlobalV::ofs_running);
    }
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
