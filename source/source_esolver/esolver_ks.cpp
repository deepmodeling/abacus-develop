#include "esolver_ks.h"

// To setup plane wave for electronic wave functions
#include "pw_setup.h"

#include "source_base/timer.h"
#include "source_base/global_variable.h"
#include "source_pw/module_pwdft/global.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_dftu/dftu.h"

#include "source_cell/cal_atoms_info.h"
#include "source_estate/elecstate_print.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_hsolver/hsolver.h"
#include "source_io/cube_io.h"

// for NSCF calculations of band structures
#include "source_io/nscf_band.h"
// for output log information
#include "source_io/output_log.h"
#include "source_io/print_info.h"
#include "source_io/write_eig_occ.h"
// for jason output information
#include "source_io/json_output/init_info.h"
#include "source_io/json_output/output_info.h"



namespace ModuleESolver
{

template <typename T, typename Device>
ESolver_KS<T, Device>::ESolver_KS()
{
}


template <typename T, typename Device>
ESolver_KS<T, Device>::~ESolver_KS()
{
    delete this->psi;
    delete this->pw_wfc;
    delete this->p_hamilt;
    delete this->p_chgmix;
    this->ppcell.release_memory();
}


template <typename T, typename Device>
void ESolver_KS<T, Device>::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_KS", "before_all_runners");
    //! 1) initialize "before_all_runniers" in ESolver_FP
    ESolver_FP::before_all_runners(ucell, inp);
    
    classname = "ESolver_KS";
    basisname = "";

    scf_thr = PARAM.inp.scf_thr;
    scf_ene_thr = PARAM.inp.scf_ene_thr;
    maxniter = PARAM.inp.scf_nmax;
    niter = maxniter;
    drho = 0.0;

    std::string fft_device = PARAM.inp.device;

    // Fast Fourier Transform
    // LCAO basis doesn't support GPU acceleration on FFT currently
    if(PARAM.inp.basis_type == "lcao")
    {
        fft_device = "cpu";
    }
    std::string fft_precision = PARAM.inp.precision;
#ifdef __ENABLE_FLOAT_FFTW
    if (PARAM.inp.cal_cond && PARAM.inp.esolver_type == "sdft")
    {
        fft_precision = "mixing";
    }
#endif

    pw_wfc = new ModulePW::PW_Basis_K_Big(fft_device, fft_precision);
    ModulePW::PW_Basis_K_Big* tmp = static_cast<ModulePW::PW_Basis_K_Big*>(pw_wfc);

    // should not use INPUT here, mohan 2024-05-12
    tmp->setbxyz(PARAM.inp.bx, PARAM.inp.by, PARAM.inp.bz);

    ///----------------------------------------------------------
    /// charge mixing
    ///----------------------------------------------------------
    p_chgmix = new Charge_Mixing();
    p_chgmix->set_rhopw(this->pw_rho, this->pw_rhod);

    // cell_factor
    this->ppcell.cell_factor = PARAM.inp.cell_factor;


    //! 3) it has been established that
    // xc_func is same for all elements, therefore
    // only the first one if used
    XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
    
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    //! 4) setup the charge mixing parameters
    p_chgmix->set_mixing(PARAM.inp.mixing_mode,
                         PARAM.inp.mixing_beta,
                         PARAM.inp.mixing_ndim,
                         PARAM.inp.mixing_gg0,
                         PARAM.inp.mixing_tau,
                         PARAM.inp.mixing_beta_mag,
                         PARAM.inp.mixing_gg0_mag,
                         PARAM.inp.mixing_gg0_min,
                         PARAM.inp.mixing_angle,
                         PARAM.inp.mixing_dmr,
                         ucell.omega,
                         ucell.tpiba);

    p_chgmix->init_mixing();

    //! 5) ESolver depends on the Symmetry module
    // symmetry analysis should be performed every time the cell is changed
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        ucell.symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    //! 6) Setup the k points according to symmetry.
    this->kv.set(ucell,ucell.symm, PARAM.inp.kpoint_file, PARAM.inp.nspin, ucell.G, ucell.latvec, GlobalV::ofs_running);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

    //! 7) print information
    ModuleIO::setup_parameters(ucell, this->kv);

    //! 8) setup plane wave for electronic wave functions
    ModuleESolver::pw_setup(inp, ucell, *this->pw_rho, this->kv, *this->pw_wfc);

    //! 9) initialize the real-space uniform grid for FFT and parallel
    //! distribution of plane waves
	Pgrid.init(this->pw_rhod->nx,
			this->pw_rhod->ny,
			this->pw_rhod->nz,
			this->pw_rhod->nplane,
			this->pw_rhod->nrxx,
			pw_big->nbz,
			pw_big->bz);

    //! 10) calculate the structure factor
    this->sf.setup_structure_factor(&ucell, Pgrid, this->pw_rhod);
}

//------------------------------------------------------------------------------
//! the 5th function of ESolver_KS: hamilt2rho_single
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename T, typename Device>
void ESolver_KS<T, Device>::hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr)
{
    ModuleBase::timer::tick(this->classname, "hamilt2rho_single");
    // Temporarily, before HSolver is constructed, it should be overrided by
    // LCAO, PW, SDFT and TDDFT.
    // After HSolver is constructed, LCAO, PW, SDFT should delete their own
    // hamilt2rho_single() and use:
    ModuleBase::timer::tick(this->classname, "hamilt2rho_single");
}

template <typename T, typename Device>
void ESolver_KS<T, Device>::hamilt2rho(UnitCell& ucell, const int istep, const int iter, const double ethr)
{
    // 1) use Hamiltonian to obtain charge density
    this->hamilt2rho_single(ucell, istep, iter, diag_ethr);

    // 2) for MPI: STOGROUP? need to rewrite
    //<Temporary> It may be changed when more clever parallel algorithm is
    // put forward.
    // When parallel algorithm for bands are adopted. Density will only be
    // treated in the first group.
    //(Different ranks should have abtained the same, but small differences
    // always exist in practice.)
    // Maybe in the future, density and wavefunctions should use different
    // parallel algorithms, in which they do not occupy all processors, for
    // example wavefunctions uses 20 processors while density uses 10.
    if (PARAM.globalv.ks_run)
    {
        // double drho = this->estate.caldr2();
        // EState should be used after it is constructed.

        drho = p_chgmix->get_drho(&this->chr, PARAM.inp.nelec);
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

                this->hamilt2rho_single(ucell, istep, iter, diag_ethr);

                drho = p_chgmix->get_drho(&this->chr, PARAM.inp.nelec);

                hsolver_error = hsolver::cal_hsolve_error(PARAM.inp.basis_type,
                                                          PARAM.inp.esolver_type,
                                                          diag_ethr,
                                                          PARAM.inp.nelec);
            }
        }
    }
}

template <typename T, typename Device>
void ESolver_KS<T, Device>::runner(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_KS", "runner");
    ModuleBase::timer::tick(this->classname, "runner");

    //----------------------------------------------------------------
    // 1) before_scf (electronic iteration loops)
    //----------------------------------------------------------------
    this->before_scf(ucell, istep);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT SCF");

    //----------------------------------------------------------------
    // 2) SCF iterations
    //----------------------------------------------------------------
    bool conv_esolver = false;
    this->niter = this->maxniter;
    this->diag_ethr = PARAM.inp.pw_diag_thr;
    for (int iter = 1; iter <= this->maxniter; ++iter)
    {
		//----------------------------------------------------------------
		// 3) initialization of SCF iterations
		//----------------------------------------------------------------
		this->iter_init(ucell, istep, iter);

		//----------------------------------------------------------------
        // 4) use Hamiltonian to obtain charge density
		//----------------------------------------------------------------
        this->hamilt2rho(ucell, istep, iter, diag_ethr);

		//----------------------------------------------------------------
        // 5) finish scf iterations
		//----------------------------------------------------------------
        this->iter_finish(ucell, istep, iter, conv_esolver);

		//----------------------------------------------------------------
        // 6) check convergence
		//----------------------------------------------------------------
        if (conv_esolver || this->oscillate_esolver)
        {
            this->niter = iter;
            if (this->oscillate_esolver)
            {
                std::cout << " !! Density oscillation is found, STOP HERE !!" << std::endl;
            }
            break;
        }
    } // end scf iterations

	//----------------------------------------------------------------
	// 7) after scf
	//----------------------------------------------------------------
    this->after_scf(ucell, istep, conv_esolver);

    ModuleBase::timer::tick(this->classname, "runner");
    return;
};

template <typename T, typename Device>
void ESolver_KS<T, Device>::before_scf(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_KS", "before_scf");
    ESolver_FP::before_scf(ucell, istep);
}

template <typename T, typename Device>
void ESolver_KS<T, Device>::iter_init(UnitCell& ucell, const int istep, const int iter)
{
    if(PARAM.inp.esolver_type != "tddft")
    {
        ModuleIO::write_head(GlobalV::ofs_running, istep, iter, this->basisname);
    }

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

    // save input charge density (rho)
    this->chr.save_rho_before_sum_band();
}

template <typename T, typename Device>
void ESolver_KS<T, Device>::iter_finish(UnitCell& ucell, const int istep, int& iter, bool &conv_esolver)
{

	if(iter % PARAM.inp.out_freq_elec == 0)
	{
		//----------------------------------------------------------------
		// 1) print out band gap 
		//----------------------------------------------------------------
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

		//----------------------------------------------------------------
		// 2) print out eigenvalues and occupations
		//----------------------------------------------------------------
		if (PARAM.inp.out_band[0] || iter == PARAM.inp.scf_nmax || conv_esolver)
		{
			ModuleIO::write_eig_iter(this->pelec->ekb,this->pelec->wg,*this->pelec->klist);
		}
	}

	//----------------------------------------------------------------
    // 2) compute magnetization, only for LSDA(spin==2)
	//----------------------------------------------------------------
    ucell.magnet.compute_mag(ucell.omega,
                                       this->chr.nrxx,
                                       this->chr.nxyz,
                                       this->chr.rho,
                                       this->pelec->nelec_spin.data());

	//----------------------------------------------------------------
    // 3) charge mixing 
	//----------------------------------------------------------------
    if (PARAM.globalv.ks_run)
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

        conv_esolver = (drho < this->scf_thr && not_restart_step && is_U_converged);

        // add energy threshold for SCF convergence
        if (this->scf_ene_thr > 0.0)
        {
            // calculate energy of output charge density
            this->update_pot(ucell, istep, iter, conv_esolver);
            this->pelec->cal_energies(2); // 2 means Kohn-Sham functional
            // now, etot_old is the energy of input density, while etot is the energy of output density
            this->pelec->f_en.etot_delta = this->pelec->f_en.etot - this->pelec->f_en.etot_old;
            // output etot_delta
            GlobalV::ofs_running << " DeltaE_womix = " << this->pelec->f_en.etot_delta * ModuleBase::Ry_to_eV << " eV"
                                 << std::endl;
            if (iter > 1 && conv_esolver == 1) // only check when density is converged
            {
                // update the convergence flag
                conv_esolver
                    = (std::abs(this->pelec->f_en.etot_delta * ModuleBase::Ry_to_eV) < this->scf_ene_thr);
            }
        }

        // If drho < hsolver_error in the first iter or drho < scf_thr, we
        // do not change rho.
        if (drho < hsolver_error || conv_esolver || PARAM.inp.calculation == "nscf")
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
                p_chgmix->mix_rho(&this->chr); // update chr->rho by mixing
            }
            if (PARAM.inp.scf_thr_type == 2)
            {
                this->chr.renormalize_rho(); // renormalize rho in R-space would
                                                  // induce a error in K-space
            }
            //----------charge mixing done-----------
        }
    }

#ifdef __MPI
    MPI_Bcast(&drho, 1, MPI_DOUBLE, 0, BP_WORLD);

    // change MPI_DOUBLE to MPI_C_BOOL, mohan 2025-04-13 
    MPI_Bcast(&conv_esolver, 1, MPI_C_BOOL, 0, BP_WORLD);
    MPI_Bcast(this->chr.rho[0], this->pw_rhod->nrxx, MPI_DOUBLE, 0, BP_WORLD);
#endif

	//----------------------------------------------------------------
    // 4) Update potentials (should be done every SF iter)
	//----------------------------------------------------------------
    // Hamilt should be used after it is constructed.
    // this->phamilt->update(conv_esolver);
    this->update_pot(ucell, istep, iter, conv_esolver);

	//----------------------------------------------------------------
    // 5) calculate energies
	//----------------------------------------------------------------
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



	//----------------------------------------------------------------
    // 6) time and meta-GGA 
	//----------------------------------------------------------------
#ifdef __MPI
    double duration = (double)(MPI_Wtime() - iter_time);
#else
    double duration
        = (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - iter_time)).count()
          / static_cast<double>(1e6);
#endif

    // get mtaGGA related parameters
    double dkin = 0.0; // for meta-GGA
    if (XC_Functional::get_ked_flag())
    {
        dkin = p_chgmix->get_dkin(&this->chr, PARAM.inp.nelec);
    }

    // pint energy
    elecstate::print_etot(ucell.magnet, *pelec,conv_esolver, iter, drho, 
    dkin, duration, diag_ethr);


#ifdef __RAPIDJSON
	//----------------------------------------------------------------
    // 7) add Json of scf mag
	//----------------------------------------------------------------
    Json::add_output_scf_mag(ucell.magnet.tot_mag,
                             ucell.magnet.abs_mag,
                             this->pelec->f_en.etot * ModuleBase::Ry_to_eV,
                             this->pelec->f_en.etot_delta * ModuleBase::Ry_to_eV,
                             drho,
                             duration);
#endif //__RAPIDJSON


	//----------------------------------------------------------------
    // 7) SCF restart information 
	//----------------------------------------------------------------
    if (PARAM.inp.mixing_restart > 0 
        && iter == this->p_chgmix->mixing_restart_step - 1 
        && iter != PARAM.inp.scf_nmax)
    {
        this->p_chgmix->mixing_restart_last = iter;
        std::cout << " SCF restart after this step!" << std::endl;
    }

	//----------------------------------------------------------------
    // 8) Iter finish 
	//----------------------------------------------------------------
    ESolver_FP::iter_finish(ucell, istep, iter, conv_esolver);
}

//! Something to do after SCF iterations when SCF is converged or comes to the max iter step.
template <typename T, typename Device>
void ESolver_KS<T, Device>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS", "after_scf");
    
    // 1) calculate the kinetic energy density tau
    if (PARAM.inp.out_elf[0] > 0)
    {
        assert(this->psi != nullptr);
        this->pelec->cal_tau(*(this->psi));
    }
    
    // 2) call after_scf() of ESolver_FP
    ESolver_FP::after_scf(ucell, istep, conv_esolver);

    // 3) write eigenvalues
    if (istep % PARAM.inp.out_interval == 0)
    {
//        elecstate::print_eigenvalue(this->pelec->ekb,this->pelec->wg,this->pelec->klist,GlobalV::ofs_running);
    }
}

template <typename T, typename Device>
void ESolver_KS<T, Device>::after_all_runners(UnitCell& ucell)
{
    // 1) write Etot information
    ESolver_FP::after_all_runners(ucell);

    // 2) write eigenvalue information
    ModuleIO::write_eig_file(this->pelec->ekb, this->pelec->wg, this->kv);

    // 3) write band information
    if (PARAM.inp.out_band[0])
    {
        const int nspin0 = (PARAM.inp.nspin == 2) ? 2 : 1;
        for (int is = 0; is < nspin0; is++)
        {
            std::stringstream ss2;
            ss2 << PARAM.globalv.global_out_dir << "eigs" << is + 1 << ".txt";
            const double eshift = 0.0;
            ModuleIO::nscf_band(is,
                                ss2.str(),
                                PARAM.inp.nbands,
                                eshift,
                                PARAM.inp.out_band[1],
                                this->pelec->ekb,
                                this->kv);
        }
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
