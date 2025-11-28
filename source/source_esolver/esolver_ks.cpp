#include "esolver_ks.h"

// for jason output information
#include "source_io/json_output/init_info.h"
#include "source_io/json_output/output_info.h"

#include "source_estate/update_pot.h" // mohan add 20251016
#include "source_estate/module_charge/chgmixing.h" // mohan add 20251018
#include "source_pw/module_pwdft/setup_pwwfc.h" // mohan add 20251018
#include "source_hsolver/hsolver.h"
#include "source_io/write_eig_occ.h"
#include "source_io/write_bands.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/output_log.h" // use write_head
#include "source_estate/elecstate_print.h" // print_etot
#include "source_io/print_info.h" // print_parameters
#include "source_psi/setup_psi.h" // mohan add 20251009
#include "source_lcao/module_dftu/dftu.h" // mohan add 2025-11-07

namespace ModuleESolver
{

template <typename T, typename Device>
ESolver_KS<T, Device>::ESolver_KS(){}


template <typename T, typename Device>
ESolver_KS<T, Device>::~ESolver_KS()
{
	//****************************************************
	// do not add any codes in this deconstructor funcion
	//****************************************************
    Setup_Psi<T>::deallocate_psi(this->psi);

    delete this->p_hamilt;
    delete this->p_chgmix;
    this->ppcell.release_memory();
    
    // mohan add 2025-10-18, should be put int clean() function
    pw::teardown_pwwfc(this->pw_wfc);
}


template <typename T, typename Device>
void ESolver_KS<T, Device>::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_KS", "before_all_runners");

    //! 1) setup "before_all_runniers" in ESolver_FP
    ESolver_FP::before_all_runners(ucell, inp);
    
    //! 2) setup some parameters
    classname = "ESolver_KS";
    basisname = "";

    this->scf_thr = inp.scf_thr;
    this->scf_ene_thr = inp.scf_ene_thr;
    this->maxniter = inp.scf_nmax;
    this->niter = maxniter;
    this->drho = 0.0;

    // cell_factor
    this->ppcell.cell_factor = inp.cell_factor;

    //! 3) setup charge mixing
    p_chgmix = new Charge_Mixing();
    p_chgmix->set_rhopw(this->pw_rho, this->pw_rhod);

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    //! 4) setup Exc for the first element '0' (all elements have same exc) 
    XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
    GlobalV::ofs_running<<XC_Functional::output_info()<<std::endl;
    
    //! 5) setup the charge mixing parameters
    p_chgmix->set_mixing(inp.mixing_mode, inp.mixing_beta, inp.mixing_ndim,
      inp.mixing_gg0, inp.mixing_tau, inp.mixing_beta_mag, inp.mixing_gg0_mag,
      inp.mixing_gg0_min, inp.mixing_angle, inp.mixing_dmr, ucell.omega, ucell.tpiba);

    p_chgmix->init_mixing();

    //! 6) symmetry analysis should be performed every time the cell is changed
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        ucell.symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    //! 7) setup k points in the Brillouin zone according to symmetry.
    this->kv.set(ucell,ucell.symm, inp.kpoint_file, inp.nspin, ucell.G, ucell.latvec, GlobalV::ofs_running);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

    //! 8) print information
    ModuleIO::print_parameters(ucell, this->kv, inp);

    //! 9) setup plane wave for electronic wave functions
    pw::setup_pwwfc(inp, ucell, *this->pw_rho, this->kv, this->pw_wfc);

    //! 10) parallel of FFT grid 
	Pgrid.init(this->pw_rhod->nx, this->pw_rhod->ny, this->pw_rhod->nz,
			this->pw_rhod->nplane, this->pw_rhod->nrxx, pw_big->nbz, pw_big->bz);

    //! 11) calculate the structure factor
    this->sf.setup(&ucell, Pgrid, this->pw_rhod);
}

template <typename T, typename Device>
void ESolver_KS<T, Device>::hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr)
{}

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
    this->scf_nmax_flag = false; // mohan add 2025-09-21
    for (int iter = 1; iter <= this->maxniter; ++iter)
	{
        // mohan add 2025-09-21
		if(iter == this->maxniter)
		{
			this->scf_nmax_flag=true;
		}

		// 3) initialization of SCF iterations
		this->iter_init(ucell, istep, iter);

        // 4) use Hamiltonian to obtain charge density
        this->hamilt2rho(ucell, istep, iter, diag_ethr);

        // 5) finish scf iterations
        this->iter_finish(ucell, istep, iter, conv_esolver);

        // 6) check convergence
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

	// 7) after scf
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
        diag_ethr = hsolver::set_diagethr_ks(PARAM.inp.basis_type, PARAM.inp.esolver_type,
          PARAM.inp.calculation, PARAM.inp.init_chg, PARAM.inp.precision, istep, iter,
          drho, PARAM.inp.pw_diag_thr, diag_ethr, PARAM.inp.nelec);
    }
    else if (PARAM.inp.esolver_type == "sdft")
    {
        diag_ethr = hsolver::set_diagethr_sdft(PARAM.inp.basis_type, PARAM.inp.esolver_type,
          PARAM.inp.calculation, PARAM.inp.init_chg, istep, iter, drho,
          PARAM.inp.pw_diag_thr, diag_ethr, PARAM.inp.nbands, esolver_KS_ne);
    }

    // save input charge density (rho)
    this->chr.save_rho_before_sum_band();
}

template <typename T, typename Device>
void ESolver_KS<T, Device>::iter_finish(UnitCell& ucell, const int istep, int& iter, bool &conv_esolver)
{

    // 1.1) print out band gap 
    if (!PARAM.globalv.two_fermi)
    {
        this->pelec->cal_bandgap();
    }
    else
    {
        this->pelec->cal_bandgap_updw();
    }

    // 1.2) print out eigenvalues and occupations
    if (PARAM.inp.out_band[0])
    {
		if (iter % PARAM.inp.out_freq_elec == 0 || iter == PARAM.inp.scf_nmax || conv_esolver)
		{
			ModuleIO::write_eig_iter(this->pelec->ekb,this->pelec->wg,*this->pelec->klist);
		}
    }

    // 2.1) compute magnetization, only for spin==2
    ucell.magnet.compute_mag(ucell.omega, this->chr.nrxx, this->chr.nxyz, this->chr.rho,
                                       this->pelec->nelec_spin.data());

    // 2.2) charge mixing 
    // SCF will continue if U is not converged for uramping calculation
	bool converged_u = true;
	// to avoid unnecessary dependence on dft+u, refactor is needed
#ifdef __LCAO
	if (PARAM.inp.dft_plus_u)
	{
		converged_u = this->dftu.u_converged();
	}
#endif

    module_charge::chgmixing_ks(iter, ucell, this->pelec, this->chr, this->p_chgmix, 
      this->pw_rhod->nrxx, this->drho, this->oscillate_esolver, conv_esolver, hsolver_error, 
      this->scf_thr, this->scf_ene_thr, converged_u, PARAM.inp);

    // 2.3) Update potentials (should be done every SF iter)
    elecstate::update_pot(ucell, this->pelec, this->chr, conv_esolver);

    // 3.1) calculate energies
    this->pelec->cal_energies(1); // Harris-Foulkes functional
    this->pelec->cal_energies(2); // Kohn-Sham functional

    if (iter == 1)
    {
        this->pelec->f_en.etot_old = this->pelec->f_en.etot;
    }
    this->pelec->f_en.etot_delta = this->pelec->f_en.etot - this->pelec->f_en.etot_old;
    this->pelec->f_en.etot_old = this->pelec->f_en.etot;

    // 4) get meta-GGA related parameters
    double dkin = 0.0; // for meta-GGA
    if (XC_Functional::get_ked_flag())
    {
        dkin = p_chgmix->get_dkin(&this->chr, PARAM.inp.nelec);
    }

    // Iter finish 
    ESolver_FP::iter_finish(ucell, istep, iter, conv_esolver);


    // the end, print time
#ifdef __MPI
    double duration = (double)(MPI_Wtime() - iter_time);
#else
    double duration
        = (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - iter_time)).count()
          / static_cast<double>(1e6);
#endif

    // print energies
    elecstate::print_etot(ucell.magnet, *pelec, conv_esolver, iter, drho, 
    dkin, duration, diag_ethr);


#ifdef __RAPIDJSON
    // add Json of scf mag
    Json::add_output_scf_mag(ucell.magnet.tot_mag, ucell.magnet.abs_mag,
                             this->pelec->f_en.etot * ModuleBase::Ry_to_eV,
                             this->pelec->f_en.etot_delta * ModuleBase::Ry_to_eV,
                             drho, duration);
#endif //__RAPIDJSON

}

//! Something to do after SCF iterations when SCF is converged or comes to the max iter step.
template <typename T, typename Device>
void ESolver_KS<T, Device>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS", "after_scf");
    
/*
    // 1) calculate the kinetic energy density tau
    if (PARAM.inp.out_elf[0] > 0)
    {
        assert(this->psi != nullptr);
        this->pelec->cal_tau(*(this->psi));
    }
*/
    
    // 2) call after_scf() of ESolver_FP
    ESolver_FP::after_scf(ucell, istep, conv_esolver);

    // 3) write eigenvalues and occupations to eig_occ.txt
    ModuleIO::write_eig_file(this->pelec->ekb, this->pelec->wg, this->kv, istep);

    // 4) write band information to band.txt
    ModuleIO::write_bands(PARAM.inp, this->pelec->ekb, this->kv);

}

template <typename T, typename Device>
void ESolver_KS<T, Device>::after_all_runners(UnitCell& ucell)
{
    // 1) write Etot information
    ESolver_FP::after_all_runners(ucell);
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
