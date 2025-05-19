#include "esolver_ks_pw.h"

#include "module_base/formatter.h"
#include "module_base/global_variable.h"
#include "module_base/kernels/math_kernel_op.h"
#include "module_base/memory.h"
#include "module_elecstate/cal_ux.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_elecstate/elecstate_pw_sdft.h"
#include "module_elecstate/elecstate_tools.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_general/module_vdw/vdw.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/elecond.h"
#include "module_hamilt_pw/hamilt_pwdft/forces.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/onsite_projector.h"
#include "module_hamilt_pw/hamilt_pwdft/stress_pw.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/hsolver_pw.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include "module_io/berryphase.h"
#include "module_io/cal_ldos.h"
#include "module_io/get_pchg_pw.h"
#include "module_io/get_wf_pw.h"
#include "module_io/numerical_basis.h"
#include "module_io/numerical_descriptor.h"
#include "module_io/to_wannier90_pw.h"
#include "module_io/winput.h"
#include "module_io/write_dos_pw.h"
#include "module_io/write_istate_info.h"
#include "module_io/write_wfc_pw.h"
#include "module_io/write_wfc_r.h"
#include "module_parameter/parameter.h"

#include <iostream>

#ifdef __MLKEDF
#include "module_hamilt_pw/hamilt_ofdft/ml_data.h"
#endif

#include <ATen/kernels/blas.h>
#include <ATen/kernels/lapack.h>

#ifdef __DSP
#include "module_base/kernels/dsp/dsp_connector.h"
#endif

#include <chrono>

namespace ModuleESolver
{

template <typename T, typename Device>
ESolver_KS_PW<T, Device>::ESolver_KS_PW()
{
    this->classname = "ESolver_KS_PW";
    this->basisname = "PW";
    this->device = base_device::get_device_type<Device>(this->ctx);

#if ((defined __CUDA) || (defined __ROCM))
    if (this->device == base_device::GpuDevice)
    {
        ModuleBase::createGpuBlasHandle();
        hsolver::createGpuSolverHandle();
        container::kernels::createGpuBlasHandle();
        container::kernels::createGpuSolverHandle();
    }
#endif

#ifdef __DSP
    std::cout << " ** Initializing DSP Hardware..." << std::endl;
    mtfunc::dspInitHandle(GlobalV::MY_RANK);
#endif
}

template <typename T, typename Device>
ESolver_KS_PW<T, Device>::~ESolver_KS_PW()
{

    // delete Hamilt
    this->deallocate_hamilt();

    if (this->pelec != nullptr)
    {
        delete reinterpret_cast<elecstate::ElecStatePW<T, Device>*>(this->pelec);
        this->pelec = nullptr;
    }

    if (this->device == base_device::GpuDevice)
    {
#if defined(__CUDA) || defined(__ROCM)
        ModuleBase::destoryBLAShandle();
        hsolver::destroyGpuSolverHandle();
        container::kernels::destroyGpuBlasHandle();
        container::kernels::destroyGpuSolverHandle();
#endif
    }

#ifdef __DSP
    std::cout << " ** Closing DSP Hardware..." << std::endl;
    mtfunc::dspDestoryHandle(GlobalV::MY_RANK);
#endif

    if (PARAM.inp.device == "gpu" || PARAM.inp.precision == "single")
    {
        delete this->kspw_psi;
    }
    if (PARAM.inp.precision == "single")
    {
        delete this->__kspw_psi;
    }

    delete this->psi;
    delete this->p_psi_init;
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::allocate_hamilt(const UnitCell& ucell)
{
    this->p_hamilt = new hamilt::HamiltPW<T, Device>(this->pelec->pot, this->pw_wfc, &this->kv, &this->ppcell, &ucell);
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
    // 1) call before_all_runners() of ESolver_KS
    ESolver_KS<T, Device>::before_all_runners(ucell, inp);

    // 2) initialize ElecState, set pelec pointer
    if (this->pelec == nullptr)
    {
        if (inp.esolver_type == "sdft")
        {
            //! SDFT only supports double precision currently
            this->pelec = new elecstate::ElecStatePW_SDFT<std::complex<double>, Device>(this->pw_wfc,
                                                                                        &(this->chr),
                                                                                        &(this->kv),
                                                                                        &ucell,
                                                                                        &(this->ppcell),
                                                                                        this->pw_rhod,
                                                                                        this->pw_rho,
                                                                                        this->pw_big);
        }
        else
        {
            this->pelec = new elecstate::ElecStatePW<T, Device>(this->pw_wfc,
                                                                &(this->chr),
                                                                &(this->kv),
                                                                &ucell,
                                                                &this->ppcell,
                                                                this->pw_rhod,
                                                                this->pw_rho,
                                                                this->pw_big);
        }
    }

    //! set the cell volume variable in pelec
    this->pelec->omega = ucell.omega;

    //! 3) inititlize the charge density.
    this->chr.allocate(PARAM.inp.nspin);

    //! 4) initialize the potential.
    if (this->pelec->pot == nullptr)
    {
        this->pelec->pot = new elecstate::Potential(this->pw_rhod,
                                                    this->pw_rho,
                                                    &ucell,
                                                    &this->locpp.vloc,
                                                    &(this->sf),
                                                    &(this->solvent),
                                                    &(this->pelec->f_en.etxc),
                                                    &(this->pelec->f_en.vtxc));
    }

    //! 5) Initalize local pseudopotential
    this->locpp.init_vloc(ucell, this->pw_rhod);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    //! 6) Initalize non-local pseudopotential
    this->ppcell.init(ucell, &this->sf, this->pw_wfc);
    this->ppcell.init_vnl(ucell, this->pw_rhod);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

    //! 7) Allocate and initialize psi
    this->p_psi_init = new psi::PSIInit<T, Device>(PARAM.inp.init_wfc,
                                                   PARAM.inp.ks_solver,
                                                   PARAM.inp.basis_type,
                                                   GlobalV::MY_RANK,
                                                   ucell,
                                                   this->sf,
                                                   this->kv,
                                                   this->ppcell,
                                                   *this->pw_wfc);

    allocate_psi(this->psi, this->kv.get_nks(), this->kv.ngk, PARAM.globalv.nbands_l, this->pw_wfc->npwk_max);

    this->p_psi_init->prepare_init(PARAM.inp.pw_seed);

    this->kspw_psi = PARAM.inp.device == "gpu" || PARAM.inp.precision == "single"
                         ? new psi::Psi<T, Device>(this->psi[0])
                         : reinterpret_cast<psi::Psi<T, Device>*>(this->psi);

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT BASIS");

    //! 8) setup occupations
    if (PARAM.inp.ocp)
    {
        elecstate::fixed_weights(PARAM.inp.ocp_kb,
                                 PARAM.inp.nbands,
                                 PARAM.inp.nelec,
                                 this->pelec->klist,
                                 this->pelec->wg,
                                 this->pelec->skip_weights);
    }

    // 9) initialize exx pw
    if (PARAM.inp.calculation == "scf" || PARAM.inp.calculation == "relax" || PARAM.inp.calculation == "cell-relax"
        || PARAM.inp.calculation == "md")
    {
        if (GlobalC::exx_info.info_global.cal_exx && GlobalC::exx_info.info_global.separate_loop == true)
        {
            XC_Functional::set_xc_first_loop(ucell);
            exx_helper.set_firstiter();
        }

        if (GlobalC::exx_info.info_global.cal_exx)
        {
            exx_helper.set_wg(&this->pelec->wg);
        }
    }
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::before_scf(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_KS_PW", "before_scf");
    ModuleBase::timer::tick("ESolver_KS_PW", "before_scf");

    //----------------------------------------------------------
    //! 1) Call before_scf() of ESolver_KS
    //----------------------------------------------------------
    ESolver_KS<T, Device>::before_scf(ucell, istep);

    //----------------------------------------------------------
    //! 2) Init variables (cell changed)
    //----------------------------------------------------------
    if (ucell.cell_parameter_updated)
    {
        this->ppcell.rescale_vnl(ucell.omega);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

        this->pw_wfc->initgrids(ucell.lat0, ucell.latvec, this->pw_wfc->nx, this->pw_wfc->ny, this->pw_wfc->nz);

        this->pw_wfc->initparameters(false, PARAM.inp.ecutwfc, this->kv.get_nks(), this->kv.kvec_d.data());

        this->pw_wfc->collect_local_pw(PARAM.inp.erf_ecut, PARAM.inp.erf_height, PARAM.inp.erf_sigma);

        this->p_psi_init->prepare_init(PARAM.inp.pw_seed);
    }

    //----------------------------------------------------------
    //! 3) init Hamiltonian (cell changed)
    //----------------------------------------------------------
    //! Operators in HamiltPW should be reallocated once cell changed
    //! delete Hamilt if not first scf
    this->deallocate_hamilt();

    // allocate HamiltPW
    this->allocate_hamilt(ucell);

    //----------------------------------------------------------
    //! 4) Exx calculations
    //----------------------------------------------------------
    if (PARAM.inp.calculation == "scf" || PARAM.inp.calculation == "relax" || PARAM.inp.calculation == "cell-relax"
        || PARAM.inp.calculation == "md")
    {
        if (GlobalC::exx_info.info_global.cal_exx && PARAM.inp.basis_type == "pw")
        {
            auto hamilt_pw = reinterpret_cast<hamilt::HamiltPW<T, Device>*>(this->p_hamilt);
            hamilt_pw->set_exx_helper(exx_helper);
        }
    }

    //----------------------------------------------------------
    //! 5) Renew local pseudopotential
    //----------------------------------------------------------
    this->pelec
        ->init_scf(istep, ucell, this->Pgrid, this->sf.strucFac, this->locpp.numeric, ucell.symm, (void*)this->pw_wfc);

    //----------------------------------------------------------
    //! 6) Symmetrize the charge density (rho)
    //----------------------------------------------------------

    //! Symmetry_rho should behind init_scf, because charge should be
    //! initialized first. liuyu comment: Symmetry_rho should be
    //! located between init_rho and v_of_rho?
    Symmetry_rho srho;
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        srho.begin(is, this->chr, this->pw_rhod, ucell.symm);
    }

    //----------------------------------------------------------
    //! 7) Calculate the effective potential with rho
    //----------------------------------------------------------
    //! liuyu move here 2023-10-09
    //! D in uspp need vloc, thus behind init_scf()
    //! calculate the effective coefficient matrix
    //! for non-local pseudopotential projectors
    ModuleBase::matrix veff = this->pelec->pot->get_effective_v();

    this->ppcell.cal_effective_D(veff, this->pw_rhod, ucell);

    //----------------------------------------------------------
    //! 8) Onsite projectors
    //----------------------------------------------------------
    if (PARAM.inp.onsite_radius > 0)
    {
        auto* onsite_p = projectors::OnsiteProjector<double, Device>::get_instance();
        onsite_p->init(PARAM.inp.orbital_dir,
                       &ucell,
                       *(this->kspw_psi),
                       this->kv,
                       *(this->pw_wfc),
                       this->sf,
                       PARAM.inp.onsite_radius,
                       PARAM.globalv.nqx,
                       PARAM.globalv.dq,
                       this->pelec->wg,
                       this->pelec->ekb);
    }

    //----------------------------------------------------------
    //! 9) Spin-constrained algorithms
    //----------------------------------------------------------
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<std::complex<double>>& sc
            = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
        sc.init_sc(PARAM.inp.sc_thr,
                   PARAM.inp.nsc,
                   PARAM.inp.nsc_min,
                   PARAM.inp.alpha_trial,
                   PARAM.inp.sccut,
                   PARAM.inp.sc_drop_thr,
                   ucell,
                   nullptr,
                   PARAM.inp.nspin,
                   this->kv,
                   this->p_hamilt,
                   this->kspw_psi,
                   this->pelec,
                   this->pw_wfc);
    }

    //----------------------------------------------------------
    //! 10) DFT+U algorithm
    //----------------------------------------------------------
    if (PARAM.inp.dft_plus_u)
    {
        auto* dftu = ModuleDFTU::DFTU::get_instance();
        dftu->init(ucell, nullptr, this->kv.get_nks());
    }

    //----------------------------------------------------------
    //! 10) Initialize wave functions
    //----------------------------------------------------------
    if (!this->already_initpsi)
    {
        this->p_psi_init->initialize_psi(this->psi, this->kspw_psi, this->p_hamilt, GlobalV::ofs_running);
        this->already_initpsi = true;
    }

    //----------------------------------------------------------
    //! 11) Exx calculations
    //----------------------------------------------------------
    if (PARAM.inp.calculation == "scf" || PARAM.inp.calculation == "relax" || PARAM.inp.calculation == "cell-relax"
        || PARAM.inp.calculation == "md")
    {
        if (GlobalC::exx_info.info_global.cal_exx && PARAM.inp.basis_type == "pw")
        {
            exx_helper.set_psi(kspw_psi);
        }
    }

    ModuleBase::timer::tick("ESolver_KS_PW", "before_scf");
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::iter_init(UnitCell& ucell, const int istep, const int iter)
{
    // call iter_init() of ESolver_KS
    ESolver_KS<T, Device>::iter_init(ucell, istep, iter);

    if (iter == 1)
    {
        this->p_chgmix->init_mixing();
        this->p_chgmix->mixing_restart_step = PARAM.inp.scf_nmax + 1;
    }

    //----------------------------------------------------------
    // for mixing restart
    // the details should move somewhere else, 2025-04-13
    //----------------------------------------------------------
    if (iter == this->p_chgmix->mixing_restart_step && PARAM.inp.mixing_restart > 0.0)
    {
        this->p_chgmix->init_mixing();
        this->p_chgmix->mixing_restart_count++;

        if (PARAM.inp.dft_plus_u)
        {
            auto* dftu = ModuleDFTU::DFTU::get_instance();
            if (dftu->uramping > 0.01 && !dftu->u_converged())
            {
                this->p_chgmix->mixing_restart_step = PARAM.inp.scf_nmax + 1;
            }
            if (dftu->uramping > 0.01)
            {
                bool do_uramping = true;
                if (PARAM.inp.sc_mag_switch)
                {
                    spinconstrain::SpinConstrain<std::complex<double>>& sc
                        = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
                    if (!sc.mag_converged()) // skip uramping if mag not converged
                    {
                        do_uramping = false;
                    }
                }
                if (do_uramping)
                {
                    dftu->uramping_update(); // update U by uramping if uramping > 0.01
                    std::cout << " U-Ramping! Current U = ";
                    for (int i = 0; i < dftu->U0.size(); i++)
                    {
                        std::cout << dftu->U[i] * ModuleBase::Ry_to_eV << " ";
                    }
                    std::cout << " eV " << std::endl;
                }
            }
        }
    }

    //----------------------------------------------------------
    // mohan move harris functional to here, 2012-06-05
    // use 'rho(in)' and 'v_h and v_xc'(in)
    //----------------------------------------------------------
    this->pelec->f_en.deband_harris = this->pelec->cal_delta_eband(ucell);

    //----------------------------------------------------------
    // update local occupations for DFT+U
    // should before lambda loop in DeltaSpin
    //----------------------------------------------------------
    if (PARAM.inp.dft_plus_u && (iter != 1 || istep != 0))
    {
        auto* dftu = ModuleDFTU::DFTU::get_instance();
        // only old DFT+U method should calculate energy correction in esolver,
        // new DFT+U method will calculate energy when evaluating the Hamiltonian
        if (dftu->omc != 2)
        {
            dftu->cal_occ_pw(iter, this->kspw_psi, this->pelec->wg, ucell, PARAM.inp.mixing_beta);
        }
        dftu->output(ucell);
    }
}

// Temporary, it should be replaced by hsolver later.
template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr)
{
    ModuleBase::timer::tick("ESolver_KS_PW", "hamilt2rho_single");

    // reset energy
    this->pelec->f_en.eband = 0.0;
    this->pelec->f_en.demet = 0.0;

    // choose if psi should be diag in subspace
    // be careful that istep start from 0 and iter start from 1
    // if (iter == 1)
    hsolver::DiagoIterAssist<T, Device>::need_subspace = ((istep == 0 || istep == 1) && iter == 1) ? false : true;

    hsolver::DiagoIterAssist<T, Device>::SCF_ITER = iter;
    hsolver::DiagoIterAssist<T, Device>::PW_DIAG_THR = ethr;

    if (PARAM.inp.calculation != "nscf")
    {
        hsolver::DiagoIterAssist<T, Device>::PW_DIAG_NMAX = PARAM.inp.pw_diag_nmax;
    }

    bool skip_charge = PARAM.inp.calculation == "nscf" ? true : false;

    // run the inner lambda loop to contrain atomic moments with the DeltaSpin method
    bool skip_solve = false;

    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<std::complex<double>>& sc
            = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
        if (!sc.mag_converged() && this->drho > 0 && this->drho < PARAM.inp.sc_scf_thr)
        {
            // optimize lambda to get target magnetic moments, but the lambda is not near target
            sc.run_lambda_loop(iter - 1);
            sc.set_mag_converged(true);
            skip_solve = true;
        }
        else if (sc.mag_converged())
        {
            // optimize lambda to get target magnetic moments, but the lambda is not near target
            sc.run_lambda_loop(iter - 1);
            skip_solve = true;
        }
    }
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
                                                     hsolver::DiagoIterAssist<T, Device>::need_subspace);

        hsolver_pw_obj.solve(this->p_hamilt,
                             this->kspw_psi[0],
                             this->pelec,
                             this->pelec->ekb.c,
                             GlobalV::RANK_IN_POOL,
                             GlobalV::NPROC_IN_POOL,
                             skip_charge,
                             ucell.tpiba,
                             ucell.nat);
    }

    // symmetrize the charge density
    Symmetry_rho srho;
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        srho.begin(is, this->chr, this->pw_rhod, ucell.symm);
    }

    ModuleBase::timer::tick("ESolver_KS_PW", "hamilt2rho_single");
}

// Temporary, it should be rewritten with Hamilt class.
template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::update_pot(UnitCell& ucell, const int istep, const int iter, const bool conv_esolver)
{
    if (!conv_esolver)
    {
        elecstate::cal_ux(ucell);
        this->pelec->pot->update_from_charge(&this->chr, &ucell);
        this->pelec->f_en.descf = this->pelec->cal_delta_escf();
#ifdef __MPI
        MPI_Bcast(&(this->pelec->f_en.descf), 1, MPI_DOUBLE, 0, BP_WORLD);
#endif
    }
    else
    {
        this->pelec->cal_converged();
    }
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver)
{
    //----------------------------------------------------------
    // Related to EXX
    //----------------------------------------------------------
    if (GlobalC::exx_info.info_global.cal_exx && !exx_helper.op_exx->first_iter)
    {
        this->pelec->set_exx(exx_helper.cal_exx_energy(kspw_psi));
    }

    //----------------------------------------------------------
    // deband is calculated from "output" charge density calculated
    // in sum_band
    // need 'rho(out)' and 'vr (v_h(in) and v_xc(in))'
    //----------------------------------------------------------
    this->pelec->f_en.deband = this->pelec->cal_delta_eband(ucell);

    //----------------------------------------------------------
    // 1) Call iter_finish() of ESolver_KS
    //----------------------------------------------------------
    ESolver_KS<T, Device>::iter_finish(ucell, istep, iter, conv_esolver);

    //----------------------------------------------------------
    // 2) Update USPP-related quantities
    // D in USPP needs vloc, thus needs update when veff updated
    // calculate the effective coefficient matrix for non-local
    // pp projectors, liuyu 2023-10-24
    //----------------------------------------------------------
    if (PARAM.globalv.use_uspp)
    {
        ModuleBase::matrix veff = this->pelec->pot->get_effective_v();
        this->ppcell.cal_effective_D(veff, this->pw_rhod, ucell);
    }

    if (GlobalC::exx_info.info_global.cal_exx)
    {
        if (GlobalC::exx_info.info_global.separate_loop)
        {
            if (conv_esolver)
            {
                auto start = std::chrono::high_resolution_clock::now();
                exx_helper.set_firstiter(false);
                exx_helper.set_psi(this->kspw_psi);

                conv_esolver = exx_helper.exx_after_converge(iter);

                if (!conv_esolver)
                {
                    auto duration = std::chrono::high_resolution_clock::now() - start;
                    std::cout << " Setting Psi for EXX PW Inner Loop took "
                              << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() / 1000.0 << "s"
                              << std::endl;
                    exx_helper.op_exx->first_iter = false;
                    XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
                    update_pot(ucell, istep, iter, conv_esolver);
                }
            }
        }
        else
        {
            exx_helper.set_psi(this->kspw_psi);
        }
    }

    //----------------------------------------------------------
    // 3) Print out electronic wavefunctions in pw basis
    //----------------------------------------------------------
	if (iter % PARAM.inp.out_freq_elec == 0 || iter == PARAM.inp.scf_nmax || conv_esolver)
	{
		ModuleIO::write_wfc_pw(GlobalV::KPAR, GlobalV::MY_POOL, GlobalV::MY_RANK, 
				PARAM.inp.nbands, PARAM.inp.nspin, PARAM.globalv.npol,
				GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL, 
				PARAM.inp.out_wfc_pw, PARAM.inp.ecutwfc, PARAM.globalv.global_out_dir,
				this->psi[0], this->kv, this->pw_wfc, GlobalV::ofs_running);
	}

    //----------------------------------------------------------
    // 4) check if oscillate for delta_spin method
    //----------------------------------------------------------
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<std::complex<double>>& sc
            = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
        if (!sc.higher_mag_prec)
        {
            sc.higher_mag_prec
                = this->p_chgmix->if_scf_oscillate(iter, this->drho, PARAM.inp.sc_os_ndim, PARAM.inp.scf_os_thr);
            if (sc.higher_mag_prec)
            { // if oscillate, increase the precision of magnetization and do mixing_restart in next iteration
                this->p_chgmix->mixing_restart_step = iter + 1;
            }
        }
    }
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS_PW", "after_scf");
    ModuleBase::timer::tick("ESolver_KS_PW", "after_scf");

    //------------------------------------------------------------------
    // 1) since ESolver_KS::psi is hidden by ESolver_KS_PW::psi,
    // we need to copy the data from ESolver_KS::psi to ESolver_KS_PW::psi.
    // This part needs to be removed when we have a better design.
    // sunliang 2025-04-10
    //------------------------------------------------------------------
    if (PARAM.inp.out_elf[0] > 0)
    {
        this->ESolver_KS<T, Device>::psi = new psi::Psi<T>(this->psi[0]);
    }

    //------------------------------------------------------------------
    // 2) call after_scf() of ESolver_KS
    //------------------------------------------------------------------
    ESolver_KS<T, Device>::after_scf(ucell, istep, conv_esolver);

    //------------------------------------------------------------------
    // 3) transfer data from GPU to CPU in pw basis
    //------------------------------------------------------------------
    if (this->device == base_device::GpuDevice)
    {
        castmem_2d_d2h_op()(this->psi[0].get_pointer() - this->psi[0].get_psi_bias(),
                            this->kspw_psi[0].get_pointer() - this->kspw_psi[0].get_psi_bias(),
                            this->psi[0].size());
    }

    //------------------------------------------------------------------
    // 4) calculate band-decomposed (partial) charge density in pw basis
    //------------------------------------------------------------------
    const std::vector<int> out_pchg = PARAM.inp.out_pchg;
    if (out_pchg.size() > 0)
    {
        ModuleIO::get_pchg_pw(out_pchg,
                              this->kspw_psi->get_nbands(),
                              PARAM.inp.nspin,
                              this->pw_rhod->nx,
                              this->pw_rhod->ny,
                              this->pw_rhod->nz,
                              this->pw_rhod->nxyz,
                              this->kv.get_nks(),
                              this->kv.isk,
                              this->kv.wk,
                              this->pw_big->bz,
                              this->pw_big->nbz,
                              this->chr.ngmc,
                              &ucell,
                              this->psi,
                              this->pw_rhod,
                              this->pw_wfc,
                              this->ctx,
                              this->Pgrid,
                              PARAM.globalv.global_out_dir,
                              PARAM.inp.if_separate_k);
    }


	// tmp 2025-05-17, mohan note
	ModuleIO::write_wfc_pw(GlobalV::KPAR, GlobalV::MY_POOL, GlobalV::MY_RANK, 
			PARAM.inp.nbands, PARAM.inp.nspin, PARAM.globalv.npol,
			GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL, 
			PARAM.inp.out_wfc_pw, PARAM.inp.ecutwfc, PARAM.globalv.global_out_dir,
			this->psi[0], this->kv, this->pw_wfc, GlobalV::ofs_running);

    //------------------------------------------------------------------
    //! 5) calculate Wannier functions in pw basis
    //------------------------------------------------------------------
    if (PARAM.inp.calculation == "nscf" && PARAM.inp.towannier90)
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Wannier functions calculation");
        toWannier90_PW wan(PARAM.inp.out_wannier_mmn,
                           PARAM.inp.out_wannier_amn,
                           PARAM.inp.out_wannier_unk,
                           PARAM.inp.out_wannier_eig,
                           PARAM.inp.out_wannier_wvfn_formatted,
                           PARAM.inp.nnkpfile,
                           PARAM.inp.wannier_spin);
        wan.set_tpiba_omega(ucell.tpiba, ucell.omega);
        wan.calculate(ucell, this->pelec->ekb, this->pw_wfc, this->pw_big, this->kv, this->psi);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Wannier functions calculation");
    }

    //------------------------------------------------------------------
    //! 6) calculate Berry phase polarization in pw basis
    //------------------------------------------------------------------
    if (PARAM.inp.calculation == "nscf" && berryphase::berry_phase_flag && ModuleSymmetry::Symmetry::symm_flag != 1)
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Berry phase polarization");
        berryphase bp;
        bp.Macroscopic_polarization(ucell, this->pw_wfc->npwk_max, this->psi, this->pw_rho, this->pw_wfc, this->kv);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Berry phase polarization");
    }

    //------------------------------------------------------------------
    // 7) write spin constrian results in pw basis
    // spin constrain calculations, write atomic magnetization and magnetic force.
    //------------------------------------------------------------------
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<std::complex<double>>& sc
            = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
        sc.cal_mi_pw();
        sc.print_Mag_Force(GlobalV::ofs_running);
    }

    //------------------------------------------------------------------
    // 8) write onsite occupations for charge and magnetizations
    //------------------------------------------------------------------
    if (PARAM.inp.onsite_radius > 0)
    { // float type has not been implemented
        auto* onsite_p = projectors::OnsiteProjector<double, Device>::get_instance();
        onsite_p->cal_occupations(reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->kspw_psi),
                                  this->pelec->wg);
    }

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

    if (this->__kspw_psi != nullptr && PARAM.inp.precision == "single")
    {
        delete reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->__kspw_psi);
    }

    // Refresh __kspw_psi
    this->__kspw_psi = PARAM.inp.precision == "single"
                           ? new psi::Psi<std::complex<double>, Device>(this->kspw_psi[0])
                           : reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->kspw_psi);

    // Calculate forces
    ff.cal_force(ucell,
                 force,
                 *this->pelec,
                 this->pw_rhod,
                 &ucell.symm,
                 &this->sf,
                 this->solvent,
                 &this->locpp,
                 &this->ppcell,
                 &this->kv,
                 this->pw_wfc,
                 this->__kspw_psi);
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::cal_stress(UnitCell& ucell, ModuleBase::matrix& stress)
{
    Stress_PW<double, Device> ss(this->pelec);

    if (this->__kspw_psi != nullptr && PARAM.inp.precision == "single")
    {
        delete reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->__kspw_psi);
    }

    // Refresh __kspw_psi
    this->__kspw_psi = PARAM.inp.precision == "single"
                           ? new psi::Psi<std::complex<double>, Device>(this->kspw_psi[0])
                           : reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->kspw_psi);
    ss.cal_stress(stress,
                  ucell,
                  this->locpp,
                  this->ppcell,
                  this->pw_rhod,
                  &ucell.symm,
                  &this->sf,
                  &this->kv,
                  this->pw_wfc,
                  this->__kspw_psi);

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
    //----------------------------------------------------------
    //! 1) after_all_runners in ESolver_KS
    //----------------------------------------------------------
    ESolver_KS<T, Device>::after_all_runners(ucell);

    //----------------------------------------------------------
    //! 2) Compute density of states (DOS)
    //----------------------------------------------------------
    if (PARAM.inp.out_dos)
    {
        ModuleIO::write_dos_pw(ucell,
                               this->pelec->ekb,
                               this->pelec->wg,
                               this->kv,
                               PARAM.inp.nbands,
                               this->pelec->eferm,
                               PARAM.inp.dos_edelta_ev,
                               PARAM.inp.dos_scale,
                               PARAM.inp.dos_sigma,
                               GlobalV::ofs_running);
    }

    //----------------------------------------------------------
    //! 3) Compute LDOS
    //----------------------------------------------------------
    if (PARAM.inp.out_ldos[0])
    {
        ModuleIO::cal_ldos_pw(reinterpret_cast<elecstate::ElecStatePW<std::complex<double>>*>(this->pelec),
                              this->psi[0],
                              this->Pgrid,
                              ucell);
    }

    //----------------------------------------------------------
    //! 4) Calculate the spillage value,
    //! which are used to generate numerical atomic orbitals
    //----------------------------------------------------------
    if (PARAM.inp.basis_type == "pw" && winput::out_spillage)
    {
        // ! Print out overlap matrices
        if (winput::out_spillage <= 2)
        {
            for (int i = 0; i < PARAM.inp.bessel_nao_rcuts.size(); i++)
            {
                if (GlobalV::MY_RANK == 0)
                {
                    std::cout << "update value: bessel_nao_rcut <- " << std::fixed << PARAM.inp.bessel_nao_rcuts[i]
                              << " a.u." << std::endl;
                }
                Numerical_Basis numerical_basis;
                numerical_basis.output_overlap(this->psi[0], this->sf, this->kv, this->pw_wfc, ucell, i);
            }
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "BASIS OVERLAP (Q and S) GENERATION.");
        }
    }

    //----------------------------------------------------------
    //! 5) Print out electronic wave functions in real space
    //----------------------------------------------------------

    //----------------------------------------------------------
    //! The write_psi_r_1 interface will be removed in the very 
    //! near future. Don't use it!
    //----------------------------------------------------------
    // if (PARAM.inp.out_wfc_r == 1) // Peize Lin add 2021.11.21
    // {
    //     ModuleIO::write_psi_r_1(ucell, this->psi[0], this->pw_wfc, "wfc_realspace", true, this->kv);
    // }

    const std::vector<int> out_wfc_norm = PARAM.inp.out_wfc_norm;
    const std::vector<int> out_wfc_re_im = PARAM.inp.out_wfc_re_im;
    if (out_wfc_norm.size() > 0 || out_wfc_re_im.size() > 0)
    {
        ModuleIO::get_wf_pw(out_wfc_norm,
                            out_wfc_re_im,
                            this->kspw_psi->get_nbands(),
                            PARAM.inp.nspin,
                            this->pw_rhod->nx,
                            this->pw_rhod->ny,
                            this->pw_rhod->nz,
                            this->pw_rhod->nxyz,
                            this->kv.get_nks(),
                            this->kv.isk,
                            this->kv.wk,
                            this->pw_big->bz,
                            this->pw_big->nbz,
                            this->chr.ngmc,
                            &ucell,
                            this->psi,
                            this->pw_rhod,
                            this->pw_wfc,
                            this->ctx,
                            this->Pgrid,
                            PARAM.globalv.global_out_dir);
    }

    //----------------------------------------------------------
    //! 6) Use Kubo-Greenwood method to compute conductivities
    //----------------------------------------------------------
    if (PARAM.inp.cal_cond)
    {
        EleCond<Real, Device> elec_cond(&ucell, &this->kv, this->pelec, this->pw_wfc, this->kspw_psi, &this->ppcell);
        elec_cond.KG(PARAM.inp.cond_smear,
                     PARAM.inp.cond_fwhm,
                     PARAM.inp.cond_wcut,
                     PARAM.inp.cond_dw,
                     PARAM.inp.cond_dt,
                     PARAM.inp.cond_nonlocal,
                     this->pelec->wg);
    }

#ifdef __MLKEDF
    //----------------------------------------------------------
    //! 7) generate training data for ML-KEDF
    //----------------------------------------------------------
    if (PARAM.inp.of_ml_gene_data == 1)
    {
        this->pelec->pot->update_from_charge(&this->chr, &ucell);

        ML_data ml_data;
        ml_data.set_para(this->chr.nrxx,
                         PARAM.inp.nelec,
                         PARAM.inp.of_tf_weight,
                         PARAM.inp.of_vw_weight,
                         PARAM.inp.of_ml_chi_p,
                         PARAM.inp.of_ml_chi_q,
                         PARAM.inp.of_ml_chi_xi,
                         PARAM.inp.of_ml_chi_pnl,
                         PARAM.inp.of_ml_chi_qnl,
                         PARAM.inp.of_ml_nkernel,
                         PARAM.inp.of_ml_kernel,
                         PARAM.inp.of_ml_kernel_scaling,
                         PARAM.inp.of_ml_yukawa_alpha,
                         PARAM.inp.of_ml_kernel_file,
                         ucell.omega,
                         this->pw_rho);

        ml_data.generateTrainData_KS(this->kspw_psi,
                                     this->pelec,
                                     this->pw_wfc,
                                     this->pw_rho,
                                     ucell,
                                     this->pelec->pot->get_effective_v(0));
    }
#endif
}

template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace ModuleESolver
