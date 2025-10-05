#include "esolver_ks_pw.h"

#include "source_estate/cal_ux.h"
#include "source_estate/elecstate_pw.h"
#include "source_estate/module_charge/symmetry_rho.h"

#include "source_hsolver/diago_iter_assist.h"
#include "source_hsolver/hsolver_pw.h"

#include "source_hsolver/kernels/hegvd_op.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_pw/module_pwdft/onsite_projector.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_pw/module_pwdft/VSep_in_pw.h"
#include "source_pw/module_pwdft/hamilt_pw.h"

#include "source_pw/module_pwdft/forces.h"
#include "source_pw/module_pwdft/stress_pw.h"

#ifdef __DSP
#include "source_base/kernels/dsp/dsp_connector.h"
#endif

#include <chrono>

#include "source_pw/module_pwdft/setup_pot.h" // mohan add 20250929
#include "source_estate/setup_estate_pw.h" // mohan add 20251005
#include "source_io/ctrl_output_pw.h" // mohan add 20250927

namespace ModuleESolver
{

template <typename T, typename Device>
ESolver_KS_PW<T, Device>::ESolver_KS_PW()
{
    this->classname = "ESolver_KS_PW";
    this->basisname = "PW";
    this->device = base_device::get_device_type<Device>(this->ctx);
}

template <typename T, typename Device>
ESolver_KS_PW<T, Device>::~ESolver_KS_PW()
{
	//****************************************************
	// do not add any codes in this deconstructor funcion
	//****************************************************
	// delete Hamilt
    this->deallocate_hamilt();

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
    //! Call before_all_runners() of ESolver_KS
    ESolver_KS<T, Device>::before_all_runners(ucell, inp);

    //! setup and allocation for pelec, charge density, potentials, etc. 
    elecstate::setup_estate_pw<T, Device>(ucell, this->kv, this->sf, this->pelec, this->chr,
      this->locpp, this->ppcell, this->vsep_cell, this->pw_wfc, this->pw_rho,
      this->pw_rhod, this->pw_big, this->solvent, inp);

    //! Allocate and initialize psi
    this->p_psi_init = new psi::PSIInit<T, Device>(inp.init_wfc,
      inp.ks_solver, inp.basis_type, GlobalV::MY_RANK, ucell,
      this->sf, this->kv, this->ppcell, *this->pw_wfc);

    allocate_psi(this->psi, this->kv.get_nks(), this->kv.ngk, PARAM.globalv.nbands_l, this->pw_wfc->npwk_max);

    this->p_psi_init->prepare_init(inp.pw_seed);

    this->kspw_psi = inp.device == "gpu" || inp.precision == "single"
                         ? new psi::Psi<T, Device>(this->psi[0])
                         : reinterpret_cast<psi::Psi<T, Device>*>(this->psi);

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT BASIS");

    //! Initialize exx pw
    if (inp.calculation == "scf" || inp.calculation == "relax" || inp.calculation == "cell-relax"
        || inp.calculation == "md")
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

    //! Call before_scf() of ESolver_KS
    ESolver_KS<T, Device>::before_scf(ucell, istep);

    //! Init variables (once the cell has changed)
    if (ucell.cell_parameter_updated)
    {
        this->ppcell.rescale_vnl(ucell.omega);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

        this->pw_wfc->initgrids(ucell.lat0, ucell.latvec, this->pw_wfc->nx, this->pw_wfc->ny, this->pw_wfc->nz);

        this->pw_wfc->initparameters(false, PARAM.inp.ecutwfc, this->kv.get_nks(), this->kv.kvec_d.data());

        this->pw_wfc->collect_local_pw(PARAM.inp.erf_ecut, PARAM.inp.erf_height, PARAM.inp.erf_sigma);

        this->p_psi_init->prepare_init(PARAM.inp.pw_seed);
    }

    //! Init Hamiltonian (cell changed)
    //! Operators in HamiltPW should be reallocated once cell changed
    //! delete Hamilt if not first scf
    this->deallocate_hamilt();

    //! Allocate HamiltPW
    this->allocate_hamilt(ucell);

    //! Setup potentials (local, non-local, sc, +U, DFT-1/2)
    pw::setup_pot(istep, ucell, this->kv, this->sf, this->pelec, this->Pgrid,
              this->chr, this->locpp, this->ppcell, this->vsep_cell,
              this->kspw_psi, this->p_hamilt, this->pw_wfc, this->pw_rhod, PARAM.inp);

    //! Initialize wave functions
    if (!this->already_initpsi)
    {
        this->p_psi_init->initialize_psi(this->psi, this->kspw_psi, this->p_hamilt, GlobalV::ofs_running);
        this->already_initpsi = true;
    }

    //! Exx calculations
    if (PARAM.inp.calculation == "scf" || PARAM.inp.calculation == "relax" 
        || PARAM.inp.calculation == "cell-relax" || PARAM.inp.calculation == "md")
    {
        if (GlobalC::exx_info.info_global.cal_exx && PARAM.inp.basis_type == "pw")
        {
            auto hamilt_pw = reinterpret_cast<hamilt::HamiltPW<T, Device>*>(this->p_hamilt);
            hamilt_pw->set_exx_helper(exx_helper);
            exx_helper.set_psi(kspw_psi);
        }
    }

    ModuleBase::timer::tick("ESolver_KS_PW", "before_scf");
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::iter_init(UnitCell& ucell, const int istep, const int iter)
{
    // Call iter_init() of ESolver_KS
    ESolver_KS<T, Device>::iter_init(ucell, istep, iter);

    if (iter == 1)
    {
        this->p_chgmix->init_mixing();
        this->p_chgmix->mixing_restart_step = PARAM.inp.scf_nmax + 1;
    }

    // For mixing restart
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

    // mohan move harris functional to here, 2012-06-05
    // use 'rho(in)' and 'v_h and v_xc'(in)
    this->pelec->f_en.deband_harris = this->pelec->cal_delta_eband(ucell);

    // update local occupations for DFT+U
    // should before lambda loop in DeltaSpin
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
                                                     hsolver::DiagoIterAssist<T, Device>::need_subspace,
                                                     PARAM.inp.use_k_continuity);

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
    // Related to EXX
    if (GlobalC::exx_info.info_global.cal_exx && !exx_helper.op_exx->first_iter)
    {
        this->pelec->set_exx(exx_helper.cal_exx_energy(kspw_psi));
    }

    // deband is calculated from "output" charge density
    this->pelec->f_en.deband = this->pelec->cal_delta_eband(ucell);

    // Call iter_finish() of ESolver_KS
    ESolver_KS<T, Device>::iter_finish(ucell, istep, iter, conv_esolver);

    // D in USPP needs vloc, thus needs update when veff updated
    // calculate the effective coefficient matrix for non-local
    // pp projectors, liuyu 2023-10-24
    if (PARAM.globalv.use_uspp)
    {
        ModuleBase::matrix veff = this->pelec->pot->get_effective_v();
        this->ppcell.cal_effective_D(veff, this->pw_rhod, ucell);
    }

    // Related to EXX
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        if (GlobalC::exx_info.info_global.separate_loop)
        {
            if (conv_esolver)
            {
                auto start = std::chrono::high_resolution_clock::now();
                exx_helper.set_firstiter(false);
                exx_helper.op_exx->first_iter = false;
                double dexx = 0.0;
                if (PARAM.inp.exx_thr_type == "energy")
                {
                    dexx = exx_helper.cal_exx_energy(this->kspw_psi);
                }
                exx_helper.set_psi(this->kspw_psi);
                if (PARAM.inp.exx_thr_type == "energy")
                {
                    dexx -= exx_helper.cal_exx_energy(this->kspw_psi);
                    // std::cout << "dexx = " << dexx << std::endl;
                }
                bool conv_ene = std::abs(dexx) < PARAM.inp.exx_ene_thr;

                conv_esolver = exx_helper.exx_after_converge(iter, conv_ene);
                if (!conv_esolver)
                {
                    auto duration = std::chrono::high_resolution_clock::now() - start;
                    std::cout << " Setting Psi for EXX PW Inner Loop took "
                              << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() / 1000.0 << "s"
                              << std::endl;
                    exx_helper.op_exx->first_iter = false;
                    XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
                    update_pot(ucell, istep, iter, conv_esolver);
                    exx_helper.iter_inc();
                }
            }
        }
        else
        {
            exx_helper.set_psi(this->kspw_psi);
        }
    }

    // check if oscillate for delta_spin method
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<std::complex<double>>& sc
            = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
        if (!sc.higher_mag_prec)
        {
            sc.higher_mag_prec = this->p_chgmix->if_scf_oscillate(iter, 
              this->drho, PARAM.inp.sc_os_ndim, PARAM.inp.scf_os_thr);
            if (sc.higher_mag_prec)
            { // if oscillate, increase the precision of magnetization and do mixing_restart in next iteration
                this->p_chgmix->mixing_restart_step = iter + 1;
            }
        }
    }

    // the output quantities
    ModuleIO::ctrl_iter_pw(istep, iter, conv_esolver, this->psi, 
              this->kv, this->pw_wfc, PARAM.inp);
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS_PW", "after_scf");
    ModuleBase::timer::tick("ESolver_KS_PW", "after_scf");

    // Since ESolver_KS::psi is hidden by ESolver_KS_PW::psi,
    // we need to copy the data from ESolver_KS::psi to ESolver_KS_PW::psi.
    // sunliang 2025-04-10
    if (PARAM.inp.out_elf[0] > 0)
    {
        this->ESolver_KS<T, Device>::psi = new psi::Psi<T>(this->psi[0]);
    }

    // Call 'after_scf' of ESolver_KS
    ESolver_KS<T, Device>::after_scf(ucell, istep, conv_esolver);

    // Transfer data from GPU to CPU in pw basis
    if (this->device == base_device::GpuDevice)
    {
        castmem_2d_d2h_op()(this->psi[0].get_pointer() - this->psi[0].get_psi_bias(),
                            this->kspw_psi[0].get_pointer() - this->kspw_psi[0].get_psi_bias(),
                            this->psi[0].size());
    }

    // Output quantities
    ModuleIO::ctrl_scf_pw<T, Device>(istep, ucell, this->pelec, this->chr, this->kv, this->pw_wfc,
              this->pw_rho, this->pw_rhod, this->pw_big, this->psi, this->kspw_psi, 
              this->__kspw_psi, this->ctx, this->Pgrid, PARAM.inp);

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
    ff.cal_force(ucell, force, *this->pelec, this->pw_rhod, &ucell.symm,
                 &this->sf, this->solvent, &this->locpp, &this->ppcell, 
                 &this->kv, this->pw_wfc, this->__kspw_psi);
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

    ss.cal_stress(stress, ucell, this->locpp, this->ppcell, this->pw_rhod,
                  &ucell.symm, &this->sf, &this->kv, this->pw_wfc, this->__kspw_psi);

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
    ESolver_KS<T, Device>::after_all_runners(ucell);

    ModuleIO::ctrl_runner_pw<T, Device>(ucell, this->pelec, this->pw_wfc, 
            this->pw_rho, this->pw_rhod, this->chr, this->kv, this->psi,
            this->kspw_psi, this->__kspw_psi, this->sf, 
            this->ppcell, this->solvent, this->ctx, this->Pgrid, PARAM.inp); 

    elecstate::teardown_estate_pw<T, Device>(this->pelec, this->vsep_cell);
    
}

template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace ModuleESolver
