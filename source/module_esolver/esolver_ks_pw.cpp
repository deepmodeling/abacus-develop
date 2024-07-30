#include "esolver_ks_pw.h"

#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/elecond.h"
#include "module_io/input_conv.h"
#include "module_io/nscf_band.h"
#include "module_io/output_log.h"
#include "module_io/write_dos_pw.h"
#include "module_io/write_istate_info.h"
#include "module_io/write_wfc_pw.h"

#include <iostream>

//--------------temporary----------------------------
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"
//-----force-------------------
#include "module_hamilt_pw/hamilt_pwdft/forces.h"
//-----stress------------------
#include "module_hamilt_pw/hamilt_pwdft/stress_pw.h"
//---------------------------------------------------
#include "module_base/memory.h"
#include "module_base/module_device/device.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_hamilt_general/module_vdw/vdw.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/hsolver_pw.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_io/berryphase.h"
#include "module_io/cube_io.h"
#include "module_io/numerical_basis.h"
#include "module_io/numerical_descriptor.h"
#include "module_io/to_wannier90_pw.h"
#include "module_io/winput.h"
#include "module_io/write_pot.h"
#include "module_io/write_wfc_r.h"
#include "module_parameter/parameter.h"
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif
#include "module_base/formatter.h"

#include <ATen/kernels/blas.h>
#include <ATen/kernels/lapack.h>

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
        hsolver::createGpuBlasHandle();
        hsolver::createGpuSolverHandle();
        container::kernels::createGpuBlasHandle();
        container::kernels::createGpuSolverHandle();
    }
#endif
}

template <typename T, typename Device>
ESolver_KS_PW<T, Device>::~ESolver_KS_PW()
{
    // delete HSolver and ElecState
    this->deallocate_hsolver();
    if (this->pelec != nullptr)
    {
        delete reinterpret_cast<elecstate::ElecStatePW<T, Device>*>(this->pelec);
        this->pelec = nullptr;
    }
    // delete Hamilt
    this->deallocate_hamilt();
    if (this->device == base_device::GpuDevice)
    {
#if defined(__CUDA) || defined(__ROCM)
        hsolver::destoryBLAShandle();
        hsolver::destroyGpuSolverHandle();
        container::kernels::destroyGpuBlasHandle();
        container::kernels::destroyGpuSolverHandle();
#endif
        delete reinterpret_cast<psi::Psi<T, Device>*>(this->kspw_psi);
    }
    if (GlobalV::precision_flag == "single")
    {
        delete reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->__kspw_psi);
    }
    delete this->psi;
    delete this->p_wf_init;
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::before_all_runners(const Input_para& inp, UnitCell& ucell) {
    // 1) call before_all_runners() of ESolver_KS
    ESolver_KS<T, Device>::before_all_runners(inp, ucell);

    // 2) initialize HSolver
    if (this->phsol == nullptr)
    {
        this->allocate_hsolver();
    }

    // 3) initialize ElecState,
    if (this->pelec == nullptr)
    {
        this->pelec = new elecstate::ElecStatePW<T, Device>(this->pw_wfc,
                                                            &(this->chr),
                                                            &(this->kv),
                                                            &ucell,
                                                            &GlobalC::ppcell,
                                                            this->pw_rhod,
                                                            this->pw_rho,
                                                            this->pw_big);
    }

    //! 4) inititlize the charge density.
    this->pelec->charge->allocate(GlobalV::NSPIN);

    //! 5) set the cell volume variable in pelec
    this->pelec->omega = ucell.omega;

    //! 6) initialize the potential.
    if (this->pelec->pot == nullptr)
    {
        this->pelec->pot = new elecstate::Potential(this->pw_rhod,
                                                    this->pw_rho,
                                                    &ucell,
                                                    &GlobalC::ppcell.vloc,
                                                    &(this->sf),
                                                    &(this->pelec->f_en.etxc),
                                                    &(this->pelec->f_en.vtxc));
    }

    //! 7) prepare some parameters for electronic wave functions initilization
    this->p_wf_init = new psi::WFInit<T, Device>(GlobalV::init_wfc,
                                                 GlobalV::KS_SOLVER,
                                                 GlobalV::BASIS_TYPE,
                                                 GlobalV::psi_initializer,
                                                 &this->wf,
                                                 this->pw_wfc);
    this->p_wf_init->prepare_init(&(this->sf),
                                  &ucell,
                                  1,
#ifdef __MPI
                                  &GlobalC::Pkpoints,
                                  GlobalV::MY_RANK,
#endif
                                  &GlobalC::ppcell);

    //! 8) setup global classes
    this->Init_GlobalC(inp, ucell, GlobalC::ppcell);

    //! 9) setup occupations
    if (PARAM.inp.ocp)
    {
        this->pelec->fixed_weights(PARAM.inp.ocp_kb, GlobalV::NBANDS, GlobalV::nelec);
    }
}


template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::before_scf(const int istep)
{
    ModuleBase::TITLE("ESolver_KS_PW", "before_scf");

    if (GlobalC::ucell.cell_parameter_updated)
    {
        this->init_after_vc(PARAM.inp, GlobalC::ucell);
    }
    if (GlobalC::ucell.ionic_position_updated)
    {
        this->CE.update_all_dis(GlobalC::ucell);
        this->CE.extrapolate_charge(
#ifdef __MPI
            &(GlobalC::Pgrid),
#endif
            GlobalC::ucell,
            this->pelec->charge,
            &this->sf);
    }

    // init Hamilt, this should be allocated before each scf loop
    // Operators in HamiltPW should be reallocated once cell changed
    // delete Hamilt if not first scf
    this->deallocate_hamilt();

    // allocate HamiltPW
    this->allocate_hamilt();

    //----------------------------------------------------------
    // about vdw, jiyy add vdwd3 and linpz add vdwd2
    //----------------------------------------------------------
    auto vdw_solver = vdw::make_vdw(GlobalC::ucell, PARAM.inp);
    if (vdw_solver != nullptr)
    {
        this->pelec->f_en.evdw = vdw_solver->get_energy();
    }

    // calculate ewald energy
    if (!PARAM.inp.test_skip_ewald)
    {
        this->pelec->f_en.ewald_energy = H_Ewald_pw::compute_ewald(GlobalC::ucell, this->pw_rhod, this->sf.strucFac);
    }

    //! cal_ux should be called before init_scf because
    //! the direction of ux is used in noncoline_rho
    if (GlobalV::NSPIN == 4 && GlobalV::DOMAG)
    {
        GlobalC::ucell.cal_ux();
    }

    //! calculate the total local pseudopotential in real space
    this->pelec->init_scf(istep, this->sf.strucFac);

    if (PARAM.inp.out_chg == 2)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            std::stringstream ss;
            ss << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG_INI.cube";
            ModuleIO::write_cube(
#ifdef __MPI
                this->pw_big->bz,
                this->pw_big->nbz,
                this->pw_rho->nplane,
                this->pw_rho->startz_current,
#endif
                this->pelec->charge->rho[is],
                is,
                GlobalV::NSPIN,
                0,
                ss.str(),
                this->pw_rho->nx,
                this->pw_rho->ny,
                this->pw_rho->nz,
                this->pelec->eferm.ef,
                &(GlobalC::ucell));
        }
    }

    ModuleIO::write_pot(GlobalV::out_pot,
                        GlobalV::NSPIN,
                        GlobalV::global_out_dir,
#ifdef __MPI
                        this->pw_big->bz,
                        this->pw_big->nbz,
                        this->pw_rho->nplane,
                        this->pw_rho->startz_current,
#endif
                        this->pw_rho->nx,
                        this->pw_rho->ny,
                        this->pw_rho->nz,
                        this->pelec->pot->get_effective_v());

    //! Symmetry_rho should behind init_scf, because charge should be
    //! initialized first. liuyu comment: Symmetry_rho should be located between
    //! init_rho and v_of_rho?
    Symmetry_rho srho;
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        srho.begin(is, *(this->pelec->charge), this->pw_rhod, GlobalC::Pgrid, GlobalC::ucell.symm);
    }

    // liuyu move here 2023-10-09
    // D in uspp need vloc, thus behind init_scf()
    // calculate the effective coefficient matrix for non-local pseudopotential
    // projectors
    ModuleBase::matrix veff = this->pelec->pot->get_effective_v();

    GlobalC::ppcell.cal_effective_D(veff, this->pw_rhod, GlobalC::ucell);

    // after init_rho (in pelec->init_scf), we have rho now.
    // before hamilt2density, we update Hk and initialize psi

    // before_scf function will be called everytime before scf. However, once
    // atomic coordinates changed, structure factor will change, therefore all
    // atomwise properties will change. So we need to reinitialize psi every
    // time before scf. But for random wavefunction, we dont, because random
    // wavefunction is not related to atomic coordinates. What the old strategy
    // does is only to initialize for once...
    if (((GlobalV::init_wfc == "random") && (istep == 0)) || (GlobalV::init_wfc != "random"))
    {
        this->p_wf_init->initialize_psi(this->psi, this->kspw_psi, this->p_hamilt, GlobalV::ofs_running);
    }
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::iter_init(const int istep, const int iter)
{
    if (iter == 1)
    {
        this->p_chgmix->init_mixing();
        this->p_chgmix->mixing_restart_step = GlobalV::SCF_NMAX + 1;
    }
    // for mixing restart
    if (iter == this->p_chgmix->mixing_restart_step && GlobalV::MIXING_RESTART > 0.0)
    {
        this->p_chgmix->init_mixing();
    }
    // mohan move harris functional to here, 2012-06-05
    // use 'rho(in)' and 'v_h and v_xc'(in)
    this->pelec->f_en.deband_harris = this->pelec->cal_delta_eband();

    //(2) save change density as previous charge,
    // prepared fox mixing.
    if (GlobalV::MY_STOGROUP == 0)
    {
        this->pelec->charge->save_rho_before_sum_band();
    }
}

// Temporary, it should be replaced by hsolver later.
template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::hamilt2density(const int istep, const int iter, const double ethr)
{
    ModuleBase::timer::tick("ESolver_KS_PW", "hamilt2density");

    if (this->phsol != nullptr)
    {
        // reset energy
        this->pelec->f_en.eband = 0.0;
        this->pelec->f_en.demet = 0.0;
        // choose if psi should be diag in subspace
        // be careful that istep start from 0 and iter start from 1
        // if (iter == 1)
        hsolver::DiagoIterAssist<T, Device>::need_subspace = ((istep == 0 || istep == 1) && iter == 1) ? false : true;
        hsolver::DiagoIterAssist<T, Device>::SCF_ITER = iter;
        hsolver::DiagoIterAssist<T, Device>::PW_DIAG_THR = ethr;
        hsolver::DiagoIterAssist<T, Device>::PW_DIAG_NMAX = GlobalV::PW_DIAG_NMAX;

        this->phsol->solve(this->p_hamilt,    // hamilt::Hamilt<T, Device>* pHamilt,
                           this->kspw_psi[0], // psi::Psi<T, Device>& psi,
                           this->pelec,       // elecstate::ElecState<T, Device>* pelec,
                           PARAM.inp.ks_solver,
                           PARAM.inp.calculation,
                           PARAM.inp.basis_type,
                           PARAM.inp.use_paw,
                           GlobalV::use_uspp,
                           GlobalV::RANK_IN_POOL,
                           GlobalV::NPROC_IN_POOL,

                           hsolver::DiagoIterAssist<T, Device>::SCF_ITER,
                           hsolver::DiagoIterAssist<T, Device>::need_subspace,
                           hsolver::DiagoIterAssist<T, Device>::PW_DIAG_NMAX,
                           hsolver::DiagoIterAssist<T, Device>::PW_DIAG_THR,

                           false);

        if (PARAM.inp.out_bandgap)
        {
            if (!GlobalV::TWO_EFERMI)
            {
                this->pelec->cal_bandgap();
            }
            else
            {
                this->pelec->cal_bandgap_updw();
            }
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW", "HSolver has not been initialed!");
    }

    // calculate the delta_harris energy
    // according to new charge density.
    // mohan add 2009-01-23
    this->pelec->cal_energies(1);

    Symmetry_rho srho;
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        srho.begin(is, *(this->pelec->charge), this->pw_rhod, GlobalC::Pgrid, GlobalC::ucell.symm);
    }

    // compute magnetization, only for LSDA(spin==2)
    GlobalC::ucell.magnet.compute_magnetization(this->pelec->charge->nrxx,
                                                this->pelec->charge->nxyz,
                                                this->pelec->charge->rho,
                                                this->pelec->nelec_spin.data());

    // deband is calculated from "output" charge density calculated
    // in sum_band
    // need 'rho(out)' and 'vr (v_h(in) and v_xc(in))'
    this->pelec->f_en.deband = this->pelec->cal_delta_eband();

    ModuleBase::timer::tick("ESolver_KS_PW", "hamilt2density");
}

// Temporary, it should be rewritten with Hamilt class.
template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::update_pot(const int istep, const int iter)
{
    if (!this->conv_elec)
    {
        if (GlobalV::NSPIN == 4)
        {
            GlobalC::ucell.cal_ux();
        }
        this->pelec->pot->update_from_charge(this->pelec->charge, &GlobalC::ucell);
        this->pelec->f_en.descf = this->pelec->cal_delta_escf();
    }
    else
    {
        this->pelec->cal_converged();
    }
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::iter_finish(const int iter)
{
    // liuyu 2023-10-24
    // D in uspp need vloc, thus needs update when veff updated
    // calculate the effective coefficient matrix for non-local pseudopotential
    // projectors
    if (GlobalV::use_uspp)
    {
        ModuleBase::matrix veff = this->pelec->pot->get_effective_v();
        GlobalC::ppcell.cal_effective_D(veff, this->pw_rhod, GlobalC::ucell);
    }

    // 1 means Harris-Foulkes functional
    // 2 means Kohn-Sham functional
    const int energy_type = 2;
    this->pelec->cal_energies(2);

    bool print = false;
    if (this->out_freq_elec && iter % this->out_freq_elec == 0)
    {
        print = true;
    }

    if (print == true)
    {
        if (PARAM.inp.out_chg > 0)
        {
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                double* data = nullptr;
                if (PARAM.inp.dm_to_rho)
                {
                    data = this->pelec->charge->rho[is];
                }
                else
                {
                    data = this->pelec->charge->rho_save[is];
                }
                std::string fn = GlobalV::global_out_dir + "/tmp_SPIN" + std::to_string(is + 1) + "_CHG.cube";
                ModuleIO::write_cube(
#ifdef __MPI
                    this->pw_big->bz,
                    this->pw_big->nbz,
                    this->pw_rhod->nplane,
                    this->pw_rhod->startz_current,
#endif
                    data,
                    is,
                    GlobalV::NSPIN,
                    iter,
                    fn,
                    this->pw_rhod->nx,
                    this->pw_rhod->ny,
                    this->pw_rhod->nz,
                    this->pelec->eferm.get_efval(is),
                    &(GlobalC::ucell),
                    3,
                    1);
                if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
                {
                    fn = GlobalV::global_out_dir + "/tmp_SPIN" + std::to_string(is + 1) + "_TAU.cube";
                    ModuleIO::write_cube(
#ifdef __MPI
                        this->pw_big->bz,
                        this->pw_big->nbz,
                        this->pw_rhod->nplane,
                        this->pw_rhod->startz_current,
#endif
                        this->pelec->charge->kin_r_save[is],
                        is,
                        GlobalV::NSPIN,
                        iter,
                        fn,
                        this->pw_rhod->nx,
                        this->pw_rhod->ny,
                        this->pw_rhod->nz,
                        this->pelec->eferm.get_efval(is),
                        &(GlobalC::ucell));
                }
            }
        }
        // output wavefunctions
        if (this->wf.out_wfc_pw == 1 || this->wf.out_wfc_pw == 2)
        {
            std::stringstream ssw;
            ssw << GlobalV::global_out_dir << "WAVEFUNC";
            // mohan update 2011-02-21
            // qianrui update 2020-10-17
            ModuleIO::write_wfc_pw(ssw.str(), this->psi[0], this->kv, this->pw_wfc);
            // ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"write wave
            // functions into file WAVEFUNC.dat");
        }
    }
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::after_scf(const int istep)
{
    // 1) call after_scf() of ESolver_FP
    ESolver_FP::after_scf(istep);

    this->create_Output_Potential(istep).write();

    if (PARAM.inp.out_chg)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            double* data = nullptr;
            if (PARAM.inp.dm_to_rho)
            {
                data = this->pelec->charge->rho[is];
            }
            else
            {
                data = this->pelec->charge->rho_save[is];
            }
            std::string fn = GlobalV::global_out_dir + "/SPIN" + std::to_string(is + 1) + "_CHG.cube";
            ModuleIO::write_cube(
#ifdef __MPI
                this->pw_big->bz,
                this->pw_big->nbz,
                this->pw_rhod->nplane,
                this->pw_rhod->startz_current,
#endif
                data,
                is,
                GlobalV::NSPIN,
                istep,
                fn,
                this->pw_rhod->nx,
                this->pw_rhod->ny,
                this->pw_rhod->nz,
                this->pelec->eferm.get_efval(is),
                &(GlobalC::ucell),
                3,
                1);
            if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
            {
                fn = GlobalV::global_out_dir + "/SPIN" + std::to_string(is + 1) + "_TAU.cube";
                ModuleIO::write_cube(
#ifdef __MPI
                    this->pw_big->bz,
                    this->pw_big->nbz,
                    this->pw_rhod->nplane,
                    this->pw_rhod->startz_current,
#endif
                    this->pelec->charge->kin_r_save[is],
                    is,
                    GlobalV::NSPIN,
                    istep,
                    fn,
                    this->pw_rhod->nx,
                    this->pw_rhod->ny,
                    this->pw_rhod->nz,
                    this->pelec->eferm.get_efval(is),
                    &(GlobalC::ucell));
            }
        }
    }

    if (this->wf.out_wfc_pw == 1 || this->wf.out_wfc_pw == 2)
    {
        std::stringstream ssw;
        ssw << GlobalV::global_out_dir << "WAVEFUNC";
        ModuleIO::write_wfc_pw(ssw.str(), this->psi[0], this->kv, this->pw_wfc);
    }

    ModuleIO::output_convergence_after_scf(this->conv_elec, this->pelec->f_en.etot);

    ModuleIO::output_efermi(this->conv_elec, this->pelec->eferm.ef);

    if (GlobalV::OUT_LEVEL != "m")
    {
        this->pelec->print_eigenvalue(GlobalV::ofs_running);
    }

    if (this->device == base_device::GpuDevice)
    {
        castmem_2d_d2h_op()(this->psi[0].get_device(),
                            this->kspw_psi[0].get_device(),
                            this->psi[0].get_pointer() - this->psi[0].get_psi_bias(),
                            this->kspw_psi[0].get_pointer() - this->kspw_psi[0].get_psi_bias(),
                            this->psi[0].size());
    }

    // Get bands_to_print through public function of INPUT (returns a const
    // pointer to string)
    const std::vector<int> bands_to_print = PARAM.inp.bands_to_print;
    if (bands_to_print.size() > 0)
    {
        // bands_picked is a vector of 0s and 1s, where 1 means the band is
        // picked to output
        std::vector<int> bands_picked;
        bands_picked.resize(this->kspw_psi->get_nbands());
        ModuleBase::GlobalFunc::ZEROS(bands_picked.data(), this->kspw_psi->get_nbands());

        // Check if length of bands_to_print is valid
        if (static_cast<int>(bands_to_print.size()) > this->kspw_psi->get_nbands())
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_PW::after_scf",
                                     "The number of bands specified by `bands_to_print` in the "
                                     "INPUT file exceeds `nbands`!");
        }

        // Check if all elements in bands_picked are 0 or 1
        for (int value: bands_to_print)
        {
            if (value != 0 && value != 1)
            {
                ModuleBase::WARNING_QUIT("ESolver_KS_PW::after_scf",
                                         "The elements of `bands_to_print` must be either 0 or 1. "
                                         "Invalid values found!");
            }
        }

        // Fill bands_picked with values from bands_to_print
        // Remaining bands are already set to 0
        int length = std::min(static_cast<int>(bands_to_print.size()), this->kspw_psi->get_nbands());
        for (int i = 0; i < length; ++i)
        {
            // bands_to_print rely on function parse_expression
            // Initially designed for ocp_set, which can be double
            bands_picked[i] = static_cast<int>(bands_to_print[i]);
        }

        std::complex<double>* wfcr = new std::complex<double>[this->pw_rho->nxyz];
        double* rho_band = new double[this->pw_rho->nxyz];

        for (int ib = 0; ib < this->kspw_psi->get_nbands(); ++ib)
        {
            // Skip the loop iteration if bands_picked[ib] is 0
            if (!bands_picked[ib])
            {
                continue;
            }

            for (int i = 0; i < this->pw_rho->nxyz; i++)
            {
                // Initialize rho_band to zero for each band
                rho_band[i] = 0.0;
            }

            for (int ik = 0; ik < this->kv.get_nks(); ik++)
            {
                this->psi->fix_k(ik);
                this->pw_wfc->recip_to_real(this->ctx, &psi[0](ib, 0), wfcr, ik);

                double w1 = static_cast<double>(this->kv.wk[ik] / GlobalC::ucell.omega);

                for (int i = 0; i < this->pw_rho->nxyz; i++)
                {
                    rho_band[i] += std::norm(wfcr[i]) * w1;
                }
            }

            // Symmetrize the charge density, otherwise the results are incorrect if the symmetry is on
            std::cout << " Symmetrizing band-decomposed charge density..." << std::endl;
            Symmetry_rho srho;
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                // Use vector instead of raw pointers
                std::vector<double*> rho_save_pointers(GlobalV::NSPIN, rho_band);
                std::vector<std::vector<std::complex<double>>> rhog(
                    GlobalV::NSPIN,
                    std::vector<std::complex<double>>(this->pelec->charge->ngmc));

                // Convert vector of vectors to vector of pointers
                std::vector<std::complex<double>*> rhog_pointers(GlobalV::NSPIN);
                for (int s = 0; s < GlobalV::NSPIN; s++)
                {
                    rhog_pointers[s] = rhog[s].data();
                }

                srho.begin(is,
                           rho_save_pointers.data(),
                           rhog_pointers.data(),
                           this->pelec->charge->ngmc,
                           nullptr,
                           this->pw_rhod,
                           GlobalC::Pgrid,
                           GlobalC::ucell.symm);
            }

            std::stringstream ssc;
            ssc << GlobalV::global_out_dir << "BAND" << ib + 1 << "_CHG.cube"; // band index starts from 1

            ModuleIO::write_cube(
#ifdef __MPI
                this->pw_big->bz,
                this->pw_big->nbz,
                this->pw_big->nplane,
                this->pw_big->startz_current,
#endif
                rho_band,
                0,
                GlobalV::NSPIN,
                0,
                ssc.str(),
                this->pw_rho->nx,
                this->pw_rho->ny,
                this->pw_rho->nz,
                0.0,
                &(GlobalC::ucell));
        }
        delete[] wfcr;
        delete[] rho_band;
    }
}

template <typename T, typename Device>
double ESolver_KS_PW<T, Device>::cal_energy()
{
    return this->pelec->f_en.etot;
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::cal_force(ModuleBase::matrix& force)
{
    Forces<double, Device> ff(GlobalC::ucell.nat);
    if (this->__kspw_psi != nullptr && GlobalV::precision_flag == "single")
    {
        delete reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->__kspw_psi);
    }

    // Refresh __kspw_psi
    this->__kspw_psi = GlobalV::precision_flag == "single"
                           ? new psi::Psi<std::complex<double>, Device>(this->kspw_psi[0])
                           : reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->kspw_psi);

    // Calculate forces
    ff.cal_force(force,
                 *this->pelec,
                 this->pw_rhod,
                 &GlobalC::ucell.symm,
                 &this->sf,
                 &this->kv,
                 this->pw_wfc,
                 this->__kspw_psi);
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::cal_stress(ModuleBase::matrix& stress)
{
    Stress_PW<double, Device> ss(this->pelec);
    if (this->__kspw_psi != nullptr && GlobalV::precision_flag == "single")
    {
        delete reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->__kspw_psi);
    }

    // Refresh __kspw_psi
    this->__kspw_psi = GlobalV::precision_flag == "single"
                           ? new psi::Psi<std::complex<double>, Device>(this->kspw_psi[0])
                           : reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->kspw_psi);
    ss.cal_stress(stress,
                  GlobalC::ucell,
                  &GlobalC::ppcell,
                  this->pw_rhod,
                  &GlobalC::ucell.symm,
                  &this->sf,
                  &this->kv,
                  this->pw_wfc,
                  this->psi,
                  this->__kspw_psi);

    // external stress
    double unit_transform = 0.0;
    unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    double external_stress[3] = {PARAM.inp.press1, PARAM.inp.press2, PARAM.inp.press3};
    for (int i = 0; i < 3; i++)
    {
        stress(i, i) -= external_stress[i] / unit_transform;
    }
    GlobalV::PRESSURE = (stress(0, 0) + stress(1, 1) + stress(2, 2)) / 3;
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::after_all_runners()
{
    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << this->pelec->f_en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;

    if (PARAM.inp.out_dos != 0 || PARAM.inp.out_band[0] != 0)
    {
        GlobalV::ofs_running << "\n\n\n\n";
        GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
                                ">>>>>>>>>>>>>>>>>>>>>>>>>"
                             << std::endl;
        GlobalV::ofs_running << " |                                            "
                                "                        |"
                             << std::endl;
        GlobalV::ofs_running << " | Post-processing of data:                   "
                                "                        |"
                             << std::endl;
        GlobalV::ofs_running << " | DOS (density of states) and bands will be "
                                "output here.             |"
                             << std::endl;
        GlobalV::ofs_running << " | If atomic orbitals are used, Mulliken "
                                "charge analysis can be done. |"
                             << std::endl;
        GlobalV::ofs_running << " | Also the .bxsf file containing fermi "
                                "surface information can be    |"
                             << std::endl;
        GlobalV::ofs_running << " | done here.                                 "
                                "                        |"
                             << std::endl;
        GlobalV::ofs_running << " |                                            "
                                "                        |"
                             << std::endl;
        GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
                                "<<<<<<<<<<<<<<<<<<<<<<<<<"
                             << std::endl;
        GlobalV::ofs_running << "\n\n\n\n";
    }

    int nspin0 = 1;
    if (GlobalV::NSPIN == 2)
    {
        nspin0 = 2;
    }
    //! print occupation in istate.info
    ModuleIO::write_istate_info(this->pelec->ekb, this->pelec->wg, this->kv, &(GlobalC::Pkpoints));

    //! compute density of states
    if (PARAM.inp.out_dos)
    {
        ModuleIO::write_dos_pw(this->pelec->ekb,
                               this->pelec->wg,
                               this->kv,
                               PARAM.inp.dos_edelta_ev,
                               PARAM.inp.dos_scale,
                               PARAM.inp.dos_sigma);

        if (nspin0 == 1)
        {
            GlobalV::ofs_running << " Fermi energy is " << this->pelec->eferm.ef << " Rydberg" << std::endl;
        }
        else if (nspin0 == 2)
        {
            GlobalV::ofs_running << " Fermi energy (spin = 1) is " << this->pelec->eferm.ef_up << " Rydberg"
                                 << std::endl;
            GlobalV::ofs_running << " Fermi energy (spin = 2) is " << this->pelec->eferm.ef_dw << " Rydberg"
                                 << std::endl;
        }
    }

    if (PARAM.inp.out_band[0]) // pengfei 2014-10-13
    {
        for (int is = 0; is < nspin0; is++)
        {
            std::stringstream ss2;
            ss2 << GlobalV::global_out_dir << "BANDS_" << is + 1 << ".dat";
            GlobalV::ofs_running << "\n Output bands in file: " << ss2.str() << std::endl;
            ModuleIO::nscf_band(is,
                                ss2.str(),
                                GlobalV::NBANDS,
                                0.0,
                                PARAM.inp.out_band[1],
                                this->pelec->ekb,
                                this->kv,
                                &(GlobalC::Pkpoints));
        }
    }

    if (GlobalV::BASIS_TYPE == "pw" && winput::out_spillage) // xiaohui add 2013-09-01
    {
        //  calculate spillage value.

        // ! Print out overlap before spillage optimization to generate atomic
        // orbitals
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
                numerical_basis.output_overlap(this->psi[0], this->sf, this->kv, this->pw_wfc, GlobalC::ucell, i);
            }
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "BASIS OVERLAP (Q and S) GENERATION.");
        }
    }

    //! Print out wave functions in real space
    if (this->wf.out_wfc_r == 1) // Peize Lin add 2021.11.21
    {
        ModuleIO::write_psi_r_1(this->psi[0], this->pw_wfc, "wfc_realspace", true, this->kv);
    }

    //! Use Kubo-Greenwood method to compute conductivities
    if (PARAM.inp.cal_cond)
    {
        EleCond elec_cond(&GlobalC::ucell, &this->kv, this->pelec, this->pw_wfc, this->psi, &GlobalC::ppcell);
        elec_cond.KG(PARAM.inp.cond_smear,
                     PARAM.inp.cond_fwhm,
                     PARAM.inp.cond_wcut,
                     PARAM.inp.cond_dw,
                     PARAM.inp.cond_dt,
                     PARAM.inp.cond_nonlocal,
                     this->pelec->wg);
    }
}

template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace ModuleESolver
