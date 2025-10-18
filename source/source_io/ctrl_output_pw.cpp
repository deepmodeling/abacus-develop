#include "source_io/ctrl_output_pw.h"

#include "source_io/write_wfc_pw.h" // use write_wfc_pw
#include "source_io/write_dos_pw.h" // use write_dos_pw
#include "source_io/to_wannier90_pw.h" // wannier90 interface
#include "source_pw/module_pwdft/onsite_projector.h" // use projector
#include "source_io/numerical_basis.h"
#include "source_io/numerical_descriptor.h"
#include "source_io/cal_ldos.h"
#include "source_io/berryphase.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_base/formatter.h"
#include "source_io/get_pchg_pw.h"
#include "source_io/get_wf_pw.h"
#include "source_pw/module_pwdft/elecond.h"

#ifdef __MLALGO
#include "source_io/write_mlkedf_descriptors.h"
#endif

void ModuleIO::ctrl_iter_pw(const int istep, 
		const int iter, 
		const double &conv_esolver,
		psi::Psi<std::complex<double>, base_device::DEVICE_CPU>* psi,
		const K_Vectors &kv,
		const ModulePW::PW_Basis_K *pw_wfc,
        const Input_para& inp)
{
    ModuleBase::TITLE("ModuleIO", "ctrl_iter_pw");
    ModuleBase::timer::tick("ModuleIO", "ctrl_iter_pw");
    //----------------------------------------------------------
    // 3) Print out electronic wavefunctions in pw basis
    // we only print information every few ionic steps
    //----------------------------------------------------------

    // if istep_in = -1, istep will not appear in file name
    // if iter_in = -1, iter will not appear in file name
    int istep_in = -1;
    int iter_in = -1;
    bool out_wfc_flag = false;
    if (inp.out_freq_ion>0) // default value of out_freq_ion is 0
    {
        if (istep % inp.out_freq_ion == 0)
        {
            if(iter % inp.out_freq_elec == 0 || iter == inp.scf_nmax || conv_esolver)
            {
                istep_in = istep;
                iter_in = iter;
                out_wfc_flag = true;
            }
        }
    }
    else if(iter == inp.scf_nmax || conv_esolver)
    {
        out_wfc_flag = true;
    }

    if (out_wfc_flag)
    {
        ModuleIO::write_wfc_pw(istep_in, iter_in,
                GlobalV::KPAR,
                GlobalV::MY_POOL,
                GlobalV::MY_RANK,
                inp.nbands,
                inp.nspin,
                PARAM.globalv.npol,
                GlobalV::RANK_IN_POOL,
                GlobalV::NPROC_IN_POOL,
                inp.out_wfc_pw,
                inp.ecutwfc,
                PARAM.globalv.global_out_dir,
                psi[0],
                kv,
                pw_wfc,
                GlobalV::ofs_running);
    }

	ModuleBase::timer::tick("ModuleIO", "ctrl_iter_pw");
	return;
}


template <typename T, typename Device>
void ModuleIO::ctrl_scf_pw(const int istep,
        UnitCell& ucell,
        elecstate::ElecState* pelec,
        const Charge &chr,
		const K_Vectors &kv,
		const ModulePW::PW_Basis_K *pw_wfc,
		const ModulePW::PW_Basis *pw_rho,
		const ModulePW::PW_Basis *pw_rhod,
		const ModulePW::PW_Basis_Big *pw_big,
        Setup_Psi<T, Device> &stp,
        const Device* ctx,
        const base_device::AbacusDevice_t &device,
        const Parallel_Grid &para_grid,
        const Input_para& inp)
{
    ModuleBase::TITLE("ModuleIO", "ctrl_scf_pw");
    ModuleBase::timer::tick("ModuleIO", "ctrl_scf_pw");

    // Transfer data from device (GPU) to host (CPU) in pw basis
    stp.copy_d2h(device);

    //----------------------------------------------------------
    //! 4) Compute density of states (DOS)
    //----------------------------------------------------------
    if (inp.out_dos)
    {
        bool out_dos_tmp = false;

        int istep_in = -1;

        // default value of out_freq_ion is 0
        if(inp.out_freq_ion==0)
        {
            out_dos_tmp = true;
        }
        else if (inp.out_freq_ion>0)
        {
            if (istep % inp.out_freq_ion == 0)
            {
                out_dos_tmp = true;
                istep_in=istep;
            }
            else
            {
                out_dos_tmp = false;
            }
        }
        else
        {
            out_dos_tmp = false;
        }

        // the above is only valid for KSDFT, not SDFT
        // Needs update in the near future
        if (inp.esolver_type == "sdft")
        {
            out_dos_tmp = false;
        }

        if(out_dos_tmp)
        {
            ModuleIO::write_dos_pw(ucell,
                    pelec->ekb,
                    pelec->wg,
                    kv,
                    inp.nbands,
                    istep_in,
                    pelec->eferm,
                    inp.dos_edelta_ev,
                    inp.dos_scale,
                    inp.dos_sigma,
                    GlobalV::ofs_running);
        }
    }


    //------------------------------------------------------------------
    // 5) calculate band-decomposed (partial) charge density in pw basis
    //------------------------------------------------------------------
    if (inp.out_pchg.size() > 0)
    {
        // update psi_d
        stp.update_psi_d();

        const int nbands = stp.psi_t->get_nbands();
        const int ngmc = chr.ngmc;

        ModuleIO::get_pchg_pw(inp.out_pchg,
                              nbands,
                              inp.nspin,
                              pw_rhod->nxyz,
                              ngmc,
                              &ucell,
                              stp.psi_d,
                              pw_rhod,
                              pw_wfc,
                              ctx,
                              para_grid,
                              PARAM.globalv.global_out_dir,
                              inp.if_separate_k,
                              kv,
                              GlobalV::KPAR,
                              GlobalV::MY_POOL,
                              &chr);
    }


    //------------------------------------------------------------------
    //! 6) calculate Wannier functions in pw basis
    //------------------------------------------------------------------
    if (inp.calculation == "nscf" && inp.towannier90)
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Wannier functions calculation");
        toWannier90_PW wan(inp.out_wannier_mmn,
                           inp.out_wannier_amn,
                           inp.out_wannier_unk,
                           inp.out_wannier_eig,
                           inp.out_wannier_wvfn_formatted,
                           inp.nnkpfile,
                           inp.wannier_spin);
        wan.set_tpiba_omega(ucell.tpiba, ucell.omega);
        wan.calculate(ucell, pelec->ekb, pw_wfc, pw_big, kv, stp.psi_cpu);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Wannier functions calculation");
    }


    //------------------------------------------------------------------
    //! 7) calculate Berry phase polarization in pw basis
    //------------------------------------------------------------------
    if (inp.calculation == "nscf" && berryphase::berry_phase_flag && ModuleSymmetry::Symmetry::symm_flag != 1)
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Berry phase polarization");
        berryphase bp;
        bp.Macroscopic_polarization(ucell, pw_wfc->npwk_max, stp.psi_cpu, pw_rho, pw_wfc, kv);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Berry phase polarization");
    }

    //------------------------------------------------------------------
    // 8) write spin constrian results in pw basis
    // spin constrain calculations, write atomic magnetization and magnetic force.
    //------------------------------------------------------------------
    if (inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<std::complex<double>>& sc
            = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
        sc.cal_mi_pw();
        sc.print_Mag_Force(GlobalV::ofs_running);
    }

    //------------------------------------------------------------------
    // 9) write onsite occupations for charge and magnetizations
    //------------------------------------------------------------------
    if (inp.onsite_radius > 0)
    { // float type has not been implemented
        auto* onsite_p = projectors::OnsiteProjector<double, Device>::get_instance();
        onsite_p->cal_occupations(reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(stp.psi_t),
                                  pelec->wg);
    }

    ModuleBase::timer::tick("ModuleIO", "ctrl_scf_pw");
    return;
}

template <typename T, typename Device>
void ModuleIO::ctrl_runner_pw(UnitCell& ucell, 
		elecstate::ElecState* pelec,	
        ModulePW::PW_Basis_K* pw_wfc,
        ModulePW::PW_Basis* pw_rho,
        ModulePW::PW_Basis* pw_rhod,
		Charge &chr,
        K_Vectors &kv,
        Setup_Psi<T, Device> &stp,
        Structure_Factor &sf,
        pseudopot_cell_vnl &ppcell,
		surchem &solvent,
        const Device* ctx,
        Parallel_Grid &para_grid,
        const Input_para& inp)
{
    ModuleBase::TITLE("ModuleIO", "ctrl_runner_pw");
    ModuleBase::timer::tick("ModuleIO", "ctrl_runner_pw");

	//----------------------------------------------------------
	//! 1) Compute LDOS
	//----------------------------------------------------------
	if (inp.out_ldos[0])
	{
		ModuleIO::cal_ldos_pw(reinterpret_cast<elecstate::ElecStatePW<std::complex<double>>*>(pelec),
			    stp.psi_cpu[0], para_grid, ucell);
	}

    //----------------------------------------------------------
    //! 2) Calculate the spillage value,
    //! which are used to generate numerical atomic orbitals
    //----------------------------------------------------------
    if (inp.basis_type == "pw" && inp.out_spillage)
    {
        // ! Print out overlap matrices
        if (inp.out_spillage <= 2)
        {
            for (int i = 0; i < inp.bessel_nao_rcuts.size(); i++)
            {
                if (GlobalV::MY_RANK == 0)
                {
                    std::cout << "update value: bessel_nao_rcut <- " << std::fixed << inp.bessel_nao_rcuts[i]
                              << " a.u." << std::endl;
                }
                Numerical_Basis numerical_basis;
                numerical_basis.output_overlap(stp.psi_cpu[0], sf, kv, pw_wfc, ucell, i);
            }
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "BASIS OVERLAP (Q and S) GENERATION.");
        }
    }

    //----------------------------------------------------------
    //! 3) Print out electronic wave functions in real space
    //----------------------------------------------------------
    if (inp.out_wfc_norm.size() > 0 || inp.out_wfc_re_im.size() > 0)
    {
        stp.update_psi_d();

        ModuleIO::get_wf_pw(inp.out_wfc_norm,
                            inp.out_wfc_re_im,
                            stp.psi_t->get_nbands(),
                            inp.nspin,
                            pw_rhod->nxyz,
                            &ucell,
                            stp.psi_d,
                            pw_wfc,
                            ctx,
                            para_grid,
                            PARAM.globalv.global_out_dir,
                            kv,
                            GlobalV::KPAR,
                            GlobalV::MY_POOL);
    }

    //----------------------------------------------------------
    //! 4) Use Kubo-Greenwood method to compute conductivities
    //----------------------------------------------------------
    if (inp.cal_cond)
    {
        using Real = typename GetTypeReal<T>::type;
        EleCond<Real, Device> elec_cond(&ucell, &kv, pelec, pw_wfc, stp.psi_t, &ppcell);
        elec_cond.KG(inp.cond_smear,
                     inp.cond_fwhm,
                     inp.cond_wcut,
                     inp.cond_dw,
                     inp.cond_dt,
                     inp.cond_nonlocal,
                     pelec->wg);
    }

#ifdef __MLALGO
    //----------------------------------------------------------
    //! 7) generate training data for ML-KEDF
    //----------------------------------------------------------
    if (inp.of_ml_gene_data == 1)
    {
        pelec->pot->update_from_charge(&chr, &ucell);

        ModuleIO::Write_MLKEDF_Descriptors write_mlkedf_desc;
        write_mlkedf_desc.cal_tool->set_para(chr.nrxx,
                                             inp.nelec,
                                             inp.of_tf_weight,
                                             inp.of_vw_weight,
                                             inp.of_ml_chi_p,
                                             inp.of_ml_chi_q,
                                             inp.of_ml_chi_xi,
                                             inp.of_ml_chi_pnl,
                                             inp.of_ml_chi_qnl,
                                             inp.of_ml_nkernel,
                                             inp.of_ml_kernel,
                                             inp.of_ml_kernel_scaling,
                                             inp.of_ml_yukawa_alpha,
                                             inp.of_ml_kernel_file,
                                             ucell.omega,
                                             pw_rho);

        write_mlkedf_desc.generateTrainData_KS(PARAM.globalv.global_mlkedf_descriptor_dir,
                                               stp.psi_t,
                                               pelec,
                                               pw_wfc,
                                               pw_rho,
                                               ucell,
                                               pelec->pot->get_effective_v(0));
    }
#endif

    ModuleBase::timer::tick("ModuleIO", "ctrl_runner_pw");
}

// complex<float> + CPU
template void ModuleIO::ctrl_scf_pw<std::complex<float>, base_device::DEVICE_CPU>(
    const int nstep,
	UnitCell& ucell,
    elecstate::ElecState* pelec,
    const Charge &chr,
    const K_Vectors &kv,
    const ModulePW::PW_Basis_K *pw_wfc,
    const ModulePW::PW_Basis *pw_rho,
    const ModulePW::PW_Basis *pw_rhod,
    const ModulePW::PW_Basis_Big *pw_big,
    Setup_Psi<std::complex<float>, base_device::DEVICE_CPU> &stp,
    const base_device::DEVICE_CPU* ctx,
    const base_device::AbacusDevice_t &device,
    const Parallel_Grid &para_grid,
    const Input_para& inp);

// complex<double> + CPU
template void ModuleIO::ctrl_scf_pw<std::complex<double>, base_device::DEVICE_CPU>(
    const int nstep,
    UnitCell& ucell,
    elecstate::ElecState* pelec,
    const Charge &chr,
    const K_Vectors &kv,
    const ModulePW::PW_Basis_K *pw_wfc,
    const ModulePW::PW_Basis *pw_rho,
    const ModulePW::PW_Basis *pw_rhod,
    const ModulePW::PW_Basis_Big *pw_big,
    Setup_Psi<std::complex<double>, base_device::DEVICE_CPU> &stp,
    const base_device::DEVICE_CPU* ctx,
    const base_device::AbacusDevice_t &device,
    const Parallel_Grid &para_grid,
    const Input_para& inp);

#if ((defined __CUDA) || (defined __ROCM))
// complex<float> + GPU
template void ModuleIO::ctrl_scf_pw<std::complex<float>, base_device::DEVICE_GPU>(
    const int nstep,
	UnitCell& ucell,
    elecstate::ElecState* pelec,
    const Charge &chr,
    const K_Vectors &kv,
    const ModulePW::PW_Basis_K *pw_wfc,
    const ModulePW::PW_Basis *pw_rho,
    const ModulePW::PW_Basis *pw_rhod,
    const ModulePW::PW_Basis_Big *pw_big,
    Setup_Psi<std::complex<float>, base_device::DEVICE_GPU> &stp,
    const base_device::DEVICE_GPU* ctx,
    const base_device::AbacusDevice_t &device,
    const Parallel_Grid &para_grid,
    const Input_para& inp);

// complex<double> + GPU
template void ModuleIO::ctrl_scf_pw<std::complex<double>, base_device::DEVICE_GPU>(
	const int nstep,
	UnitCell& ucell,
    elecstate::ElecState* pelec,
    const Charge &chr,
    const K_Vectors &kv,
    const ModulePW::PW_Basis_K *pw_wfc,
    const ModulePW::PW_Basis *pw_rho,
    const ModulePW::PW_Basis *pw_rhod,
    const ModulePW::PW_Basis_Big *pw_big,
    Setup_Psi<std::complex<double>, base_device::DEVICE_GPU> &stp,
    const base_device::DEVICE_GPU* ctx,
    const base_device::AbacusDevice_t &device,
    const Parallel_Grid &para_grid,
    const Input_para& inp);
#endif

// complex<float> + CPU
template void ModuleIO::ctrl_runner_pw<std::complex<float>, base_device::DEVICE_CPU>(
	UnitCell& ucell, 
	elecstate::ElecState* pelec,	
    ModulePW::PW_Basis_K* pw_wfc,
    ModulePW::PW_Basis* pw_rho,
	ModulePW::PW_Basis* pw_rhod,
	Charge &chr,
	K_Vectors &kv,
    Setup_Psi<std::complex<float>, base_device::DEVICE_CPU> &stp,
    Structure_Factor &sf,
    pseudopot_cell_vnl &ppcell,
	surchem &solvent,
    const base_device::DEVICE_CPU* ctx,
    Parallel_Grid &para_grid,
    const Input_para& inp);

// complex<double> + CPU
template void ModuleIO::ctrl_runner_pw<std::complex<double>, base_device::DEVICE_CPU>(
	UnitCell& ucell, 
	elecstate::ElecState* pelec,	
    ModulePW::PW_Basis_K* pw_wfc,
    ModulePW::PW_Basis* pw_rho,
	ModulePW::PW_Basis* pw_rhod,
	Charge &chr,
	K_Vectors &kv,
    Setup_Psi<std::complex<double>, base_device::DEVICE_CPU> &stp,
    Structure_Factor &sf,
    pseudopot_cell_vnl &ppcell,
	surchem &solvent,
    const base_device::DEVICE_CPU* ctx,
    Parallel_Grid &para_grid,
    const Input_para& inp);

#if ((defined __CUDA) || (defined __ROCM))
// complex<float> + GPU
template void ModuleIO::ctrl_runner_pw<std::complex<float>, base_device::DEVICE_GPU>(
	UnitCell& ucell, 
	elecstate::ElecState* pelec,	
    ModulePW::PW_Basis_K* pw_wfc,
    ModulePW::PW_Basis* pw_rho,
    ModulePW::PW_Basis* pw_rhod,
	Charge &chr,
	K_Vectors &kv,
    Setup_Psi<std::complex<float>, base_device::DEVICE_GPU> &stp,
    Structure_Factor &sf,
    pseudopot_cell_vnl &ppcell,
	surchem &solvent,
    const base_device::DEVICE_GPU* ctx,
    Parallel_Grid &para_grid,
    const Input_para& inp);

// complex<double> + GPU
template void ModuleIO::ctrl_runner_pw<std::complex<double>, base_device::DEVICE_GPU>(
	UnitCell& ucell, 
	elecstate::ElecState* pelec,	
    ModulePW::PW_Basis_K* pw_wfc,
    ModulePW::PW_Basis* pw_rho,
    ModulePW::PW_Basis* pw_rhod,
	Charge &chr,
	K_Vectors &kv,
    Setup_Psi<std::complex<double>, base_device::DEVICE_GPU> &stp,
    Structure_Factor &sf,
    pseudopot_cell_vnl &ppcell,
	surchem &solvent,
    const base_device::DEVICE_GPU* ctx,
    Parallel_Grid &para_grid,
    const Input_para& inp);
#endif
