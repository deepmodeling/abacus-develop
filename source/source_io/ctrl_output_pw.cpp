/*
#include "source_io/ctrl_output_fp.h" // use ctrl_output_fp() 
#include "source_estate/module_charge/symmetry_rho.h" // use Symmetry_rho
#include "source_io/write_elecstat_pot.h" // use write_elecstat_pot 
#include "source_io/write_elf.h"
#include "cube_io.h"  // use write_vdata_palgrid
#include "source_hamilt/module_xc/xc_functional.h" // use XC_Functional

#ifdef USE_LIBXC
#include "source_io/write_libxc_r.h"
#endif
*/

namespace ModuleIO
{

void ctrl_iter_pw()
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
    if (PARAM.inp.out_freq_ion>0) // default value of out_freq_ion is 0
    {
        if (istep % PARAM.inp.out_freq_ion == 0)
        {
            if(iter % PARAM.inp.out_freq_elec == 0 || iter == PARAM.inp.scf_nmax || conv_esolver)
            {
                istep_in = istep;
                iter_in = iter;
                out_wfc_flag = true;
            }
        }
    }
    else if(iter == PARAM.inp.scf_nmax || conv_esolver)
    {
        out_wfc_flag = true;
    }

    if (out_wfc_flag)
    {
        ModuleIO::write_wfc_pw(istep_in, iter_in,
                GlobalV::KPAR,
                GlobalV::MY_POOL,
                GlobalV::MY_RANK,
                PARAM.inp.nbands,
                PARAM.inp.nspin,
                PARAM.globalv.npol,
                GlobalV::RANK_IN_POOL,
                GlobalV::NPROC_IN_POOL,
                PARAM.inp.out_wfc_pw,
                PARAM.inp.ecutwfc,
                PARAM.globalv.global_out_dir,
                this->psi[0],
                this->kv,
                this->pw_wfc,
                GlobalV::ofs_running);
    }

	ModuleBase::timer::tick("ModuleIO", "ctrl_iter_pw");
	return;
}


void ctrl_scf_pw()
{
    ModuleBase::TITLE("ModuleIO", "ctrl_scf_pw");
    ModuleBase::timer::tick("ModuleIO", "ctrl_scf_pw");

    //----------------------------------------------------------
    //! 4) Compute density of states (DOS)
    //----------------------------------------------------------
    if (PARAM.inp.out_dos)
    {
        bool out_dos_tmp = false;

        int istep_in = -1;

        // default value of out_freq_ion is 0
        if(PARAM.inp.out_freq_ion==0)
        {
            out_dos_tmp = true;
        }
        else if (PARAM.inp.out_freq_ion>0)
        {
            if (istep % PARAM.inp.out_freq_ion == 0)
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
        // this part needs update in the near future
        if (PARAM.inp.esolver_type == "sdft")
        {
            out_dos_tmp = false;
        }

        if(out_dos_tmp)
        {
            ModuleIO::write_dos_pw(ucell,
                    this->pelec->ekb,
                    this->pelec->wg,
                    this->kv,
                    PARAM.inp.nbands,
                    istep_in,
                    this->pelec->eferm,
                    PARAM.inp.dos_edelta_ev,
                    PARAM.inp.dos_scale,
                    PARAM.inp.dos_sigma,
                    GlobalV::ofs_running);
        }
    }


    //------------------------------------------------------------------
    // 5) calculate band-decomposed (partial) charge density in pw basis
    //------------------------------------------------------------------
    if (PARAM.inp.out_pchg.size() > 0)
    {
        if (this->__kspw_psi != nullptr && PARAM.inp.precision == "single")
        {
            delete reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->__kspw_psi);
        }

        // Refresh __kspw_psi
        this->__kspw_psi = PARAM.inp.precision == "single"
                               ? new psi::Psi<std::complex<double>, Device>(this->kspw_psi[0])
                               : reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->kspw_psi);

        ModuleIO::get_pchg_pw(PARAM.inp.out_pchg,
                              this->kspw_psi->get_nbands(),
                              PARAM.inp.nspin,
                              this->pw_rhod->nxyz,
                              this->chr.ngmc,
                              &ucell,
                              this->__kspw_psi,
                              this->pw_rhod,
                              this->pw_wfc,
                              this->ctx,
                              this->Pgrid,
                              PARAM.globalv.global_out_dir,
                              PARAM.inp.if_separate_k,
                              this->kv,
                              GlobalV::KPAR,
                              GlobalV::MY_POOL,
                              &this->chr);
    }


    //------------------------------------------------------------------
    //! 6) calculate Wannier functions in pw basis
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
    //! 7) calculate Berry phase polarization in pw basis
    //------------------------------------------------------------------
    if (PARAM.inp.calculation == "nscf" && berryphase::berry_phase_flag && ModuleSymmetry::Symmetry::symm_flag != 1)
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Berry phase polarization");
        berryphase bp;
        bp.Macroscopic_polarization(ucell, this->pw_wfc->npwk_max, this->psi, this->pw_rho, this->pw_wfc, this->kv);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Berry phase polarization");
    }

    //------------------------------------------------------------------
    // 8) write spin constrian results in pw basis
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
    // 9) write onsite occupations for charge and magnetizations
    //------------------------------------------------------------------
    if (PARAM.inp.onsite_radius > 0)
    { // float type has not been implemented
        auto* onsite_p = projectors::OnsiteProjector<double, Device>::get_instance();
        onsite_p->cal_occupations(reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->kspw_psi),
                                  this->pelec->wg);
    }

    ModuleBase::timer::tick("ModuleIO", "ctrl_scf_pw");
    return;
}


void ctrl_runner_pw(UnitCell& ucell, 
		elecstate::ElecState* pelec,	
        ModulePW::PW_Basis_Big* pw_big,
        ModulePW::PW_Basis* pw_rhod,
        Charge &chr,
        surchem &solvent,
        Parallel_Grid &para_grid,
		const int istep)
{
    ModuleBase::TITLE("ModuleIO", "ctrl_runner_pw");
    ModuleBase::timer::tick("ModuleIO", "ctrl_runner_pw");

	//----------------------------------------------------------
	//! 1) Compute LDOS
	//----------------------------------------------------------
	if (PARAM.inp.out_ldos[0])
	{
		ModuleIO::cal_ldos_pw(reinterpret_cast<elecstate::ElecStatePW<std::complex<double>>*>(this->pelec),
				this->psi[0],
				this->Pgrid,
				ucell);
	}

    //----------------------------------------------------------
    //! 2) Calculate the spillage value,
    //! which are used to generate numerical atomic orbitals
    //----------------------------------------------------------
    if (PARAM.inp.basis_type == "pw" && PARAM.inp.out_spillage)
    {
        // ! Print out overlap matrices
        if (PARAM.inp.out_spillage <= 2)
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
    //! 3) Print out electronic wave functions in real space
    //----------------------------------------------------------
    if (PARAM.inp.out_wfc_norm.size() > 0 || PARAM.inp.out_wfc_re_im.size() > 0)
    {
        if (this->__kspw_psi != nullptr && PARAM.inp.precision == "single")
        {
            delete reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->__kspw_psi);
        }

        // Refresh __kspw_psi
        this->__kspw_psi = PARAM.inp.precision == "single"
                               ? new psi::Psi<std::complex<double>, Device>(this->kspw_psi[0])
                               : reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->kspw_psi);

        ModuleIO::get_wf_pw(PARAM.inp.out_wfc_norm,
                            PARAM.inp.out_wfc_re_im,
                            this->kspw_psi->get_nbands(),
                            PARAM.inp.nspin,
                            this->pw_rhod->nxyz,
                            &ucell,
                            this->__kspw_psi,
                            this->pw_wfc,
                            this->ctx,
                            this->Pgrid,
                            PARAM.globalv.global_out_dir,
                            this->kv,
                            GlobalV::KPAR,
                            GlobalV::MY_POOL);
    }

    //----------------------------------------------------------
    //! 4) Use Kubo-Greenwood method to compute conductivities
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

#ifdef __MLALGO
    //----------------------------------------------------------
    //! 7) generate training data for ML-KEDF
    //----------------------------------------------------------
    if (PARAM.inp.of_ml_gene_data == 1)
    {
        this->pelec->pot->update_from_charge(&this->chr, &ucell);

        ModuleIO::Write_MLKEDF_Descriptors write_mlkedf_desc;
        write_mlkedf_desc.cal_tool->set_para(this->chr.nrxx,
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

        write_mlkedf_desc.generateTrainData_KS(PARAM.globalv.global_mlkedf_descriptor_dir,
                                               this->kspw_psi,
                                               this->pelec,
                                               this->pw_wfc,
                                               this->pw_rho,
                                               ucell,
                                               this->pelec->pot->get_effective_v(0));
    }
#endif

    ModuleBase::timer::tick("ModuleIO", "ctrl_runner_pw");
}

} // End ModuleIO
