#include "esolver_ks_lcao.h"

#include "module_base/formatter.h"
#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_io/berryphase.h"
#include "module_io/cube_io.h"
#include "module_io/dos_nao.h"
#include "module_io/io_dmk.h"
#include "module_io/io_npz.h"
#include "module_io/nscf_band.h"
#include "module_io/output_dmk.h"
#include "module_io/output_log.h"
#include "module_io/output_mat_sparse.h"
#include "module_io/output_mulliken.h"
#include "module_io/output_sk.h"
#include "module_io/to_qo.h"
#include "module_io/to_wannier90_lcao.h"
#include "module_io/to_wannier90_lcao_in_pw.h"
#include "module_io/write_HS.h"
#include "module_io/write_dmr.h"
#include "module_io/write_elecstat_pot.h"
#include "module_io/write_istate_info.h"
#include "module_io/write_proj_band_lcao.h"
#include "module_io/write_wfc_nao.h"
#include "module_parameter/parameter.h"


//--------------temporary----------------------------
#include "module_base/global_function.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/cal_ux.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h" // need DeePKS_init
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"


//mohan add 20250302
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/ekinetic_new.h"


#include <memory>
#ifdef __EXX
#include "module_io/restart_exx_csr.h"
#include "module_ri/RPA_LRI.h"
#endif

#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#include "module_hamilt_lcao/module_deepks/LCAO_deepks_interface.h"
#endif
//-----force& stress-------------------
#include "module_hamilt_lcao/hamilt_lcaodft/FORCE_STRESS.h"

//-----HSolver ElecState Hamilt--------
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hsolver/hsolver_lcao.h"
// function used by deepks
// #include "module_elecstate/cal_dm.h"
//---------------------------------------------------

// test RDMFT
#include "module_rdmft/rdmft.h"

#include <iostream>

namespace ModuleESolver
{

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "after_scf");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "after_scf");

    //------------------------------------------------------------------
    //! 1) calculate the kinetic energy density tau in LCAO basis
    //！sunliang 2024-09-18
    //------------------------------------------------------------------
    if (PARAM.inp.out_elf[0] > 0)
    {
        assert(this->psi != nullptr);
        this->pelec->cal_tau(*(this->psi));
    }

    //------------------------------------------------------------------
    //! 2) call after_scf() of ESolver_KS
    //------------------------------------------------------------------
    ESolver_KS<TK>::after_scf(ucell, istep, conv_esolver);


    //------------------------------------------------------------------
    //! 3) write density matrix for sparse matrix in LCAO basis
    //------------------------------------------------------------------
    ModuleIO::write_dmr(dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->get_DMR_vector(),
                        this->pv,
                        PARAM.inp.out_dm1,
                        false,
                        PARAM.inp.out_app_flag,
                        ucell.get_iat2iwt(),
                        &ucell.nat,
                        istep);

    //------------------------------------------------------------------
    //! 4) write density matrix in LCAO basis
    //------------------------------------------------------------------
    if (PARAM.inp.out_dm)
    {
        std::vector<double> efermis(PARAM.inp.nspin == 2 ? 2 : 1);
        for (int ispin = 0; ispin < efermis.size(); ispin++)
        {
            efermis[ispin] = this->pelec->eferm.get_efval(ispin);
        }
        const int precision = 3;
        ModuleIO::write_dmk(dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->get_DMK_vector(),
                            precision,
                            efermis,
                            &(ucell),
                            this->pv);
    }

#ifdef __EXX
    //------------------------------------------------------------------
    //! 5) write Hexx matrix in LCAO basis
    // (see `out_chg` in docs/advanced/input_files/input-main.md)
    //------------------------------------------------------------------
    if (PARAM.inp.calculation != "nscf")
    {
        if (GlobalC::exx_info.info_global.cal_exx && 
            PARAM.inp.out_chg[0] && 
            istep % PARAM.inp.out_interval == 0) // Peize Lin add if 2022.11.14
        {
			const std::string file_name_exx = PARAM.globalv.global_out_dir + 
				"HexxR" + 
				std::to_string(GlobalV::MY_RANK);
            if (GlobalC::exx_info.info_ri.real_number)
            {
                ModuleIO::write_Hexxs_csr(file_name_exx, ucell, this->exd->get_Hexxs());
            }
            else
            {
                ModuleIO::write_Hexxs_csr(file_name_exx, ucell, this->exc->get_Hexxs());
            }
        }
    }
#endif

	//------------------------------------------------------------------
	// 6) write Hamiltonian and Overlap matrix in LCAO basis
	//------------------------------------------------------------------
	for (int ik = 0; ik < this->kv.get_nks(); ++ik)
	{
		if (PARAM.inp.out_mat_hs[0])
		{
			this->p_hamilt->updateHk(ik);
		}
		bool bit = false; // LiuXh, 2017-03-21
		// if set bit = true, there would be error in soc-multi-core
		// calculation, noted by zhengdy-soc
		if (this->psi != nullptr && (istep % PARAM.inp.out_interval == 0))
		{
			hamilt::MatrixBlock<TK> h_mat;
			hamilt::MatrixBlock<TK> s_mat;

			this->p_hamilt->matrix(h_mat, s_mat);

			if (PARAM.inp.out_mat_hs[0])
			{
				ModuleIO::save_mat(istep,
						h_mat.p,
						PARAM.globalv.nlocal,
						bit,
						PARAM.inp.out_mat_hs[1],
						1,
						PARAM.inp.out_app_flag,
						"H",
						"data-" + std::to_string(ik),
						this->pv,
						GlobalV::DRANK);
				ModuleIO::save_mat(istep,
						s_mat.p,
						PARAM.globalv.nlocal,
						bit,
						PARAM.inp.out_mat_hs[1],
						1,
						PARAM.inp.out_app_flag,
						"S",
						"data-" + std::to_string(ik),
						this->pv,
						GlobalV::DRANK);
			}
		}
	}

    //------------------------------------------------------------------
    // 7) write electronic wavefunctions in LCAO basis
    //------------------------------------------------------------------
    if (elecstate::ElecStateLCAO<TK>::out_wfc_lcao && (istep % PARAM.inp.out_interval == 0))
    {
        ModuleIO::write_wfc_nao(elecstate::ElecStateLCAO<TK>::out_wfc_lcao,
                                this->psi[0],
                                this->pelec->ekb,
                                this->pelec->wg,
                                this->pelec->klist->kvec_c,
                                this->pv,
                                istep);
    }

    //------------------------------------------------------------------
    //! 8) write DeePKS information in LCAO basis
    //------------------------------------------------------------------
#ifdef __DEEPKS
    if (this->psi != nullptr && (istep % PARAM.inp.out_interval == 0))
    {
        hamilt::HamiltLCAO<TK, TR>* p_ham_deepks = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);
        std::shared_ptr<LCAO_Deepks<TK>> ld_shared_ptr(&ld, [](LCAO_Deepks<TK>*) {});
        LCAO_Deepks_Interface<TK, TR> deepks_interface(ld_shared_ptr);

        deepks_interface.out_deepks_labels(this->pelec->f_en.etot,
                              this->pelec->klist->get_nks(),
                              ucell.nat,
                              PARAM.globalv.nlocal,
                              this->pelec->ekb,
                              this->pelec->klist->kvec_d,
                              ucell,
                              orb_,
                              this->gd,
                              &(this->pv),
                              *(this->psi),
                              dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                              p_ham_deepks,
                              GlobalV::MY_RANK);
    }
#endif


    //------------------------------------------------------------------
    //! 9) Perform RDMFT calculations
    // rdmft, added by jghan, 2024-10-17
    //------------------------------------------------------------------
    if (PARAM.inp.rdmft == true)
    {
        ModuleBase::matrix occ_num(this->pelec->wg);
        for (int ik = 0; ik < occ_num.nr; ++ik)
        {
            for (int inb = 0; inb < occ_num.nc; ++inb)
            {
                occ_num(ik, inb) /= this->kv.wk[ik];
            }
        }
        this->rdmft_solver.update_elec(ucell, occ_num, *(this->psi));

        //! initialize the gradients of Etotal with respect to occupation numbers and wfc,
        //! and set all elements to 0.
        //! dedocc = d E/d Occ_Num
        ModuleBase::matrix dedocc(this->pelec->wg.nr, this->pelec->wg.nc, true);

        //! dedwfc = d E/d wfc
        psi::Psi<TK> dedwfc(this->psi->get_nk(), this->psi->get_nbands(), this->psi->get_nbasis(), this->kv.ngk, true);
        dedwfc.zero_out();

        double etot_rdmft = this->rdmft_solver.run(dedocc, dedwfc);
    }


#ifdef __EXX
    //------------------------------------------------------------------
    // 10) Write RPA information in LCAO basis
    //------------------------------------------------------------------
    if (PARAM.inp.rpa)
    {
        RPA_LRI<TK, double> rpa_lri_double(GlobalC::exx_info.info_ri);
        rpa_lri_double.cal_postSCF_exx(*dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                                       MPI_COMM_WORLD,
                                       ucell,
                                       this->kv,
                                       orb_);
        rpa_lri_double.init(MPI_COMM_WORLD, this->kv, orb_.cutoffs());
        rpa_lri_double.out_for_RPA(ucell, this->pv, *(this->psi), this->pelec);
    }
#endif


    //------------------------------------------------------------------
    // 11) write HR in npz format in LCAO basis
    //------------------------------------------------------------------
    if (PARAM.inp.out_hr_npz)
    {
        this->p_hamilt->updateHk(0); // first k point, up spin
        hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao
            = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(this->p_hamilt);
        std::string zipname = "output_HR0.npz";
        ModuleIO::output_mat_npz(ucell, zipname, *(p_ham_lcao->getHR()));

        if (PARAM.inp.nspin == 2)
        {
            this->p_hamilt->updateHk(this->kv.get_nks() / 2); // the other half of k points, down spin
            hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao
                = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(this->p_hamilt);
            zipname = "output_HR1.npz";
            ModuleIO::output_mat_npz(ucell, zipname, *(p_ham_lcao->getHR()));
        }
    }

	//------------------------------------------------------------------
	// 12) write density matrix in the 'npz' format in LCAO basis
	//------------------------------------------------------------------
	if (PARAM.inp.out_dm_npz)
	{
		const elecstate::DensityMatrix<TK, double>* dm
			= dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
		std::string zipname = "output_DM0.npz";
		ModuleIO::output_mat_npz(ucell, zipname, *(dm->get_DMR_pointer(1)));

		if (PARAM.inp.nspin == 2)
		{
			zipname = "output_DM1.npz";
			ModuleIO::output_mat_npz(ucell, zipname, *(dm->get_DMR_pointer(2)));
		}
	}

	//------------------------------------------------------------------
	//! 13) Print out information every 'out_interval' steps.
	//------------------------------------------------------------------
	if (PARAM.inp.calculation != "md" || istep % PARAM.inp.out_interval == 0)
	{
		//! Print out sparse matrix
		ModuleIO::output_mat_sparse(PARAM.inp.out_mat_hs2,
				PARAM.inp.out_mat_dh,
				PARAM.inp.out_mat_t,
				PARAM.inp.out_mat_r,
				istep,
				this->pelec->pot->get_effective_v(),
				this->pv,
				this->GK,
				two_center_bundle_,
				orb_,
				ucell,
				this->gd,
				this->kv,
				this->p_hamilt);

		//! Perform Mulliken charge analysis in LCAO basis
		if (PARAM.inp.out_mul)
		{
			ModuleIO::cal_mag(&(this->pv),
					this->p_hamilt,
					this->kv,
					this->pelec,
					this->two_center_bundle_,
					this->orb_,
					ucell,
					this->gd,
					istep,
					true);
		}
	}

    //------------------------------------------------------------------
    //! 14) Print out atomic magnetization in LCAO basis
    //! only when 'spin_constraint' is on.
    //------------------------------------------------------------------
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        sc.cal_mi_lcao(istep);
        sc.print_Mi(GlobalV::ofs_running);
        sc.print_Mag_Force(GlobalV::ofs_running);
    }

    //------------------------------------------------------------------
    //! 15) Print out kinetic matrix in LCAO basis
    //------------------------------------------------------------------
    if (PARAM.inp.out_mat_tk[0])
	{
		hamilt::HS_Matrix_K<TK> hsk(&pv, true);
		hamilt::HContainer<TR> hR(&pv);
		hamilt::Operator<TK>* ekinetic
			= new hamilt::EkineticNew<hamilt::OperatorLCAO<TK, TR>>(&hsk,
					this->kv.kvec_d,
					&hR,
					&ucell,
					orb_.cutoffs(),
					&this->gd,
					two_center_bundle_.kinetic_orb.get());

		const int nspin_k = (PARAM.inp.nspin == 2 ? 2 : 1);
		for (int ik = 0; ik < this->kv.get_nks() / nspin_k; ++ik)
		{
			ekinetic->init(ik);
			ModuleIO::save_mat(0,
					hsk.get_hk(),
					PARAM.globalv.nlocal,
					false,
					PARAM.inp.out_mat_tk[1],
					1,
					PARAM.inp.out_app_flag,
					"T",
					"data-" + std::to_string(ik),
					this->pv,
					GlobalV::DRANK);
		}

		delete ekinetic;
	}

    //------------------------------------------------------------------
    //! 16) wannier90 interface in LCAO basis 
    // added by jingan in 2018.11.7
    //------------------------------------------------------------------
    if (PARAM.inp.calculation == "nscf" && PARAM.inp.towannier90)
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Wave function to Wannier90");
        if (PARAM.inp.wannier_method == 1)
		{
			toWannier90_LCAO_IN_PW wan(PARAM.inp.out_wannier_mmn,
					PARAM.inp.out_wannier_amn,
					PARAM.inp.out_wannier_unk,
					PARAM.inp.out_wannier_eig,
					PARAM.inp.out_wannier_wvfn_formatted,
					PARAM.inp.nnkpfile,
					PARAM.inp.wannier_spin);
			wan.set_tpiba_omega(ucell.tpiba, ucell.omega);
			wan.calculate(ucell,
					this->pelec->ekb,
					this->pw_wfc,
					this->pw_big,
					this->sf,
					this->kv,
					this->psi,
					&(this->pv));
		}
		else if (PARAM.inp.wannier_method == 2)
		{
			toWannier90_LCAO wan(PARAM.inp.out_wannier_mmn,
					PARAM.inp.out_wannier_amn,
					PARAM.inp.out_wannier_unk,
					PARAM.inp.out_wannier_eig,
					PARAM.inp.out_wannier_wvfn_formatted,
					PARAM.inp.nnkpfile,
					PARAM.inp.wannier_spin,
					orb_);

			wan.calculate(ucell, this->gd, this->pelec->ekb, this->kv, *(this->psi), &(this->pv));
		}
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Wave function to Wannier90");
    }

    //------------------------------------------------------------------
    //! 17) berry phase calculations in LCAO basis 
    // added by jingan
    //------------------------------------------------------------------
    if (PARAM.inp.calculation == "nscf" && 
        berryphase::berry_phase_flag && 
        ModuleSymmetry::Symmetry::symm_flag != 1)
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Berry phase calculation");
        berryphase bp(&(this->pv));
        bp.lcao_init(ucell, this->gd, this->kv, this->GridT, orb_);
        // additional step before calling macroscopic_polarization 
        bp.Macroscopic_polarization(ucell, this->pw_wfc->npwk_max, this->psi, this->pw_rho, this->pw_wfc, this->kv);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Berry phase calculation");
    }

    //------------------------------------------------------------------
    //! 18) calculate quasi-orbitals in LCAO basis
    //------------------------------------------------------------------
	if (PARAM.inp.qo_switch)
	{
		toQO tqo(PARAM.inp.qo_basis, PARAM.inp.qo_strategy, PARAM.inp.qo_thr, PARAM.inp.qo_screening_coeff);
		tqo.initialize(PARAM.globalv.global_out_dir,
				PARAM.inp.pseudo_dir,
				PARAM.inp.orbital_dir,
				&ucell,
				this->kv.kvec_d,
				GlobalV::ofs_running,
				GlobalV::MY_RANK,
				GlobalV::NPROC);
		tqo.calculate();
	}

    //------------------------------------------------------------------
    //! 19) Clean up RA, which is used to serach for adjacent atoms
    //------------------------------------------------------------------
    if (!PARAM.inp.cal_force && !PARAM.inp.cal_stress)
    {
        RA.delete_grid();
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "after_scf");
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
