#ifndef CTRL_OUTPUT_LCAO_H 
#define CTRL_OUTPUT_LCAO_H 

#include <complex>

#include "module_io/write_dmr.h" // use ModuleIO::write_dmr() 

namespace ModuleIO
{

template <typename TK, typename TR>
void ctrl_output_lcao()
{
    ModuleBase::TITLE("ModuleIO", "ctrl_output_lcao");
    ModuleBase::timer::tick("ModuleIO", "ctrl_output_lcao");

/*
    const bool out_app_flag = PARAM.inp.out_app_flag;

	const auto* estate = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec);
	if(!estate)
	{
		WARNING_QUIT("ModuleIO::ctrl_output_lcao","pelec does not exist");
	}

	//------------------------------------------------------------------
	//! 1) write density matrix DM(R) for sparse matrix in LCAO basis
	//------------------------------------------------------------------
    if(PARAM.inp.out_dm1)
	{
		const auto& dmr_vector = estate->get_DM()->get_DMR_vector();
		ModuleIO::write_dmr(dmr_vector, pv,	out_app_flag,
				ucell.get_iat2iwt(), ucell.nat, istep);
	}
*/

/*
	//------------------------------------------------------------------
	//! 2) write density matrix DM(k) in LCAO basis
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
    //! 3) write Hexx matrix in LCAO basis
    // (see `out_chg` in docs/advanced/input_files/input-main.md)
    //------------------------------------------------------------------
    if (PARAM.inp.out_chg[0])
    {
        if (GlobalC::exx_info.info_global.cal_exx && PARAM.inp.calculation != "nscf") // Peize Lin add if 2022.11.14
        {
			const std::string file_name_exx = PARAM.globalv.global_out_dir 
				+ "HexxR" + std::to_string(GlobalV::MY_RANK);
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
    // 4) write H(k) and S(k) in NAO basis set 
    //------------------------------------------------------------------
	if (PARAM.inp.out_mat_hs[0])
	{
		if (this->psi != nullptr)
		{
			ModuleIO::write_hsk(PARAM.globalv.global_out_dir,
					PARAM.inp.nspin,
					this->kv.get_nks(), 
					this->kv.get_nkstot(), 
					this->kv.ik2iktot, 
                    this->kv.isk,
					this->p_hamilt, 
					this->pv, 
					PARAM.globalv.gamma_only_local,
					PARAM.inp.out_app_flag,
					istep,
					GlobalV::ofs_running);
		}
	}

    //------------------------------------------------------------------
    // 5) write electronic wavefunctions in LCAO basis
    //------------------------------------------------------------------
    if (elecstate::ElecStateLCAO<TK>::out_wfc_lcao)
    {
		ModuleIO::write_wfc_nao(elecstate::ElecStateLCAO<TK>::out_wfc_lcao,
				PARAM.inp.out_app_flag,
				this->psi[0],
				this->pelec->ekb,
				this->pelec->wg,
				this->kv.kvec_c,
				this->kv.ik2iktot,
				this->kv.get_nkstot(),
				this->pv,
				PARAM.inp.nspin,
				istep);
	}

    //------------------------------------------------------------------
    //! 6) write DeePKS information in LCAO basis
    //------------------------------------------------------------------
#ifdef __DEEPKS
    if (this->psi != nullptr)
    {
        hamilt::HamiltLCAO<TK, TR>* p_ham_deepks = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);
        std::shared_ptr<LCAO_Deepks<TK>> ld_shared_ptr(&ld, [](LCAO_Deepks<TK>*) {});
        LCAO_Deepks_Interface<TK, TR> deepks_interface(ld_shared_ptr);

        deepks_interface.out_deepks_labels(this->pelec->f_en.etot,
                                           this->kv.get_nks(),
                                           ucell.nat,
                                           PARAM.globalv.nlocal,
                                           this->pelec->ekb,
                                           this->kv.kvec_d,
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
    // 7) write HR in npz format in LCAO basis
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
    // 8) write density matrix in the 'npz' format in LCAO basis
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
    //! 9) Print out <phi_i|O|phi_j>, where O is H, S, dH, dS, T, r 
    //------------------------------------------------------------------
	ModuleIO::output_mat_sparse(PARAM.inp.out_mat_hs2,
			PARAM.inp.out_mat_dh,
			PARAM.inp.out_mat_ds,
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

    //------------------------------------------------------------------
	//! 10) Perform Mulliken charge analysis in LCAO basis
    //------------------------------------------------------------------
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

    //------------------------------------------------------------------
    //! 11) Print out kinetic matrix in LCAO basis
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

            const int out_label = 1; // 1: .txt, 2: .dat

			std::string t_fn = ModuleIO::filename_output(PARAM.globalv.global_out_dir,
					"tk","nao",ik,this->kv.ik2iktot,
					PARAM.inp.nspin,this->kv.get_nkstot(),
					out_label,PARAM.inp.out_app_flag,
                    PARAM.globalv.gamma_only_local,istep);

            ModuleIO::save_mat(istep,
                               hsk.get_hk(),
                               PARAM.globalv.nlocal,
                               false, // bit
                               PARAM.inp.out_mat_tk[1],
                               1, // true for upper triangle matrix
                               PARAM.inp.out_app_flag,
                               t_fn, 
                               this->pv,
                               GlobalV::DRANK);
        }

        delete ekinetic;
    }

    //------------------------------------------------------------------
    //! 12) calculate expectation of angular momentum operator in LCAO basis
    //------------------------------------------------------------------
    if (PARAM.inp.out_mat_l[0])
    {
        ModuleIO::AngularMomentumCalculator mylcalculator(
            PARAM.inp.orbital_dir,
            ucell,
            PARAM.inp.search_radius,
            PARAM.inp.test_deconstructor,
            PARAM.inp.test_grid,
            PARAM.inp.test_atom_input,
            PARAM.globalv.search_pbc,
            &GlobalV::ofs_running,
            GlobalV::MY_RANK
        );
        mylcalculator.calculate(PARAM.inp.suffix,
                                PARAM.globalv.global_out_dir,
                                ucell,
                                PARAM.inp.out_mat_l[1],
                                GlobalV::MY_RANK);
    }

    //------------------------------------------------------------------
    //! 13) Print out atomic magnetization in LCAO basis
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
    //! 14) Berry phase calculations in LCAO basis, added by jingan
    //------------------------------------------------------------------
    if (PARAM.inp.calculation == "nscf" && berryphase::berry_phase_flag && ModuleSymmetry::Symmetry::symm_flag != 1)
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Berry phase calculation");
        berryphase bp(&(this->pv));
        bp.lcao_init(ucell, this->gd, this->kv, this->GridT, orb_);
        // additional step before calling macroscopic_polarization
        bp.Macroscopic_polarization(ucell, this->pw_wfc->npwk_max, this->psi, this->pw_rho, this->pw_wfc, this->kv);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Berry phase calculation");
    }

    //------------------------------------------------------------------
    //! 15) Calculate quasi-orbitals in LCAO basis
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

#ifdef __EXX
    //------------------------------------------------------------------
    // 17) Write RPA information in LCAO basis
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
    //! 18) Perform RDMFT calculations, added by jghan, 2024-10-17
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
*/

    ModuleBase::timer::tick("ModuleIO", "ctrl_output_lcao");
}

}

template void ModuleIO::ctrl_output_lcao<double, double>();
template void ModuleIO::ctrl_output_lcao<std::complex<double>, double>();
template void ModuleIO::ctrl_output_lcao<std::complex<double>, std::complex<double>>();

#endif
