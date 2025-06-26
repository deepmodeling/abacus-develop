#include <complex>

#include "source_estate/elecstate_lcao.h" // use elecstate::ElecState
#include "source_io/ctrl_output_lcao.h" // use ctrl_output_lcao() 
#include "source_lcao/hamilt_lcaodft/hamilt_lcao.h" // use hamilt::HamiltLCAO<TK, TR>
#include "source_hamilt/hamilt.h" // use Hamilt<T>  

// functions
#include "source_io/write_dos_lcao.h" // use ModuleIO::write_dos_lcao() 
#include "source_io/write_dmr.h" // use ModuleIO::write_dmr() 
#include "source_io/io_dmk.h" // use ModuleIO::write_dmk()
#include "source_io/write_HS.h" // use ModuleIO::write_hsk()
#include "source_io/write_wfc_nao.h" // use ModuleIO::write_wfc_nao() 
#include "source_io/output_mat_sparse.h" // use ModuleIO::output_mat_sparse() 
#include "source_io/output_mulliken.h" // use cal_mag()
#include "source_lcao/hamilt_lcaodft/operator_lcao/ekinetic_new.h" // use hamilt::EkineticNew
#include "source_io/cal_pLpR.h" // use AngularMomentumCalculator()
#include "source_lcao/module_deltaspin/spin_constrain.h" // use spinconstrain::SpinConstrain<TK>
#include "source_io/berryphase.h" // use berryphase
#include "source_io/to_wannier90_lcao.h" // use toWannier90_LCAO
#include "source_io/to_wannier90_lcao_in_pw.h" // use toWannier90_LCAO_IN_PW
#ifdef __MLALGO
#include "source_lcao/module_deepks/LCAO_deepks.h"
#include "source_lcao/module_deepks/LCAO_deepks_interface.h"
#endif
#ifdef __EXX
#include "module_ri/Exx_LRI_interface.h" // use EXX codes
#include "module_ri/RPA_LRI.h" // use RPA code
#endif
#include "module_rdmft/rdmft.h" // use RDMFT codes
#include "source_io/to_qo.h" // use toQO

namespace ModuleIO
{

template <typename TK, typename TR>
void ctrl_output_lcao(UnitCell& ucell, 
		K_Vectors& kv,
		elecstate::ElecStateLCAO<TK>* pelec, 
		Parallel_Orbitals& pv,
		Grid_Driver& gd,
		psi::Psi<TK>* psi,
		hamilt::HamiltLCAO<TK, TR>* p_hamilt,
		TwoCenterBundle &two_center_bundle,
		Gint_k &gk,
		LCAO_Orbitals &orb,
		const ModulePW::PW_Basis_K* pw_wfc, // for berryphase
		const ModulePW::PW_Basis* pw_rho, // for berryphase
		Grid_Technique &gt, // for berryphase
		const ModulePW::PW_Basis_Big* pw_big, // for Wannier90
		const Structure_Factor& sf, // for Wannier90
        rdmft::RDMFT<TK, TR> &rdmft_solver, // for RDMFT
#ifdef __MLALGO
		LCAO_Deepks<TK>& ld,
#endif
#ifdef __EXX
		Exx_LRI_Interface<TK, double>& exd,
		Exx_LRI_Interface<TK, std::complex<double>>& exc,
#endif
		const int istep)
{
    ModuleBase::TITLE("ModuleIO", "ctrl_output_lcao");
    ModuleBase::timer::tick("ModuleIO", "ctrl_output_lcao");

    const bool out_app_flag = PARAM.inp.out_app_flag;
    const bool gamma_only = PARAM.globalv.gamma_only_local;
    const int nspin = PARAM.inp.nspin;
    const std::string global_out_dir = PARAM.globalv.global_out_dir;

	//------------------------------------------------------------------
    // print out density of states (DOS)
	//------------------------------------------------------------------
    if (PARAM.inp.out_dos)
    {
        ModuleIO::write_dos_lcao(psi,
                p_hamilt,
                pv,
                ucell,
                kv,
                PARAM.inp.nbands,
                pelec->eferm,
                pelec->ekb,
                pelec->wg,
                PARAM.inp.dos_edelta_ev,
                PARAM.inp.dos_scale,
                PARAM.inp.dos_sigma,
                out_app_flag,
                istep,
                GlobalV::ofs_running);
    }

	//------------------------------------------------------------------
	//! 1) Output density matrix DM(R)
	//------------------------------------------------------------------
    if(PARAM.inp.out_dmr)
	{
		const auto& dmr_vector = pelec->get_DM()->get_DMR_vector();
		ModuleIO::write_dmr(dmr_vector, pv,	out_app_flag,
				ucell.get_iat2iwt(), ucell.nat, istep);
	}

	//------------------------------------------------------------------
	//! 2) Output density matrix DM(k)
	//------------------------------------------------------------------
	if (PARAM.inp.out_dmk)
	{
		std::vector<double> efermis(nspin == 2 ? 2 : 1);
		for (int ispin = 0; ispin < efermis.size(); ispin++)
		{
			efermis[ispin] = pelec->eferm.get_efval(ispin);
		}
		const int precision = 3;
		ModuleIO::write_dmk(pelec->get_DM()->get_DMK_vector(),
				precision, efermis, &(ucell), pv);
	}

    //------------------------------------------------------------------
    // 3) Output H(k) and S(k) matrices for each k-point
    //------------------------------------------------------------------
	if (PARAM.inp.out_mat_hs[0])
	{
		ModuleIO::write_hsk(global_out_dir,
				nspin,
				kv.get_nks(), 
				kv.get_nkstot(), 
				kv.ik2iktot, 
				kv.isk,
				p_hamilt, 
				pv, 
				gamma_only,
				out_app_flag,
				istep,
				GlobalV::ofs_running);
	}

    //------------------------------------------------------------------
    // 4) Output electronic wavefunctions Psi(k)
    //------------------------------------------------------------------
    if (elecstate::ElecStateLCAO<TK>::out_wfc_lcao)
    {
		ModuleIO::write_wfc_nao(elecstate::ElecStateLCAO<TK>::out_wfc_lcao,
				out_app_flag,
				psi[0],
				pelec->ekb,
				pelec->wg,
				kv.kvec_c,
				kv.ik2iktot,
				kv.get_nkstot(),
				pv,
				nspin,
				istep);
	}

    //------------------------------------------------------------------
    //! 5) Output DeePKS information
    //------------------------------------------------------------------
#ifdef __MLALGO
    // need control parameter
	hamilt::HamiltLCAO<TK, TR>* p_ham_deepks = p_hamilt;
	std::shared_ptr<LCAO_Deepks<TK>> ld_shared_ptr(&ld, [](LCAO_Deepks<TK>*) {});
	LCAO_Deepks_Interface<TK, TR> deepks_interface(ld_shared_ptr);

	deepks_interface.out_deepks_labels(pelec->f_en.etot,
			kv.get_nks(),
			ucell.nat,
			PARAM.globalv.nlocal,
			pelec->ekb,
			kv.kvec_d,
			ucell,
			orb,
			gd,
			&pv,
			*psi,
			pelec->get_DM(),
			p_ham_deepks,
            -1, // -1 when called in after scf
            true, // no used when after scf
			GlobalV::MY_RANK,
            GlobalV::ofs_running);
#endif

    //------------------------------------------------------------------
    //! 6) Output <phi_i|O|phi_j> matrices, where O can be chosen as
    //!    H, S, dH, dS, T, r. The format is CSR format. 
    //------------------------------------------------------------------
    hamilt::Hamilt<TK>* p_ham_tk = static_cast<hamilt::Hamilt<TK>*>(p_hamilt);

	ModuleIO::output_mat_sparse(PARAM.inp.out_mat_hs2,
			PARAM.inp.out_mat_dh,
			PARAM.inp.out_mat_ds,
			PARAM.inp.out_mat_t,
			PARAM.inp.out_mat_r,
			istep,
			pelec->pot->get_effective_v(),
			pv,
			gk,
			two_center_bundle,
			orb,
			ucell,
			gd,
			kv,
			p_ham_tk);

    //------------------------------------------------------------------
    //! 7) Output kinetic matrix
    //------------------------------------------------------------------
    if (PARAM.inp.out_mat_tk[0])
    {
        hamilt::HS_Matrix_K<TK> hsk(&pv, true);
        hamilt::HContainer<TR> hR(&pv);
        hamilt::Operator<TK>* ekinetic
			= new hamilt::EkineticNew<hamilt::OperatorLCAO<TK, TR>>(&hsk,
					kv.kvec_d,
					&hR,
					&ucell,
					orb.cutoffs(),
					&gd,
					two_center_bundle.kinetic_orb.get());

        const int nspin_k = (nspin == 2 ? 2 : 1);
        for (int ik = 0; ik < kv.get_nks() / nspin_k; ++ik)
        {
            ekinetic->init(ik);

            const int out_label = 1; // 1: .txt, 2: .dat

			std::string t_fn = ModuleIO::filename_output(global_out_dir,
					"tk","nao",ik,kv.ik2iktot,
					PARAM.inp.nspin,kv.get_nkstot(),
					out_label,out_app_flag,
                    gamma_only,istep);

            ModuleIO::save_mat(istep,
                               hsk.get_hk(),
                               PARAM.globalv.nlocal,
                               false, // bit
                               PARAM.inp.out_mat_tk[1],
                               1, // true for upper triangle matrix
                               PARAM.inp.out_app_flag,
                               t_fn, 
                               pv,
                               GlobalV::DRANK);
        }

        delete ekinetic;
    }

    //------------------------------------------------------------------
    //! 8) Output expectation of angular momentum operator
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
                                global_out_dir,
                                ucell,
                                PARAM.inp.out_mat_l[1],
                                GlobalV::MY_RANK);
    }

    //------------------------------------------------------------------
    //! 9) Output Mulliken charge
    //------------------------------------------------------------------
    if (PARAM.inp.out_mul)
    {
        ModuleIO::cal_mag(&pv,
                p_hamilt,
                kv,
                pelec,
                two_center_bundle,
                orb,
                ucell,
                gd,
                istep,
                true);
    }

    //------------------------------------------------------------------
    //! 10) Output atomic magnetization by using 'spin_constraint'
    //------------------------------------------------------------------
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        sc.cal_mi_lcao(istep);
        sc.print_Mi(GlobalV::ofs_running);
        sc.print_Mag_Force(GlobalV::ofs_running);
    }

    //------------------------------------------------------------------
    //! 11) Output Berry phase
    //------------------------------------------------------------------
    if (PARAM.inp.calculation == "nscf" && berryphase::berry_phase_flag && ModuleSymmetry::Symmetry::symm_flag != 1)
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Berry phase calculation");
        berryphase bp(&pv);
        bp.lcao_init(ucell, gd, kv, gt, orb);
        // additional step before calling macroscopic_polarization
        bp.Macroscopic_polarization(ucell, pw_wfc->npwk_max, psi, pw_rho, pw_wfc, kv);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Berry phase calculation");
    }

    //------------------------------------------------------------------
    //! 12) Wannier90 interface in LCAO basis
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
			wan.calculate(ucell,pelec->ekb,pw_wfc,pw_big,
					sf,kv,psi,&pv);
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
					orb);

			wan.calculate(ucell, gd, pelec->ekb, kv, *psi, &pv);
		}
		std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Wave function to Wannier90");
	}


#ifdef __EXX
    //------------------------------------------------------------------
    //! 13) Output Hexx matrix in LCAO basis
    // (see `out_chg` in docs/advanced/input_files/input-main.md)
    //------------------------------------------------------------------
    if (PARAM.inp.out_chg[0])
    {
        if (GlobalC::exx_info.info_global.cal_exx && PARAM.inp.calculation != "nscf") // Peize Lin add if 2022.11.14
        {
            const std::string file_name_exx = global_out_dir
                + "HexxR" + std::to_string(GlobalV::MY_RANK);
            if (GlobalC::exx_info.info_ri.real_number)
            {
                ModuleIO::write_Hexxs_csr(file_name_exx, ucell, exd.get_Hexxs());
            }
            else
            {
                ModuleIO::write_Hexxs_csr(file_name_exx, ucell, exc.get_Hexxs());
            }
        }
    }

    //------------------------------------------------------------------
    //! 14) Write RPA information in LCAO basis
    //------------------------------------------------------------------
    if (PARAM.inp.rpa)
    {
        RPA_LRI<TK, double> rpa_lri_double(GlobalC::exx_info.info_ri);
        rpa_lri_double.cal_postSCF_exx(*dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(pelec)->get_DM(),
                                       MPI_COMM_WORLD,
                                       ucell,
                                       kv,
                                       orb);
        rpa_lri_double.init(MPI_COMM_WORLD, kv, orb.cutoffs());
        rpa_lri_double.out_for_RPA(ucell, pv, *psi, pelec);
    }
#endif

    //------------------------------------------------------------------
    //! 18) Perform RDMFT calculations, added by jghan, 2024-10-17
    //------------------------------------------------------------------
    if (PARAM.inp.rdmft == true)
    {
        ModuleBase::matrix occ_num(pelec->wg);
        for (int ik = 0; ik < occ_num.nr; ++ik)
        {
            for (int inb = 0; inb < occ_num.nc; ++inb)
            {
                occ_num(ik, inb) /= kv.wk[ik];
            }
        }
        rdmft_solver.update_elec(ucell, occ_num, *psi);

        //! initialize the gradients of Etotal with respect to occupation numbers and wfc,
        //! and set all elements to 0.
        //! dedocc = d E/d Occ_Num
        ModuleBase::matrix dedocc(pelec->wg.nr, pelec->wg.nc, true);

        //! dedwfc = d E/d wfc
        psi::Psi<TK> dedwfc(psi->get_nk(), psi->get_nbands(), psi->get_nbasis(), kv.ngk, true);
        dedwfc.zero_out();

        double etot_rdmft = rdmft_solver.run(dedocc, dedwfc);
    }

    //------------------------------------------------------------------
    //! Output quasi orbitals
    //------------------------------------------------------------------
    if (PARAM.inp.qo_switch)
    {
        toQO tqo(PARAM.inp.qo_basis, PARAM.inp.qo_strategy, PARAM.inp.qo_thr, PARAM.inp.qo_screening_coeff);
        tqo.initialize(global_out_dir,
                       PARAM.inp.pseudo_dir,
                       PARAM.inp.orbital_dir,
                       &ucell,
                       kv.kvec_d,
                       GlobalV::ofs_running,
                       GlobalV::MY_RANK,
                       GlobalV::NPROC);
        tqo.calculate();
    }


    ModuleBase::timer::tick("ModuleIO", "ctrl_output_lcao");
}

} // End ModuleIO


// For gamma only
template void ModuleIO::ctrl_output_lcao<double, double>(UnitCell& ucell, 
		K_Vectors& kv,
		elecstate::ElecStateLCAO<double>* pelec, 
		Parallel_Orbitals& pv,
		Grid_Driver& gd,
		psi::Psi<double>* psi,
		hamilt::HamiltLCAO<double, double>* p_hamilt,
		TwoCenterBundle &two_center_bundle,
		Gint_k &gk,
		LCAO_Orbitals &orb,
		const ModulePW::PW_Basis_K* pw_wfc, // for berryphase
		const ModulePW::PW_Basis* pw_rho, // for berryphase
		Grid_Technique &gt, // for berryphase
		const ModulePW::PW_Basis_Big* pw_big, // for Wannier90
		const Structure_Factor& sf, // for Wannier90
		rdmft::RDMFT<double, double> &rdmft_solver, // for RDMFT
#ifdef __MLALGO
		LCAO_Deepks<double>& ld,
#endif
#ifdef __EXX
		Exx_LRI_Interface<double, double>& exd,
		Exx_LRI_Interface<double, std::complex<double>>& exc,
#endif
		const int istep);

// For multiple k-points
template void ModuleIO::ctrl_output_lcao<std::complex<double>, double>(UnitCell& ucell, 
		K_Vectors& kv,
		elecstate::ElecStateLCAO<std::complex<double>>* pelec, 
		Parallel_Orbitals& pv,
		Grid_Driver& gd,
		psi::Psi<std::complex<double>>* psi,
		hamilt::HamiltLCAO<std::complex<double>, double>* p_hamilt,
		TwoCenterBundle &two_center_bundle,
		Gint_k &gk,
		LCAO_Orbitals &orb,
		const ModulePW::PW_Basis_K* pw_wfc, // for berryphase
		const ModulePW::PW_Basis* pw_rho, // for berryphase
		Grid_Technique &gt, // for berryphase
		const ModulePW::PW_Basis_Big* pw_big, // for Wannier90
		const Structure_Factor& sf, // for Wannier90
		rdmft::RDMFT<std::complex<double>, double> &rdmft_solver, // for RDMFT
#ifdef __MLALGO
		LCAO_Deepks<std::complex<double>>& ld,
#endif
#ifdef __EXX
		Exx_LRI_Interface<std::complex<double>, double>& exd,
		Exx_LRI_Interface<std::complex<double>, std::complex<double>>& exc,
#endif
		const int istep);

template void ModuleIO::ctrl_output_lcao<std::complex<double>, std::complex<double>>(UnitCell& ucell, 
		K_Vectors& kv,
		elecstate::ElecStateLCAO<std::complex<double>>* pelec, 
		Parallel_Orbitals& pv,
		Grid_Driver& gd,
		psi::Psi<std::complex<double>>* psi,
		hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>* p_hamilt,
		TwoCenterBundle &two_center_bundle,
		Gint_k &gk,
		LCAO_Orbitals &orb,
		const ModulePW::PW_Basis_K* pw_wfc, // for berryphase
		const ModulePW::PW_Basis* pw_rho, // for berryphase
		Grid_Technique &gt, // for berryphase
		const ModulePW::PW_Basis_Big* pw_big, // for Wannier90
		const Structure_Factor& sf, // for Wannier90
		rdmft::RDMFT<std::complex<double>, std::complex<double>> &rdmft_solver, // for RDMFT
#ifdef __MLALGO
		LCAO_Deepks<std::complex<double>>& ld,
#endif
#ifdef __EXX
		Exx_LRI_Interface<std::complex<double>, double>& exd,
		Exx_LRI_Interface<std::complex<double>, std::complex<double>>& exc,
#endif
		const int istep);

