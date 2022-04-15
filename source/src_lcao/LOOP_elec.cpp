#include "LOOP_elec.h"
#include "LCAO_diago.h"
#include "../src_pw/global.h"
#include "../src_pw/symmetry_rho.h"
#include "input_update.h"
#include "../src_io/chi0_hilbert.h"
#include "LCAO_evolve.h"
#include "dftu.h"
//
#include "../module_neighbor/sltk_atom_arrange.h"
#include "../src_io/istate_charge.h"
#include "../src_io/istate_envelope.h"
#include "ELEC_scf.h"
#include "ELEC_nscf.h"
#include "ELEC_cbands_gamma.h"
#include "ELEC_cbands_k.h"
#include "ELEC_evolve.h"
//
#include "../src_ri/exx_abfs.h"
#include "../src_ri/exx_opt_orb.h"
#include "../src_pw/vdwd2.h"
#include "../src_pw/vdwd3.h"
#include "../module_base/timer.h"
#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h"
#endif

void LOOP_elec::solve_elec_stru(const int& istep,
    Record_adj &ra,
    Local_Orbital_Charge& loc,
    Local_Orbital_wfc& lowf,
    LCAO_Hamilt& uhm_in)
{
    ModuleBase::TITLE("LOOP_elec","solve_elec_stru"); 
    ModuleBase::timer::tick("LOOP_elec", "solve_elec_stru");

    this->UHM = &uhm_in;

	// prepare HS matrices, prepare grid integral
	this->set_matrix_grid(ra);
	// density matrix extrapolation and prepare S,T,VNL matrices 
	this->before_solver(istep, loc, lowf);
	// do self-interaction calculations / nscf/ tddft, etc. 
	this->solver(istep, loc, lowf);

    ModuleBase::timer::tick("LOOP_elec","solve_elec_stru"); 
	return;
}


void LOOP_elec::set_matrix_grid(Record_adj &ra)
{
    ModuleBase::TITLE("LOOP_elec","set_matrix_grid"); 
    ModuleBase::timer::tick("LOOP_elec","set_matrix_grid"); 

	// (1) Find adjacent atoms for each atom.
	GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(
		GlobalV::ofs_running,
		GlobalV::OUT_LEVEL,
		GlobalC::ORB.get_rcutmax_Phi(), 
		GlobalC::ucell.infoNL.get_rcutmax_Beta(), 
		GlobalV::GAMMA_ONLY_LOCAL);

	atom_arrange::search(
		GlobalV::SEARCH_PBC,
		GlobalV::ofs_running,
		GlobalC::GridD, 
		GlobalC::ucell, 
		GlobalV::SEARCH_RADIUS, 
		GlobalV::test_atom_input);

	//ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"SEARCH ADJACENT ATOMS");

	// (3) Periodic condition search for each grid.
	GlobalC::GridT.set_pbc_grid(
			GlobalC::pw.ncx, GlobalC::pw.ncy, GlobalC::pw.ncz,
			GlobalC::pw.bx, GlobalC::pw.by, GlobalC::pw.bz,
			GlobalC::pw.nbx, GlobalC::pw.nby, GlobalC::pw.nbz,
			GlobalC::pw.nbxx, GlobalC::pw.nbzp_start, GlobalC::pw.nbzp);

    // (2)For each atom, calculate the adjacent atoms in different cells
    // and allocate the space for H(R) and S(R).
    // If k point is used here, allocate HlocR after atom_arrange.
    Parallel_Orbitals* pv = this->UHM->LM->ParaV;
    ra.for_2d(*pv, GlobalV::GAMMA_ONLY_LOCAL);
	if(!GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->UHM->LM->allocate_HS_R(pv->nnr);
#ifdef __DEEPKS
		GlobalC::ld.allocate_V_deltaR(pv->nnr);
#endif

		// need to first calculae lgd.
		// using GlobalC::GridT.init.
		GlobalC::GridT.cal_nnrg();
	}

    ModuleBase::timer::tick("LOOP_elec","set_matrix_grid"); 
	return;
}


void LOOP_elec::before_solver(const int& istep,
    Local_Orbital_Charge& loc,
    Local_Orbital_wfc& lowf)
{
    ModuleBase::TITLE("LOOP_elec","before_solver"); 
    ModuleBase::timer::tick("LOOP_elec","before_solver"); 

	// set the augmented orbitals index.
	// after ParaO and GridT, 
	// this information is used to calculate
	// the force.

	// init density kernel and wave functions.
	loc.allocate_dm_wfc(GlobalC::GridT, lowf);

	//======================================
	// do the charge extrapolation before the density matrix is regenerated.
	// mohan add 2011-04-08
	// because once atoms are moving out of this processor,
	// the density matrix will not std::map the new atomic configuration,
	//======================================
	// THIS IS A BUG, BECAUSE THE INDEX GlobalC::GridT.trace_lo
	// HAS BEEN REGENERATED, SO WE NEED TO
	// REALLOCATE DENSITY MATRIX FIRST, THEN READ IN DENSITY MATRIX,
	// AND USE DENSITY MATRIX TO DO RHO GlobalV::CALCULATION.-- mohan 2013-03-31
	//======================================
	if(GlobalC::pot.chg_extrap=="dm" && istep>1)//xiaohui modify 2015-02-01
	{
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[is], GlobalC::pw.nrxx);
			std::stringstream ssd;
			ssd << GlobalV::global_out_dir << "SPIN" << is + 1 << "_DM" ;
			// reading density matrix,
			loc.read_dm(is, ssd.str() );
		}

		// calculate the charge density
		if(GlobalV::GAMMA_ONLY_LOCAL)
		{
			this->UHM->GG.cal_rho(loc.DM);
		}
		else
		{
			this->UHM->GK.cal_rho_k(loc.DM_R);
		}

		// renormalize the charge density
		GlobalC::CHR.renormalize_rho();

		// initialize the potential
		GlobalC::pot.init_pot( istep-1, GlobalC::pw.strucFac );
	}


	// (9) compute S, T, Vnl, Vna matrix.
    this->UHM->set_lcao_matrices();

#ifdef __DEEPKS
    //for each ionic step, the overlap <psi|alpha> must be rebuilt
    //since it depends on ionic positions
    if (GlobalV::deepks_setorb)
    {
        const Parallel_Orbitals* pv = this->UHM->LM->ParaV;
        //build and save <psi(0)|alpha(R)> at beginning
        GlobalC::ld.build_psialpha(GlobalV::CAL_FORCE,
			GlobalC::ucell,
			GlobalC::ORB,
			GlobalC::GridD,
			pv->trace_loc_row,
			pv->trace_loc_col,
			GlobalC::UOT);

		if(GlobalV::deepks_out_unittest)
		{
			GlobalC::ld.check_psialpha(GlobalV::CAL_FORCE,
					GlobalC::ucell,
					GlobalC::ORB,
					GlobalC::GridD,
					pv->trace_loc_row,
					pv->trace_loc_col,
					GlobalC::UOT);
		}
    }
#endif

    ModuleBase::timer::tick("LOOP_elec","before_solver"); 
	return;
}

void LOOP_elec::solver(const int& istep,
    Local_Orbital_Charge& loc,
    Local_Orbital_wfc& lowf)
{
    ModuleBase::TITLE("LOOP_elec","solver"); 
    ModuleBase::timer::tick("LOOP_elec","solver"); 

	// self consistent calculations for electronic ground state
	if (GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="md"
			|| GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax") //pengfei 2014-10-13
	{
	#ifdef __MPI
		//Peize Lin add 2016-12-03
		if( Exx_Global::Hybrid_Type::HF==GlobalC::exx_lcao.info.hybrid_type 
			|| Exx_Global::Hybrid_Type::PBE0==GlobalC::exx_lcao.info.hybrid_type 
			|| Exx_Global::Hybrid_Type::HSE==GlobalC::exx_lcao.info.hybrid_type )
		{
			GlobalC::exx_lcao.cal_exx_ions(*lowf.ParaV);
		}

		// No exx
		if( Exx_Global::Hybrid_Type::No==GlobalC::exx_global.info.hybrid_type  )
		{
			ELEC_scf es;
            es.scf(istep - 1, loc, lowf, *this->UHM);
        }
		else if( Exx_Global::Hybrid_Type::Generate_Matrix == GlobalC::exx_global.info.hybrid_type )
		{
			Exx_Opt_Orb exx_opt_orb;
			exx_opt_orb.generate_matrix();
		}
		else    // Peize Lin add 2016-12-03
		{
		#endif // __MPI
			ELEC_scf es;
            es.scf(istep - 1, loc, lowf, *this->UHM);
		#ifdef __MPI
            if (GlobalC::exx_global.info.separate_loop)
			{
				for( size_t hybrid_step=0; hybrid_step!=GlobalC::exx_global.info.hybrid_step; ++hybrid_step )
				{
					XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].xc_func);
					GlobalC::exx_lcao.cal_exx_elec(loc, lowf.wfc_k_grid);
					
					ELEC_scf es;
					es.scf(istep-1, loc, lowf, *this->UHM);
					if(ELEC_scf::iter==1)     // exx converge
					{
						break;
					}
				}
			}
			else
			{
				XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].xc_func);
				ELEC_scf es;
				es.scf(istep-1, loc, lowf, *this->UHM);
			}
		}
		#endif // __MPI
	}
	else if (GlobalV::CALCULATION=="nscf")
	{
		ELEC_nscf::nscf(*this->UHM, loc.dm_gamma, loc.dm_k, lowf);
	}
	else if (GlobalV::CALCULATION=="istate")
	{
		IState_Charge ISC(lowf.wfc_gamma, loc);
		ISC.begin(this->UHM->GG);
	}
	else if (GlobalV::CALCULATION=="ienvelope")
	{
		IState_Envelope IEP;
        IEP.begin(lowf, this->UHM->GG, INPUT.out_wfc_pw, GlobalC::wf.out_wfc_r);
    }
	else
	{
		ModuleBase::WARNING_QUIT("LOOP_elec::solver","CALCULATION type not supported");
	}

    ModuleBase::timer::tick("LOOP_elec","solver"); 
	return;
}

