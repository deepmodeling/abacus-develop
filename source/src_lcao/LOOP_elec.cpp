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
#include "LCAO_nnr.h"
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
#ifdef __DEEPKS
#include "LCAO_descriptor.h"
#endif

void LOOP_elec::solve_elec_stru(const int &istep)
{
    ModuleBase::TITLE("LOOP_elec","solve_elec_stru"); 
    ModuleBase::timer::tick("LOOP_elec","solve_elec_stru"); 

	// prepare HS matrices, prepare grid integral
	this->set_matrix_grid();
	// density matrix extrapolation and prepare S,T,VNL matrices 
	this->before_solver(istep);
	// do self-interaction calculations / nscf/ tddft, etc. 
	this->solver(istep);

    ModuleBase::timer::tick("LOOP_elec","solve_elec_stru"); 
	return;
}


void LOOP_elec::set_matrix_grid(void)
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

	// (2) If k point is used here, allocate HlocR after atom_arrange.
	if(!GlobalV::GAMMA_ONLY_LOCAL)
	{
		// For each atom, calculate the adjacent atoms in different cells
		// and allocate the space for H(R) and S(R).
		GlobalC::LNNR.cal_nnr();
		GlobalC::LM.allocate_HS_R(GlobalC::LNNR.nnr);

		// need to first calculae lgd.
		// using GlobalC::GridT.init.
		GlobalC::LNNR.cal_nnrg(GlobalC::GridT);
	}

    ModuleBase::timer::tick("LOOP_elec","set_matrix_grid"); 
	return;
}


void LOOP_elec::before_solver(const int &istep)
{
    ModuleBase::TITLE("LOOP_elec","before_solver"); 
    ModuleBase::timer::tick("LOOP_elec","before_solver"); 

	// set the augmented orbitals index.
	// after ParaO and GridT, 
	// this information is used to calculate
	// the force.
	GlobalC::LOWF.set_trace_aug(GlobalC::GridT);

	// init density kernel and wave functions.
	GlobalC::LOC.allocate_dm_wfc(GlobalC::GridT);

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
	if(GlobalC::pot.extra_pot=="dm" && istep>1)//xiaohui modify 2015-02-01
	{
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[is], GlobalC::pw.nrxx);
			std::stringstream ssd;
			ssd << GlobalV::global_out_dir << "SPIN" << is + 1 << "_DM" ;
			// reading density matrix,
			GlobalC::LOC.read_dm(is, ssd.str() );
		}

		// calculate the charge density
		if(GlobalV::GAMMA_ONLY_LOCAL)
		{
			GlobalC::UHM.GG.cal_rho(GlobalC::LOC.DM);
		}
		else
		{
			GlobalC::UHM.GK.cal_rho_k();
		}

		// renormalize the charge density
		GlobalC::CHR.renormalize_rho();

		// initialize the potential
		GlobalC::pot.init_pot( istep-1, GlobalC::pw.strucFac );
	}


	// (9) compute S, T, Vnl, Vna matrix.
	GlobalC::UHM.set_lcao_matrices();

    ModuleBase::timer::tick("LOOP_elec","before_solver"); 
	return;
}

void LOOP_elec::solver(const int &istep)
{
    ModuleBase::TITLE("LOOP_elec","solver"); 
    ModuleBase::timer::tick("LOOP_elec","solver"); 

	// self consistent calculations for electronic ground state
	if (GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="md"
			|| GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax") //pengfei 2014-10-13
	{
		//Peize Lin add 2016-12-03
		switch(GlobalC::exx_lcao.info.hybrid_type)
		{
			case Exx_Global::Hybrid_Type::HF:
			case Exx_Global::Hybrid_Type::PBE0:
			case Exx_Global::Hybrid_Type::HSE:
				GlobalC::exx_lcao.cal_exx_ions();
				break;
			case Exx_Global::Hybrid_Type::No:
			case Exx_Global::Hybrid_Type::Generate_Matrix:
				break;
			default:
				throw std::invalid_argument(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}

		// No exx
		if( Exx_Global::Hybrid_Type::No==GlobalC::exx_global.info.hybrid_type  )
		{
			ELEC_scf es;
			es.scf(istep-1);
		}
		else if( Exx_Global::Hybrid_Type::Generate_Matrix == GlobalC::exx_global.info.hybrid_type )
		{
			Exx_Opt_Orb exx_opt_orb;
			exx_opt_orb.generate_matrix();
		}
		else    // Peize Lin add 2016-12-03
		{
			ELEC_scf es;
			es.scf(istep-1);
			if( GlobalC::exx_global.info.separate_loop )
			{
				for( size_t hybrid_step=0; hybrid_step!=GlobalC::exx_global.info.hybrid_step; ++hybrid_step )
				{
					GlobalC::exx_global.info.set_xcfunc(GlobalC::xcf);
					GlobalC::exx_lcao.cal_exx_elec();
					
					ELEC_scf es;
					es.scf(istep-1);
					if(ELEC_scf::iter==1)     // exx converge
					{
						break;
					}
				}
			}
			else
			{
				GlobalC::exx_global.info.set_xcfunc(GlobalC::xcf);

				ELEC_scf es;
				es.scf(istep-1);
				
			}
		}
	}
	else if (GlobalV::CALCULATION=="nscf")
	{
		ELEC_nscf::nscf(GlobalC::UHM);
	}
	else if (GlobalV::CALCULATION=="istate")
	{
		IState_Charge ISC;
		ISC.begin();
	}
	else if (GlobalV::CALCULATION=="ienvelope")
	{
		IState_Envelope IEP;
		IEP.begin();
	}
	else
	{
		ModuleBase::WARNING_QUIT("LOOP_elec::solver","CALCULATION type not supported");
	}

    ModuleBase::timer::tick("LOOP_elec","solver"); 
	return;
}

