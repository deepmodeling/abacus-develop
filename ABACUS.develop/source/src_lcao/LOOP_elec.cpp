#include "LOOP_elec.h"
#include "LCAO_diago.h"
#include "src_pw/global.h"
#include "src_pw/symmetry_rho.h"
#include "input_update.h"
#include "src_io/chi0_hilbert.h"
#include "LCAO_evolve.h"
#include "dftu.h"
//
#include "module_neighbor/sltk_atom_arrange.h"
#include "src_lcao/LCAO_nnr.h"
#include "src_io/istate_charge.h"
#include "src_io/istate_envelope.h"
#include "ELEC_scf.h"
#include "ELEC_nscf.h"
#include "ELEC_cbands_gamma.h"
#include "ELEC_cbands_k.h"
#include "ELEC_evolve.h"
//
#include "src_ri/exx_abfs.h"
#include "src_ri/exx_opt_orb.h"
#include "src_pw/vdwd2.h"
#include "src_pw/vdwd3.h"


void LOOP_elec::solve_elec_stru(const int &istep)
{
    TITLE("LOOP_elec","solve_elec_stru"); 
    timer::tick("LOOP_elec","solve_elec_stru",'C'); 

	// prepare HS matrices, prepare grid integral
	this->set_matrix_grid();
	// density matrix extrapolation and prepare S,T,VNL matrices 
	this->before_solver(istep);
	// do self-interaction calculations / nscf/ tddft, etc. 
	this->solver(istep);

    timer::tick("LOOP_elec","solve_elec_stru",'C'); 
	return;
}


void LOOP_elec::set_matrix_grid(void)
{
    TITLE("LOOP_elec","set_matrix_grid"); 
    timer::tick("LOOP_elec","set_matrix_grid",'D'); 

	// (1) Find adjacent atoms for each atom.
	SEARCH_RADIUS = atom_arrange::set_sr_NL(
		ofs_running,
		OUT_LEVEL,
		ORB.get_rcutmax_Phi(), 
		ORB.get_rcutmax_Beta(), 
		GAMMA_ONLY_LOCAL);

	atom_arrange::search(
		SEARCH_PBC,
		ofs_running,
		GridD, 
		ucell, 
		SEARCH_RADIUS, 
		test_atom_input);

	//DONE(ofs_running,"SEARCH ADJACENT ATOMS");

	// (3) Periodic condition search for each grid.
	GridT.set_pbc_grid(
			pw.ncx, pw.ncy, pw.ncz,
			pw.bx, pw.by, pw.bz,
			pw.nbx, pw.nby, pw.nbz,
			pw.nbxx, pw.nbzp_start, pw.nbzp);

	// (2) If k point is used here, allocate HlocR after atom_arrange.
	if(!GAMMA_ONLY_LOCAL)
	{
		// For each atom, calculate the adjacent atoms in different cells
		// and allocate the space for H(R) and S(R).
		LNNR.cal_nnr();
		LM.allocate_HS_R(LNNR.nnr);

		// need to first calculae lgd.
		// using GridT.init.
		LNNR.cal_nnrg(GridT);
	}

    timer::tick("LOOP_elec","set_matrix_grid",'D'); 
	return;
}


void LOOP_elec::before_solver(const int &istep)
{
    TITLE("LOOP_elec","before_solver"); 
    timer::tick("LOOP_elec","before_solver",'D'); 

	// set the augmented orbitals index.
	// after ParaO and GridT, 
	// this information is used to calculate
	// the force.
	LOWF.set_trace_aug(GridT);

	// init density kernel and wave functions.
	LOC.allocate_dm_wfc(GridT);

	//======================================
	// do the charge extrapolation before the density matrix is regenerated.
	// mohan add 2011-04-08
	// because once atoms are moving out of this processor,
	// the density matrix will not map the new atomic configuration,
	//======================================
	// THIS IS A BUG, BECAUSE THE INDEX GridT.trace_lo
	// HAS BEEN REGENERATED, SO WE NEED TO
	// REALLOCATE DENSITY MATRIX FIRST, THEN READ IN DENSITY MATRIX,
	// AND USE DENSITY MATRIX TO DO RHO CALCULATION.-- mohan 2013-03-31
	//======================================
	if(pot.extra_pot=="dm" && istep>1)//xiaohui modify 2015-02-01
	{
		for(int is=0; is<NSPIN; is++)
		{
			ZEROS(CHR.rho[is], pw.nrxx);
			stringstream ssd;
			ssd << global_out_dir << "SPIN" << is + 1 << "_DM" ;
			// reading density matrix,
			LOC.read_dm(is, ssd.str() );
		}

		// calculate the charge density
		if(GAMMA_ONLY_LOCAL)
		{
			UHM.GG.cal_rho();
		}
		else
		{
			UHM.GK.calculate_charge();
		}

		// renormalize the charge density
		CHR.renormalize_rho();

		// initialize the potential
		pot.init_pot( istep-1, pw.strucFac );
	}


	// (9) compute S, T, Vnl, Vna matrix.
	UHM.set_lcao_matrices();

    timer::tick("LOOP_elec","before_solver",'D'); 
	return;
}

void LOOP_elec::solver(const int &istep)
{
    TITLE("LOOP_elec","solver"); 
    timer::tick("LOOP_elec","solver",'D'); 

	// Peize Lin add 2014.04.04, update 2021.03.09
	if(vdwd2_para.flag_vdwd2)
	{
		Vdwd2 vdwd2(ucell,vdwd2_para);
		vdwd2.cal_energy();
		en.evdw = vdwd2.get_energy();
	}
	// jiyy add 2019-05-18, update 2021.05.02
	else if(vdwd3_para.flag_vdwd3)
	{
		Vdwd3 vdwd3(ucell,vdwd3_para);
		vdwd3.cal_energy();
		en.evdw = vdwd3.get_energy();
	}

	// self consistent calculations for electronic ground state
	if (CALCULATION=="scf" || CALCULATION=="md"
			|| CALCULATION=="relax" || CALCULATION=="cell-relax") //pengfei 2014-10-13
	{
		//Peize Lin add 2016-12-03
		switch(exx_lcao.info.hybrid_type)
		{
			case Exx_Global::Hybrid_Type::HF:
			case Exx_Global::Hybrid_Type::PBE0:
			case Exx_Global::Hybrid_Type::HSE:
				exx_lcao.cal_exx_ions();
				break;
			case Exx_Global::Hybrid_Type::No:
			case Exx_Global::Hybrid_Type::Generate_Matrix:
				break;
			default:
				throw invalid_argument(TO_STRING(__FILE__)+TO_STRING(__LINE__));
		}

		// No exx
		if( Exx_Global::Hybrid_Type::No==exx_global.info.hybrid_type  )
		{
			ELEC_scf es;
			es.scf(istep-1);
		}
		else if( Exx_Global::Hybrid_Type::Generate_Matrix == exx_global.info.hybrid_type )
		{
			Exx_Opt_Orb exx_opt_orb;
			exx_opt_orb.generate_matrix();
		}
		else    // Peize Lin add 2016-12-03
		{
			ELEC_scf es;
			es.scf(istep-1);
			if( exx_global.info.separate_loop )
			{
				for( size_t hybrid_step=0; hybrid_step!=exx_global.info.hybrid_step; ++hybrid_step )
				{
					exx_global.info.set_xcfunc(xcf);
					exx_lcao.cal_exx_elec();
					
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
				exx_global.info.set_xcfunc(xcf);

				ELEC_scf es;
				es.scf(istep-1);
				
			}
		}
	}
	else if (CALCULATION=="nscf")
	{
		ELEC_nscf::nscf(UHM);
	}
	else if (CALCULATION=="istate")
	{
		IState_Charge ISC;
		ISC.begin();
	}
	else if (CALCULATION=="ienvelope")
	{
		IState_Envelope IEP;
		IEP.begin();
	}
	else
	{
		WARNING_QUIT("LOOP_elec::solver","CALCULATION type not supported");
	}

    timer::tick("LOOP_elec","solver",'D'); 
	return;
}

