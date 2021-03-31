#include "ELEC_cbands_gamma.h"
#include "LOOP_elec.h"
#include "LCAO_diago.h"
#include "src_pw/global.h"
#include "src_pw/symmetry_rho.h"
#include "LCAO_evolve.h"
#include "dftu.h"


ELEC_cbands_gamma::ELEC_cbands_gamma(){};
ELEC_cbands_gamma::~ELEC_cbands_gamma(){};


void ELEC_cbands_gamma::cal_bands(const int &istep, LCAO_Hamilt &uhm)
{
	TITLE("ELEC_cbands_gamma","cal_bands");
	timer::tick("ELEC_cband_gamma","cal_bands",'E');

	assert(NSPIN == kv.nks);
						
	// pool parallization in future -- mohan note 2021-02-09
	for(int ik=0; ik<kv.nks; ik++)
	{	
		//-----------------------------------------
		//(1) prepare data for this k point.
		// copy the local potential from array.
		//-----------------------------------------
		if(NSPIN==2) 
		{
			CURRENT_SPIN = kv.isk[ik];
		}
		wf.npw = kv.ngk[ik];

		for(int ir=0; ir<pw.nrxx; ir++)
		{
			pot.vr_eff1[ir] = pot.vr_eff( CURRENT_SPIN, ir);
		}
		
		if(!uhm.init_s)
    	{
    	    WARNING_QUIT("Hamilt_Linear::solve_using_cg","Need init S matrix firstly");
    	}

		//--------------------------------------------
		// (3) folding matrix, 
		// and diagonalize the H matrix (T+Vl+Vnl).
		//--------------------------------------------

		// Peize Lin add ik 2016-12-03
		uhm.calculate_Hgamma(ik);

		// Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
		if(INPUT.dft_plus_u) 
		{
			dftu.cal_eff_pot_mat(ik, istep);

			const int spin = kv.isk[ik];
			for(int irc=0; irc<ParaO.nloc; irc++)
			{
				LM.Hloc[irc] += dftu.pot_eff_gamma.at(spin).at(irc);
			}
		}

		// Peize Lin add at 2020.04.04
		if(restart.info_load.load_H && !restart.info_load.load_H_finish)
		{
			restart.load_disk("H", ik);
			restart.info_load.load_H_finish = true;
		}			
		if(restart.info_save.save_H)
		{
			restart.save_disk("H", ik);
		}

		// SGO: sub_grid_operation
		SGO.cal_totwfc();

		//--------------------------------------
		// DIAG GROUP OPERATION HERE
		//--------------------------------------
		if(DCOLOR==0)
		{
			Diago_LCAO_Matrix DLM;
			// the temperary array totwfc only have one spin direction.
			DLM.solve_double_matrix(ik, SGO.totwfc[0], LOC.wfc_dm_2d.wfc_gamma[ik]);
		}
		else
		{
#ifdef __MPI
			ofs_running << " no diagonalization." << endl;
#else
			cout << " DCOLOR=" << DCOLOR << endl;
			WARNING_QUIT("ELEC_cbands_gamma::cal_bands","no diagonalization");
#endif

		}
#ifdef __MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		// distribute the wave functions again.
		// delete the function -- mohan 2021-02-09
		SGO.dis_subwfc();
	}// end k points
			
	timer::tick("ELEC_cband_gamma","cal_bands",'E');
	return;	
}

