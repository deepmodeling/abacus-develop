#include "ELEC_cbands_k.h"
#include "LOOP_elec.h"
#include "LCAO_diago.h"
#include "src_pw/global.h"
#include "src_pw/symmetry_rho.h"
#include "LCAO_evolve.h"
#include "dftu.h"

ELEC_cbands_k::ELEC_cbands_k(){};
ELEC_cbands_k::~ELEC_cbands_k(){};


void ELEC_cbands_k::cal_bands(const int &istep, LCAO_Hamilt &uhm)
{
	TITLE("ELEC_cbands_k","cal_bands");
	timer::tick("ELEC_cbands_k","cal_bands",'E');

	int start_spin = -1;
	uhm.GK.reset_spin(start_spin);
	uhm.GK.allocate_pvpR();
						
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
		
		//--------------------------------------------
		//(2) check if we need to calculate 
		// pvpR = < phi0 | v(spin) | phiR> for a new spin.
		//--------------------------------------------
		if(CURRENT_SPIN == uhm.GK.get_spin() )
		{
			//ofs_running << " Same spin, same vlocal integration." << endl;
		}
		else
		{
			//ofs_running << " (spin change)" << endl;
			uhm.GK.reset_spin( CURRENT_SPIN );

			// if you change the place of the following code,
			// rememeber to delete the #include	
			if(VL_IN_H)
			{
				// vlocal = Vh[rho] + Vxc[rho] + Vl(pseudo)
				uhm.GK.cal_vlocal_k(pot.vr_eff1,GridT);
				// added by zhengdy-soc, for non-collinear case
				// integral 4 times, is there any method to simplify?
				if(NSPIN==4)
				{
					for(int is=1;is<4;is++)
					{
						for(int ir=0; ir<pw.nrxx; ir++)
						{
							pot.vr_eff1[ir] = pot.vr_eff( is, ir);
						}
						uhm.GK.cal_vlocal_k(pot.vr_eff1, GridT, is);
					}
				}
			}
		}


		if(!uhm.init_s)
    	{
    	    WARNING_QUIT("Hamilt_Linear::solve_using_cg","Need init S matrix firstly");
    	}

		//--------------------------------------------
		// (3) folding matrix, 
		// and diagonalize the H matrix (T+Vl+Vnl).
		//--------------------------------------------

		// with k points
		timer::tick("Efficience","each_k");
		timer::tick("Efficience","H_k");
		uhm.calculate_Hk(ik);

		// Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
		if(INPUT.dft_plus_u)
		{		
			dftu.cal_eff_pot_mat(ik, istep);

			for(int irc=0; irc<ParaO.nloc; irc++)
			{
				LM.Hloc2[irc] += dftu.pot_eff_k.at(ik).at(irc);
			}							
		}

		timer::tick("Efficience","H_k");

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

		// write the wave functions into LOWF.WFC_K[ik].
		timer::tick("Efficience","diago_k");
		Diago_LCAO_Matrix DLM;
		DLM.solve_complex_matrix(ik, LOWF.WFC_K[ik], LOC.wfc_dm_2d.wfc_k[ik]);
		timer::tick("Efficience","diago_k");

		timer::tick("Efficience","each_k");
	} // end k
			
	// LiuXh modify 2019-07-15*/
	if(!ParaO.out_hsR)
	{
		uhm.GK.destroy_pvpR();
	}

	timer::tick("ELEC_cbands_k","cal_bands",'E');
	return;	
}

