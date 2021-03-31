#include "ELEC_evolve.h"
#include "LOOP_elec.h"
#include "LCAO_diago.h"
#include "src_pw/global.h"
#include "src_pw/symmetry_rho.h"
#include "LCAO_evolve.h"
#include "dftu.h"

ELEC_evolve::ELEC_evolve(){};
ELEC_evolve::~ELEC_evolve(){};

// this routine only serves for TDDFT using LCAO basis set
void ELEC_evolve::evolve_psi(const int &istep, LCAO_Hamilt &uhm, complex<double> ***wfc)
{
	TITLE("ELEC_evolve","eveolve_psi");
	timer::tick("ELEC_evolve","evolve_psi",'E');

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
			uhm.GK.reset_spin( CURRENT_SPIN );

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
						pot.vr_eff1[ir] = pot.vr_eff(is, ir);
					}
					uhm.GK.cal_vlocal_k(pot.vr_eff1, GridT, is);
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

		bool diago = true;
		if (istep >= 1) diago = false;
		if(diago)
		{
			timer::tick("Efficience","diago_k");
			Diago_LCAO_Matrix DLM;
			DLM.solve_complex_matrix(ik, LOWF.WFC_K[ik], LOC.wfc_dm_2d.wfc_k[ik]);
			timer::tick("Efficience","diago_k");
		}
		else
		{
			timer::tick("Efficience","evolve_k");
			Evolve_LCAO_Matrix ELM;
			ELM.evolve_complex_matrix(ik, LOWF.WFC_K[ik], wfc[ik]);
			timer::tick("Efficience","evolve_k");
		}
	} // end k
			
	// LiuXh modify 2019-07-15*/
	if(!ParaO.out_hsR)
	{
		uhm.GK.destroy_pvpR();
	}

	timer::tick("ELEC_evolve","evolve_psi",'E');
	return;	
}

