#include "ELEC_cbands_gamma.h"
#include "LOOP_elec.h"
#include "LCAO_diago.h"
#include "../src_pw/global.h"
#include "../src_pw/symmetry_rho.h"
#include "LCAO_evolve.h"
#include "dftu.h"
#include "../module_base/timer.h"


ELEC_cbands_gamma::ELEC_cbands_gamma(){};
ELEC_cbands_gamma::~ELEC_cbands_gamma(){};


void ELEC_cbands_gamma::cal_bands(const int& istep, LCAO_Hamilt& uhm,
    Local_Orbital_wfc &lowf,
    std::vector<ModuleBase::matrix>& dm_gamma)
{
	ModuleBase::TITLE("ELEC_cbands_gamma","cal_bands");
	ModuleBase::timer::tick("ELEC_cband_gamma","cal_bands");

	assert(GlobalV::NSPIN == GlobalC::kv.nks);
						
	// pool parallization in future -- mohan note 2021-02-09
	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{	
		//-----------------------------------------
		//(1) prepare data for this k point.
		// copy the local potential from array.
		//-----------------------------------------
		if(GlobalV::NSPIN==2) 
		{
			GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
		}
		GlobalC::wf.npw = GlobalC::kv.ngk[ik];

		for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
		{
			GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff( GlobalV::CURRENT_SPIN, ir);
		}
		
		if(!uhm.init_s)
    	{
    	    ModuleBase::WARNING_QUIT("Hamilt_Linear::solve_using_cg","Need init S matrix firstly");
    	}

		//--------------------------------------------
		// (3) folding matrix, 
		// and diagonalize the H matrix (T+Vl+Vnl).
		//--------------------------------------------

		// Peize Lin add ik 2016-12-03
		uhm.calculate_Hgamma(ik, dm_gamma);

    // Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
		if(INPUT.dft_plus_u) 
		{
      std::vector<double> eff_pot(lowf.ParaV->nloc);
			GlobalC::dftu.cal_eff_pot_mat_real(ik, istep, &eff_pot[0]);

			const int spin = GlobalC::kv.isk[ik];
			for(int irc=0; irc<lowf.ParaV->nloc; irc++)
				uhm.LM->Hloc[irc] += eff_pot[irc];
        
		}

		// Peize Lin add at 2020.04.04
		if(GlobalC::restart.info_load.load_H && !GlobalC::restart.info_load.load_H_finish)
		{
			GlobalC::restart.load_disk( *uhm.LM, "H", ik);
			GlobalC::restart.info_load.load_H_finish = true;
		}			
		if(GlobalC::restart.info_save.save_H)
		{
			GlobalC::restart.save_disk( *uhm.LM, "H", ik);
		}

		// SGO: sub_grid_operation
		//GlobalC::SGO.cal_totwfc(); //LiuXh modify 2021-09-06, clear memory, totwfc not used now

		//--------------------------------------
		// DIAG GROUP OPERATION HERE
		//--------------------------------------
		if(GlobalV::DCOLOR==0)
		{
			Diago_LCAO_Matrix DLM(uhm.LM);
			// the temperary array totwfc only have one spin direction.
			//DLM.solve_double_matrix(ik, GlobalC::SGO.totwfc[0], wfc_gamma[ik]);
			DLM.solve_double_matrix(ik, lowf); //LiuXh modify 2021-09-06, clear memory, totwfc not used now
		}
		else
		{
#ifdef __MPI
			GlobalV::ofs_running << " no diagonalization." << std::endl;
#else
			std::cout << " DCOLOR=" << GlobalV::DCOLOR << std::endl;
			ModuleBase::WARNING_QUIT("ELEC_cbands_gamma::cal_bands","no diagonalization");
#endif

		}
#ifdef __MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		// distribute the wave functions again.
		// delete the function -- mohan 2021-02-09
		//GlobalC::SGO.dis_subwfc(); //LiuXh modify 2021-09-06, clear memory, totwfc and WFC_GAMMA not used now
	}// end k points
			
	ModuleBase::timer::tick("ELEC_cband_gamma","cal_bands");
	return;	
}

