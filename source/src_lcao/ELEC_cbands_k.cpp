#include "ELEC_cbands_k.h"
#include "LCAO_diago.h"
#include "../src_pw/global.h"
#include "../src_pw/symmetry_rho.h"
#include "LCAO_evolve.h"
#include "dftu.h"
#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h"
#endif
#include "../module_base/timer.h"
#include "src_lcao/local_orbital_wfc.h"

ELEC_cbands_k::ELEC_cbands_k(){};
ELEC_cbands_k::~ELEC_cbands_k(){};


void ELEC_cbands_k::cal_bands(const int& istep, LCAO_Hamilt& uhm,
    Local_Orbital_wfc &lowf,
    std::vector<ModuleBase::ComplexMatrix>& dm_k)
{
	ModuleBase::TITLE("ELEC_cbands_k","cal_bands");
    ModuleBase::timer::tick("ELEC_cbands_k", "cal_bands");
    const Parallel_Orbitals* pv = lowf.ParaV;

	int start_spin = -1;
	uhm.GK.reset_spin(start_spin);
	uhm.GK.allocate_pvpR();

#ifdef __DEEPKS
	if (GlobalV::deepks_scf)
    {
		GlobalC::ld.cal_projected_DM_k(dm_k,
			GlobalC::ucell,
            GlobalC::ORB,
            GlobalC::GridD,
            pv->trace_loc_row,
			pv->trace_loc_col,
			GlobalC::kv.nks,
			GlobalC::kv.kvec_d);
    	GlobalC::ld.cal_descriptor();
		//calculate dE/dD
		GlobalC::ld.cal_gedm(GlobalC::ucell.nat);

		//calculate H_V_deltaR from saved <alpha(0)|psi(R)>
		GlobalC::ld.add_v_delta_k(GlobalC::ucell,
            GlobalC::ORB,
            GlobalC::GridD,
            pv->trace_loc_row,
			pv->trace_loc_col,
			pv->nnr);
	}
#endif

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
		
		//--------------------------------------------
		//(2) check if we need to calculate 
		// pvpR = < phi0 | v(spin) | phiR> for a new spin.
		//--------------------------------------------
		if(GlobalV::CURRENT_SPIN == uhm.GK.get_spin() )
		{
			//GlobalV::ofs_running << " Same spin, same vlocal integration." << std::endl;
		}
		else
		{
			//GlobalV::ofs_running << " (spin change)" << std::endl;
			uhm.GK.reset_spin( GlobalV::CURRENT_SPIN );

			// if you change the place of the following code,
			// rememeber to delete the #include	
			if(GlobalV::VL_IN_H)
			{
				// vlocal = Vh[rho] + Vxc[rho] + Vl(pseudo)
				uhm.GK.cal_vlocal_k(GlobalC::pot.vr_eff1,GlobalC::GridT);
				// added by zhengdy-soc, for non-collinear case
				// integral 4 times, is there any method to simplify?
				if(GlobalV::NSPIN==4)
				{
					for(int is=1;is<4;is++)
					{
						for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
						{
							GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff( is, ir);
						}
						uhm.GK.cal_vlocal_k(GlobalC::pot.vr_eff1, GlobalC::GridT, is);
					}
				}
			}
		}


		if(!uhm.init_s)
    	{
    	    ModuleBase::WARNING_QUIT("Hamilt_Linear::solve_using_cg","Need init S matrix firstly");
    	}

		//--------------------------------------------
		// (3) folding matrix, 
		// and diagonalize the H matrix (T+Vl+Vnl).
		//--------------------------------------------

		// with k points
		ModuleBase::timer::tick("Efficience","each_k");
		ModuleBase::timer::tick("Efficience","H_k");
		uhm.calculate_Hk(ik);

		// Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
		if(INPUT.dft_plus_u)
		{
      std::vector<std::complex<double>> eff_pot(lowf.ParaV->nloc);
			GlobalC::dftu.cal_eff_pot_mat_complex(ik, istep, &eff_pot[0]);
      
			for(int irc=0; irc<lowf.ParaV->nloc; irc++)
				uhm.LM->Hloc2[irc] += eff_pot[irc];					
		}

		ModuleBase::timer::tick("Efficience","H_k");

		// Peize Lin add at 2020.04.04
		if(GlobalC::restart.info_load.load_H && !GlobalC::restart.info_load.load_H_finish)
		{
			GlobalC::restart.load_disk(*uhm.LM, "H", ik);
			GlobalC::restart.info_load.load_H_finish = true;
		}
		if(GlobalC::restart.info_save.save_H)
		{
			GlobalC::restart.save_disk(*uhm.LM, "H", ik);
		}

		// write the wave functions into wfc_k_grid[ik].
		ModuleBase::timer::tick("Efficience","diago_k");
		Diago_LCAO_Matrix DLM(uhm.LM);
		DLM.solve_complex_matrix(ik, lowf);
		ModuleBase::timer::tick("Efficience","diago_k");

		ModuleBase::timer::tick("Efficience","each_k");
	} // end k
			
	// LiuXh modify 2019-07-15*/
	if(!Pdiag_Double::out_mat_hsR)
	{
		uhm.GK.destroy_pvpR();
	}

	ModuleBase::timer::tick("ELEC_cbands_k","cal_bands");
	return;	
}

