#include "ELEC_evolve.h"
#include "LOOP_elec.h"
#include "LCAO_diago.h"
#include "../src_pw/global.h"
#include "../src_pw/symmetry_rho.h"
#include "LCAO_evolve.h"
#include "dftu.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"

ELEC_evolve::ELEC_evolve(){};
ELEC_evolve::~ELEC_evolve(){};

int ELEC_evolve::tddft;
double ELEC_evolve::td_scf_thr;
double ELEC_evolve::td_dt;
double ELEC_evolve::td_force_dt;
int ELEC_evolve::td_val_elec_01;
int ELEC_evolve::td_val_elec_02;
int ELEC_evolve::td_val_elec_03;
int ELEC_evolve::td_vext;
int ELEC_evolve::td_vext_dire;
double ELEC_evolve::td_timescale;
int ELEC_evolve::td_vexttype;
int ELEC_evolve::td_vextout;
int ELEC_evolve::td_dipoleout;

// this routine only serves for TDDFT using LCAO basis set
void ELEC_evolve::evolve_psi(
	const int &istep,
    LCAO_Hamilt& uhm,
    Local_Orbital_wfc &lowf)
{
	ModuleBase::TITLE("ELEC_evolve","eveolve_psi");
	ModuleBase::timer::tick("ELEC_evolve","evolve_psi");

	int start_spin = -1;
	uhm.GK.reset_spin(start_spin);
	uhm.GK.allocate_pvpR();

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
			uhm.GK.reset_spin( GlobalV::CURRENT_SPIN );

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
						GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(is, ir);
					}
					uhm.GK.cal_vlocal_k(GlobalC::pot.vr_eff1, GlobalC::GridT, is);
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
		uhm.calculate_Hk(ik);

		// Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
		if(INPUT.dft_plus_u)
		{
      std::vector<std::complex<double>> eff_pot(lowf.ParaV->nloc);
			GlobalC::dftu.cal_eff_pot_mat_complex(ik, istep, &eff_pot[0]);

			for(int irc=0; irc<lowf.ParaV->nloc; irc++)
				uhm.LM->Hloc2[irc] += eff_pot[irc];
		}

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
		ModuleBase::timer::tick("Efficience","evolve_k");
		Evolve_LCAO_Matrix ELM(uhm.LM);
		ELM.evolve_complex_matrix(ik, lowf);
		ModuleBase::timer::tick("Efficience","evolve_k");
	} // end k

	// LiuXh modify 2019-07-15*/
	if(!Pdiag_Double::out_mat_hsR)
	{
		uhm.GK.destroy_pvpR();
	}

	ModuleBase::timer::tick("ELEC_evolve","evolve_psi");
	return;
}
