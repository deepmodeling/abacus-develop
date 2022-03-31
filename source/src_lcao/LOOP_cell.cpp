#include "LOOP_cell.h"
#include "LOOP_ions.h"

#include "dftu.h"   //Quxin add for DFT+U on 20201029
#include "dmft.h"

// delete in near future
#include "../src_pw/global.h"

LOOP_cell::LOOP_cell(Parallel_Orbitals &pv)
{
    // * allocate H and S matrices according to computational resources
    // * set the 'trace' between local H/S and global H/S
    this->LM.divide_HS_in_frag(GlobalV::GAMMA_ONLY_LOCAL, pv);
}
LOOP_cell::~LOOP_cell() {}

void LOOP_cell::opt_cell(ORB_control &orb_con, ModuleESolver::ESolver *p_esolver)
{
	ModuleBase::TITLE("LOOP_cell","opt_cell");

    // Initialize the local wave functions.
    // npwx, eigenvalues, and weights
    // npwx may change according to cell change
    // this function belongs to cell LOOP
    GlobalC::wf.allocate_ekb_wg(GlobalC::kv.nks);

    // Initialize the FFT.
    // this function belongs to cell LOOP
    GlobalC::UFFT.allocate();

    // output is GlobalC::ppcell.vloc 3D local pseudopotentials
	// without structure factors
    // this function belongs to cell LOOP
    GlobalC::ppcell.init_vloc(GlobalC::pw.nggm, GlobalC::ppcell.vloc);

    // Initialize the sum of all local potentials.
    // if ion_step==0, read in/initialize the potentials
    // this function belongs to ions LOOP
    int ion_step=0;
    GlobalC::pot.init_pot(ion_step, GlobalC::pw.strucFac);

#ifdef __MPI  
	// PLEASE simplify the Exx_Global interface
	// mohan add 2021-03-25
	// Peize Lin 2016-12-03
	if (GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")
	{
		switch(GlobalC::exx_global.info.hybrid_type)
		{
			case Exx_Global::Hybrid_Type::HF:
			case Exx_Global::Hybrid_Type::PBE0:
			case Exx_Global::Hybrid_Type::HSE:
				GlobalC::exx_lcao.init();
				break;
			case Exx_Global::Hybrid_Type::No:
			case Exx_Global::Hybrid_Type::Generate_Matrix:
				break;
			default:
				throw std::invalid_argument(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
	}	
#endif

	// PLEASE do not use INPUT global variable
	// mohan add 2021-03-25
	// Quxin added for DFT+U
	if(INPUT.dft_plus_u) 
	{
		GlobalC::dftu.init(GlobalC::ucell, this->LM);
	}

  if(INPUT.dft_plus_dmft) GlobalC::dmft.init(INPUT, GlobalC::ucell);

	LOOP_ions ions(this->LM); 
    ions.opt_ions(p_esolver);

	// mohan update 2021-02-10
    orb_con.clear_after_ions(GlobalC::UOT, GlobalC::ORB, GlobalV::deepks_setorb, GlobalC::ucell.infoNL.nproj);
	
	return;
}

