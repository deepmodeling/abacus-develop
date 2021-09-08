#include "LOOP_cell.h"
#include "LOOP_ions.h"

#include "dftu.h"   //Quxin add for DFT+U on 20201029

// delete in near future
#include "../src_pw/global.h"

LOOP_cell::LOOP_cell(){}
LOOP_cell::~LOOP_cell(){}

void LOOP_cell::opt_cell(void)
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

	// PLEASE do not use INPUT global variable
	// mohan add 2021-03-25
	// Quxin added for DFT+U
	if(INPUT.dft_plus_u) 
	{
		GlobalC::dftu.init(GlobalC::ucell, GlobalC::ParaO);
	}

	LOOP_ions ions;
	ions.opt_ions();

	// mohan update 2021-02-10
    GlobalC::LOWF.orb_con.clear_after_ions(GlobalC::UOT, GlobalC::ORB, INPUT.out_descriptor, GlobalC::ucell.infoNL.nproj);
	
	return;
}

