#include "LOOP_cell.h"
#include "LOOP_ions.h"

#include "src_lcao/dftu.h"   //Quxin add for DFT+U on 20201029

// delete in near future
#include "src_pw/global.h"

LOOP_cell::LOOP_cell(){}
LOOP_cell::~LOOP_cell(){}

void LOOP_cell::opt_cell(void)
{
	TITLE("LOOP_cell","opt_cell");

    // Initialize the local wave functions.
    // npwx, eigenvalues, and weights
    // npwx may change according to cell change
    // this function belongs to cell LOOP
    wf.allocate_ekb_wg(kv.nks);

    // Initialize the FFT.
    // this function belongs to cell LOOP
    UFFT.allocate();

    // output is ppcell.vloc 3D local pseudopotentials
	// without structure factors
    // this function belongs to cell LOOP
    ppcell.init_vloc(pw.nggm, ppcell.vloc);

    // Initialize the sum of all local potentials.
    // if ion_step==0, read in/initialize the potentials
    // this function belongs to ions LOOP
    int ion_step=0;
    pot.init_pot(ion_step, pw.strucFac);


	// PLEASE simplify the Exx_Global interface
	// mohan add 2021-03-25
	// Peize Lin 2016-12-03
	if (CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax")
	{
		switch(exx_global.info.hybrid_type)
		{
			case Exx_Global::Hybrid_Type::HF:
			case Exx_Global::Hybrid_Type::PBE0:
			case Exx_Global::Hybrid_Type::HSE:
				exx_lcao.init();
				break;
			case Exx_Global::Hybrid_Type::No:
			case Exx_Global::Hybrid_Type::Generate_Matrix:
				break;
			default:
				throw invalid_argument(TO_STRING(__FILE__)+TO_STRING(__LINE__));
		}
	}	

	// PLEASE do not use INPUT global variable
	// mohan add 2021-03-25
	// Quxin added for DFT+U
	if(INPUT.dft_plus_u) 
	{
		dftu.init(ucell, ParaO);
	}

	LOOP_ions ions;
	ions.opt_ions();
	
	return;
}

