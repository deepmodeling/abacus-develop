#include "LOOP_cell.h"
#include "LOOP_ions.h"

#include "dftu.h"   //Quxin add for DFT+U on 20201029
#include "dmft.h"

// delete in near future
#include "../src_pw/global.h"

LOOP_cell::LOOP_cell(){}
LOOP_cell::~LOOP_cell() {}

void LOOP_cell::opt_cell(ModuleESolver::ESolver *p_esolver)
{
	ModuleBase::TITLE("LOOP_cell","opt_cell");

	LOOP_ions ions; 
    ions.opt_ions(p_esolver);
	
	return;
}

