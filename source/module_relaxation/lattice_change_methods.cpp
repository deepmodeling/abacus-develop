#include "lattice_change_methods.h"
#include "lattice_change_basic.h"
#include "../src_pw/global.h"

Lattice_Change_Methods::Lattice_Change_Methods(){}
Lattice_Change_Methods::~Lattice_Change_Methods(){}

void Lattice_Change_Methods::allocate()
{
	Lattice_Change_Basic::dim = 9;
	lccg.allocate();
	
	return;
}

void Lattice_Change_Methods::cal_lattice_change(const int &istep, const ModuleBase::matrix &stress, const double &etot)
{
	ModuleBase::TITLE("Lattice_Change_Methods","lattice_change_init");
	Lattice_Change_Basic::istep = istep;
	
	lccg.start(stress,etot);
	
	return;	
}


