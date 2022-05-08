#ifndef LOOP_IONS_H
#define LOOP_IONS_H

#include "../src_ions/ions_move_methods.h"
#include "../src_pw/charge_extra.h"
#include "../src_ions/lattice_change_methods.h"
#include "src_lcao/local_orbital_wfc.h"
#include "module_orbital/ORB_control.h"
#include "src_lcao/LCAO_hamilt.h"
#include "module_esolver/esolver.h"

#include <fstream>

class LOOP_ions
{

	public:

	LOOP_ions();
	~LOOP_ions();

	void opt_ions(ModuleESolver::ESolver *p_esolver); //output for dos

	private:

	Ions_Move_Methods IMM;

    Lattice_Change_Methods LCM;
	
	// PLEASE move 'force_stress()'  function to other places, such as FORCE_STRESS.cpp or
	// you might think to create a new file, it is because 'force_stress' do not
	// belong to 'LOOP_ions', 'GlobalC::pot.init_pot' also do not belong to force_stress()
	// the renew of structure factors, etc. should be ran in other places
	// the 'IMM' and 'LCM' objects should be passed to force_stress() via parameters list
	// mohan note 2021-03-23
	bool force_stress(const int &istep, int &force_step, int &stress_step,  ModuleESolver::ESolver* p_esolver);

	int istep;

	// electron charge density extropolation method
	Charge_Extra CE;

};

#endif
