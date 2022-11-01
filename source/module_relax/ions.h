#ifndef IONS_H
#define IONS_H

#include "../src_pw/electrons.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../src_pw/charge_extra.h"
#include "relax_old/ions_move_methods.h"
#include "relax_old/lattice_change_methods.h"
#include "relax_new/relax.h"
#include "module_esolver/esolver.h"

class Ions
{

	public:

    Ions(){};
    ~Ions(){};

    void opt_ions(ModuleESolver::ESolver *p_esolver);

	private:

	// mohan add 2021-01-28
    Electrons elec;

	// mohan add 2021-01-28
	// mohan moved this variable from electrons.h to ions.h
    int istep;

	//old relaxation method
	Ions_Move_Methods IMM;
	Lattice_Change_Methods LCM;
	//new relaxation method
	Relax rl;

	//seperate force_stress function first
	bool relaxation(ModuleBase::matrix force,ModuleBase::matrix stress,const int &istep, int &force_step, int &stress_step);
	bool if_do_relax();
	bool if_do_cellrelax();
	bool do_relax(const int& istep, int& jstep, const ModuleBase::matrix& ionic_force, const double& total_energy);
	bool do_cellrelax(const int& istep, const int& stress_step, const ModuleBase::matrix& stress, const double& total_energy);
	void reset_after_relax(const int& istep);
	void reset_after_cellrelax(int& force_step, int& stress_step);

    void update_pot(void);

};

#endif
