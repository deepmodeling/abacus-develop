#ifndef RELAX_OLD_H
#define RELAX_OLD_H

#include "ions_move_methods.h"
#include "lattice_change_methods.h"

class Relax_old
{
    public:

    void init_relax();
    bool relax_step(ModuleBase::matrix force,ModuleBase::matrix stress,const int &istep, int &force_step, int &stress_step);

    private:

	Ions_Move_Methods IMM;
	Lattice_Change_Methods LCM;

	//seperate force_stress function first
	bool if_do_relax();
	bool if_do_cellrelax();
	bool do_relax(const int& istep, int& jstep, const ModuleBase::matrix& ionic_force, const double& total_energy);
	bool do_cellrelax(const int& istep, const int& stress_step, const ModuleBase::matrix& stress, const double& total_energy);
	void reset_after_relax(const int& istep);
	void reset_after_cellrelax(int& force_step, int& stress_step);
};


#endif