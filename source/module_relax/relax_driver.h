#ifndef RELAX_DRIVER_H
#define RELAX_DRIVER_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_elecstate/module_charge/charge_extra.h"
#include "relax_old/relax_old.h"
#include "relax_new/relax.h"
#include "module_esolver/esolver.h"

class Relax_Driver
{

	public:

    Relax_Driver(){};
    ~Relax_Driver(){};

    void relax_driver(ModuleESolver::ESolver *p_esolver);

	private:

	// mohan add 2021-01-28
	// mohan moved this variable from electrons.h to relax_driver.h
    int istep;

	//new relaxation method
	Relax rl;

	//old relaxation method
	Relax_old rl_old;

};

#endif
