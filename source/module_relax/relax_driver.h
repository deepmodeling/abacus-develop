#ifndef RELAX_DRIVER_H
#define RELAX_DRIVER_H

#include "module_esolver/esolver.h"
#include "relax_new/relax.h"
#include "relax_old/relax_old.h"

class Relax_Driver
{

  public:
    Relax_Driver(){};
    ~Relax_Driver(){};

    void relax_driver(ModuleESolver::ESolver *p_esolver);

  private:
    // mohan add 2021-01-28
    // mohan moved this variable from electrons.h to relax_driver.h
    int istep = 0;
    double etot = 0;

    // new relaxation method
    Relax rl;

    // old relaxation method
    Relax_old rl_old;
};

#endif
