#ifndef NVE_H
#define NVE_H

#include "verlet.h"

class NVE : public Verlet
{
public:
    NVE(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in);
    ~NVE();

    void setup(ModuleESolver::ESolver *p_ensolve);
    void first_half();
    void second_half();
    void outputMD();
    void write_restart();
    void restart();

};

#endif