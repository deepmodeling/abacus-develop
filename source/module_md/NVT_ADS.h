#ifndef NVT_ADS_H
#define NVT_ADS_H

#include "verlet.h"

class NVT_ADS : public Verlet
{
public:
    NVT_ADS(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in);
    ~NVT_ADS();

    void setup(ModuleESolver::ESolver *p_ensolve);
    void first_half();
    void second_half();
    void outputMD();
    void write_restart();
    void restart();

    double nraise;    // nu*dt, 1/nraise is the collision probability 

};

#endif