#ifndef LANGEVIN_H
#define LANGEVIN_H

#include "verlet.h"

class Langevin : public Verlet
{
public:
    Langevin(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in);
    ~Langevin();

    void setup(ModuleESolver::ESolver *p_ensolve);
    void first_half();
    void second_half();
    void outputMD();
    void write_restart();
    void restart();
    void post_force();
    void temp_target();

    double t_target;       // target temperature

};

#endif