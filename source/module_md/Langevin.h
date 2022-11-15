#ifndef LANGEVIN_H
#define LANGEVIN_H

#include "mdrun.h"

class Langevin : public MDrun
{
public:
    Langevin(MD_parameters& MD_para_in, UnitCell &unit_in);
    ~Langevin();

    void setup(ModuleESolver::ESolver *p_ensolve);
    void first_half();
    void second_half();
    void outputMD(std::ofstream &ofs, bool cal_stress);
    void write_restart();
    void restart();
    void post_force();

    ModuleBase::Vector3<double> *fictitious_force;  // Langevin fictitious_force

};

#endif