#ifndef VERLET_H
#define VERLET_H

#include "mdrun.h"

class Verlet : public MDrun
{
public:
    Verlet(MD_parameters& MD_para_in, UnitCell &unit_in);
    ~Verlet();

    void setup(ModuleESolver::ESolver *p_ensolve);
    void first_half();
    void second_half();
    void apply_thermostat();
    void thermalize(const int &nraise, const double &current_temp, const double &target_temp);
    void outputMD(std::ofstream &ofs, bool cal_stress);
    void write_restart();
    void restart();

};

#endif