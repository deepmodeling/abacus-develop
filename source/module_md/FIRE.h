#ifndef FIRE_H
#define FIRE_H

#include "mdrun.h"

class FIRE : public MDrun
{
public:
    FIRE(MD_parameters& MD_para_in, UnitCell &unit_in);
    ~FIRE();

    void setup(ModuleESolver::ESolver *p_ensolve);
    void first_half();
    void second_half();
    void outputMD(std::ofstream &ofs, bool cal_stress);
    void write_restart();
    void restart();
    void check_force();
    void check_FIRE();

    double max;            // max force
    double alpha_start;    // alpha_start begin
    double alpha;          // alpha begin
    double finc;           // finc begin
    double fdec;           // fdec begin
    double f_alpha;
    int N_min;
    double dt_max;         // dt_max
    int negative_count;    // Negative count

};

#endif