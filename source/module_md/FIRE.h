#ifndef FIRE_H
#define FIRE_H

#include "verlet.h"

class FIRE : public Verlet
{
public:
    FIRE(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in);
    ~FIRE();

    void setup(ModuleESolver::ESolver *p_ensolve);
    void first_half();
    void second_half();
    void outputMD(std::ofstream &ofs);
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