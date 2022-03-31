#ifndef MSST_H
#define MSST_H

#include "verlet.h"

class MSST : public Verlet
{
public:
    MSST(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in);
    ~MSST();

    void setup(ModuleESolver::ESolver *p_ensolve);
    void first_half();
    void second_half();
    void outputMD();
    void write_restart();
    void restart();
    double extra_term();
    double vel_sum();
    void rescale(double volume);
    void propagate_vel();
    void propagate_voldot();

    ModuleBase::Vector3<double> *old_v;
    ModuleBase::Vector3<double> dilation;      // dilation scale
    ModuleBase::Vector3<double> omega;         // time derivative of volume
    double p0;               // initial pressure
    double v0;               // initial volume
    double e0;               // initial energy
    double totmass;          // total mass of the cell
    double lag_pos;          // Lagrangian location of cell
    double vsum;             // sum over v^2

};

#endif