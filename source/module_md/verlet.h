#ifndef VERLET_H
#define VERLET_H

#include "MD_parameters.h"
#include "../module_cell/unitcell_pseudo.h"
#include "../module_base/matrix.h"
#include "module_esolver/esolver.h"

class Verlet
{
public:
    Verlet(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in);
    virtual ~Verlet();

    virtual void setup(ModuleESolver::ESolver *p_esolver);
    virtual void first_half();
    virtual void second_half();
    virtual void outputMD(std::ofstream &ofs);
    virtual void write_restart();
    virtual void restart();

    MD_parameters &mdp;
	UnitCell_pseudo &ucell;
    bool stop;                           // MD stop or not

    // All parameters are in a.u. unit.
    double temperature_;
    int step_;
    int step_rst_;
    double energy_;
	int frozen_freedom_;

    double *allmass;                     // atom mass 
    ModuleBase::Vector3<double> *pos;    // atom position
    ModuleBase::Vector3<double> *vel;    // atom velocity
    ModuleBase::Vector3<int> *ionmbl;    // atom is frozen or not
    ModuleBase::Vector3<double> *force;  // force of each atom
    ModuleBase::matrix virial;           // virial for this lattice
	ModuleBase::matrix stress;           // stress for this lattice
    double potential;                    // potential energy
    double kinetic;                      // kinetic energy

};

#endif