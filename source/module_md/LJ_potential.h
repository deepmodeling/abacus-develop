#ifndef LJ_POTENTIAL_H
#define LJ_POTENTIAL_H

#include "cmd_neighbor.h"
#include "../module_base/vector3.h"
#include "../module_base/matrix.h"
#include "../module_cell/unitcell_pseudo.h"
#include "../module_neighbor/sltk_grid_driver.h"

class LJ_potential
{
public:

    LJ_potential();
    ~LJ_potential();

    static double Lennard_Jones(const UnitCell_pseudo &ucell_c, Grid_Driver &grid_neigh, ModuleBase::Vector3<double> *force, ModuleBase::matrix &stress);
    static double Lennard_Jones(const UnitCell_pseudo &ucell_c, CMD_neighbor &cmd_neigh, ModuleBase::Vector3<double> *force, ModuleBase::matrix &stress);
    static double LJ_energy(const double d);
    static ModuleBase::Vector3<double> LJ_force(const double d, const ModuleBase::Vector3<double> dr);
    static void LJ_virial(double *virial, const ModuleBase::Vector3<double> &force, const ModuleBase::Vector3<double> &dtau);

};

#endif