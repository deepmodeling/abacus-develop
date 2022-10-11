#ifndef GATEFIELD_H
#define GATEFIELD_H

#include "../module_cell/unitcell.h"
#include "../module_pw/pw_basis.h"

class Gatefield
{
public:
    Gatefield();
    ~Gatefield();

    static void add_gatefield(double *vltot,
                            const UnitCell &cell, 
                            const ModulePW::PW_Basis *rho_basis,
                            const bool &linear, 
                            const bool &quadratic);

    static double mopopla(double &zgate, double z, bool flag);

    static void compute_force(const UnitCell &cell, ModuleBase::matrix &fgate);

    static double etotgatefield;           // energy for gatefield
    static double rho_surface;             // surface charge density of charged plate
    static double zgate;                   // position of charged plate
    static bool relax;                     // 
    static bool block;                     // add a block potential or not
    static double block_down;              // low bound of the block
    static double block_up;                // high bound of the block
    static double block_height;            // height of the block
};

#endif