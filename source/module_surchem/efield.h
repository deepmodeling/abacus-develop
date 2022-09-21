#ifndef EFIELD_H
#define EFIELD_H

#include "../module_cell/unitcell.h"
#include "../module_pw/pw_basis.h"
#include "surchem.h"

class Efield
{
public:
    Efield();
    ~Efield();

    static ModuleBase::matrix add_efield(const UnitCell &cell, 
                                            ModulePW::PW_Basis *rho_basis, 
                                            const int &nspin, 
                                            const double *const *const rho,
                                            surchem &solvent);

    static double cal_elec_dipole(const UnitCell &cell, 
                                ModulePW::PW_Basis *rho_basis, 
                                const int &nspin, 
                                const double *const *const rho,
                                const double &bmod);

    static double cal_ion_dipole(const UnitCell &cell, const double &bmod);

    static double cal_induced_dipole(const UnitCell &cell, 
                                ModulePW::PW_Basis *rho_basis, 
                                surchem &solvent,
                                const double &bmod);

    static double saw_function(const double &a, const double &b, const double &x);

    static void compute_force(const UnitCell &cell, ModuleBase::matrix &fdip);

    static void prepare(const UnitCell &cell, double &latvec, double &area);


    static double etotefield;           // dipole energy
    static double tot_dipole;           // total dipole
    static int efield_dir;                    // 0, 1, 2 denotes x, y, z direction for dipole correction
    static double efield_pos_max;              // the maximum position of the saw function
    static double efield_pos_dec;               // the decrease region length of the saw function
    static double efield_amp ;                 // field amplitude (in a.u.) (1 a.u. = 51.44 10^10 V/m)
    static double bvec[3];
    static double bmod;
};





#endif