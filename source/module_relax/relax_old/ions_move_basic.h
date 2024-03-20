#ifndef IONS_MOVE_BASIC_H
#define IONS_MOVE_BASIC_H

#include "module_base/matrix.h"
#include "module_cell/unitcell.h"

namespace Ions_Move_Basic
{
extern int dim;             // dimension of the free variables,
extern bool converged;      // converged force or not,
extern double largest_grad; // largest gradient among the forces,
extern int update_iter;     // number of sucesfully updated iterations,
extern int istep;           // index of ionic steps,
extern double ediff;        // energy difference compared to last step,
extern double etot;         // total energy of this step,
extern double etot_p;       // total energy of last step,

extern double trust_radius;     // trust radius now,
extern double trust_radius_old; // old trust radius,
extern double relax_bfgs_rmax;  // max value of trust radius,
extern double relax_bfgs_rmin;  // min value of trust radius,
extern double relax_bfgs_init;  // initial value of trust radius,
extern double best_xxx;         // the last step length of cg , we use it as  bfgs`s initial step length

extern int out_stru; // output the structure or not
// funny way to pass this parameter, but nevertheless

//----------------------------------------------------------------------------
// setup the gradient, all the same for any geometry optimization methods.
//----------------------------------------------------------------------------
void setup_gradient(const UnitCell &ucell, const ModuleBase::matrix &force, double *pos, double *grad);

//----------------------------------------------------------------------------
// move the atom positions, considering the periodic boundary condition.
//----------------------------------------------------------------------------
void move_atoms(UnitCell &ucell, double *move, double *pos);

//----------------------------------------------------------------------------
// check the converged conditions ( if largest gradient is smaller than
// the threshold)
//----------------------------------------------------------------------------
void check_converged(const UnitCell &ucell, const double *grad);

//----------------------------------------------------------------------------
// terminate the geometry optimization.
//----------------------------------------------------------------------------
void terminate(const UnitCell &ucell);

//----------------------------------------------------------------------------
// setup the total energy, keep the new energy or not.
//----------------------------------------------------------------------------
void setup_etot(const double &energy_in, const bool judgement);

double dot_func(const double *a, const double *b, const int &dim);

//----------------------------------------------------------------------------
// third order interpolation scheme,
//----------------------------------------------------------------------------
void third_order();

} // namespace Ions_Move_Basic
#endif
