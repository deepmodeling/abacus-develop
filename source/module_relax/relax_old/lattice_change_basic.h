#ifndef LATTICE_CHANGE_BASIC_H
#define LATTICE_CHANGE_BASIC_H

#include "module_base/matrix.h"
#include "module_cell/unitcell.h"

namespace Lattice_Change_Basic
{
extern int dim;             // dimension of the free variables,
extern bool converged;      // converged force or not,
extern double largest_grad; // largest gradient among the forces,
extern int update_iter;     // number of sucesfully updated iterations,
extern int istep;           // index of ionic steps,
extern int stress_step;     // index of stress step
extern double ediff;        // energy difference compared to last step,
extern double etot;         // total energy of this step,
extern double etot_p;       // total energy of last step,

extern double lattice_change_ini; // initial value of trust radius,
extern std::string fixed_axes;    // convert from INPUT.fixed_axes

//----------------------------------------------------------------------------
// setup the gradient, all the same for any geometry optimization methods.
//----------------------------------------------------------------------------
void setup_gradient(const UnitCell &ucell, double *lat, double *grad, ModuleBase::matrix &stress);

//----------------------------------------------------------------------------
// move the atom positions, considering the periodic boundary condition.
//----------------------------------------------------------------------------
void change_lattice(UnitCell &ucell, double *move, double *lat);

//----------------------------------------------------------------------------
// check the converged conditions ( if largest gradient is smaller than
// the threshold)
//----------------------------------------------------------------------------
void check_converged(const UnitCell &ucell, ModuleBase::matrix &stress, double *grad);

//----------------------------------------------------------------------------
// terminate the geometry optimization.
//----------------------------------------------------------------------------
void terminate(void);

//----------------------------------------------------------------------------
// setup the total energy, keep the new energy or not.
//----------------------------------------------------------------------------
void setup_etot(const double &energy_in, const bool judgement);
} // namespace Lattice_Change_Basic
#endif
