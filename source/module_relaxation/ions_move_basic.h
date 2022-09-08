#ifndef IONS_MOVE_BASIC_H
#define IONS_MOVE_BASIC_H

using namespace std;
#include "../module_base/matrix.h"

namespace Ions_Move_Basic
{
	extern int dim; // dimension of the free variables,
	extern bool converged; // converged force or not,
	extern double largest_grad; // largest gradient among the forces,
	extern int update_iter; // number of sucesfully updated iterations,
	extern int istep; // index of ionic steps,
	extern double ediff; // energy difference compared to last step,
	extern double etot; // total energy of this step,
	extern double etot_p; // total energy of last step,

	extern double trust_radius; // trust radius now,
	extern double trust_radius_old; // old trust radius,
	extern double relax_bfgs_rmax; // max value of trust radius,
	extern double relax_bfgs_rmin; // min value of trust radius,
	extern double relax_bfgs_init; // initial value of trust radius,
        extern double best_xxx;         // the last step length of cg , we use it as  bfgs`s initial step length

	extern int out_stru; // output the structure or not,

	//----------------------------------------------------------------------------
	// setup the gradient, all the same for any geometry optimization methods.
	//----------------------------------------------------------------------------
	void setup_gradient(double *pos, double* grad, const ModuleBase::matrix &force);

	//----------------------------------------------------------------------------
	// move the atom positions, considering the periodic boundary condition.
	//----------------------------------------------------------------------------
	void move_atoms(double* move, double *pos);

	//----------------------------------------------------------------------------
	// check the converged conditions ( if largest gradient is smaller than
	// the threshold)
	//----------------------------------------------------------------------------
	void check_converged(const double *grad);

	//----------------------------------------------------------------------------
	// terminate the geometry optimization.
	//----------------------------------------------------------------------------
	void terminate(void);

	//----------------------------------------------------------------------------
	// setup the total energy, keep the new energy or not.
	//----------------------------------------------------------------------------
	void setup_etot(const double &energy_in, const bool judgement);

	double dot_func(const double* a, const double* b, const int &dim);


	//----------------------------------------------------------------------------
	// second order interpolation scheme,
	//----------------------------------------------------------------------------
	void second_order(
		const double &e0, // energy of previous step
		const double &e1, // energy at this step
		const double *f0, // force at first step
		const double *x, // movement
		const int &dim,
		double &best_x,
		double &best_e);

	//----------------------------------------------------------------------------
	// third order interpolation scheme,
	//----------------------------------------------------------------------------
	void third_order();

}
#endif
