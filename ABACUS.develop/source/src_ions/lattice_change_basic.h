#ifndef LATTICE_CHANGE_BASIC_H
#define LATTICE_CHANGE_BASIC_H

using namespace std;
#include "../module_base/matrix.h"

namespace Lattice_Change_Basic
{
	extern int dim; // dimension of the free variables,
	extern bool converged; // converged force or not,
	extern double largest_grad; // largest gradient among the forces,
	extern int update_iter; // number of sucesfully updated iterations, 
	extern int istep; // index of ionic steps, 
	extern double ediff; // energy difference compared to last step,
	extern double etot; // total energy of this step,
	extern double etot_p; // total energy of last step,

	extern double lattice_change_ini; // initial value of trust radius,  

	extern int out_stru; // output the structure or not,

	//----------------------------------------------------------------------------
	// setup the gradient, all the same for any geometry optimization methods.
	//----------------------------------------------------------------------------
	void setup_gradient(double* lat, double *grad, matrix &stress);

	//----------------------------------------------------------------------------
	// move the atom positions, considering the periodic boundary condition. 
	//----------------------------------------------------------------------------
	void change_lattice(double *move, double *lat);

	//----------------------------------------------------------------------------
	// check the converged conditions ( if largest gradient is smaller than
	// the threshold)
	//----------------------------------------------------------------------------
	void check_converged(matrix &stress, double *grad);

	//----------------------------------------------------------------------------
	// terminate the geometry optimization.
	//----------------------------------------------------------------------------
	void terminate(void);

	//----------------------------------------------------------------------------
	// setup the total energy, keep the new energy or not.
	//----------------------------------------------------------------------------
	void setup_etot(const double &energy_in, const bool judgement); 
}
#endif
