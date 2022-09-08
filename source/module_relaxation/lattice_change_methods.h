#ifndef LATTICE_CHANGE_METHODS_H
#define LATTICE_CHANGE_METHODS_H

using namespace std;
#include "lattice_change_basic.h"
#include "lattice_change_cg.h"

class Lattice_Change_Methods
{
	public:

	Lattice_Change_Methods();

	~Lattice_Change_Methods();

	void allocate(void);

	void cal_lattice_change(const int &istep, const ModuleBase::matrix &stress, const double &etot);

	bool get_converged(void)const {return Lattice_Change_Basic::converged;}

	double get_ediff(void)const {return Lattice_Change_Basic::ediff;}

	double get_largest_grad(void)const {return Lattice_Change_Basic::largest_grad;}

	private:

	Lattice_Change_CG lccg;	
};
#endif
