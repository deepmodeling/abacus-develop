#ifndef H_HARTREE_PW_H
#define H_HARTREE_PW_H

#include "tools.h"
#include "unitcell.h"
#include "pw_basis.h"
#include "use_fft.h"

class H_Hartree_pw 
{
	public:

    H_Hartree_pw();
    ~H_Hartree_pw();

	// the Hartree energy
    static double hartree_energy;

	// compute the Hartree energy
    static matrix v_hartree(
		const UnitCell &cell, 
		PW_Basis &pwb, 
		const int &nspin,
		const double*const*const rho);

	private:


};

#endif //Hartree energy
