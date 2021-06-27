#ifndef H_EWALD_PW_H
#define H_EWALD_PW_H

#include "tools.h"
#include "module_cell/unitcell.h"
#include "pw_basis.h"

class H_Ewald_pw 
{
	public:

	// need to be updated in near future 2021-02-25
	friend class Stress_Func; // Ewald stress
	friend class Forces; // Ewald forces
	friend class Force_LCAO; // Ewald forces

    H_Ewald_pw();
    ~H_Ewald_pw();

	// the Ewald energy
    static double ewald_energy;

	// compute the Ewald energy
    static void compute_ewald(const UnitCell &cell, const PW_Basis &pwb);

	private:

    static void rgen(
        const Vector3<double> &dtau,
        const double &rmax,
        int *irr,
        const Matrix3 &at,
        const Matrix3 &bg,
        Vector3<double> *r,
        double *r2,
        int  &nrm
    );

	// the coefficient of ewald method
	static double alpha;
    static int mxr;

};

#endif //ewald energy
