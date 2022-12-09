#ifndef H_EWALD_PW_H
#define H_EWALD_PW_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_cell/unitcell.h"
#include "../module_pw/pw_basis.h"
#include "src_pw/forces.h"
#include "src_pw/stress_func.h"

class H_Ewald_pw 
{
	public:

	// need to be updated in near future 2021-02-25
    friend class Stress_Func<double, psi::DEVICE_CPU>; // Ewald stress
    friend class Stress_Func<double, psi::DEVICE_GPU>; // Ewald stress
    friend class Forces<double, psi::DEVICE_CPU>; // Ewald forces
    friend class Forces<double, psi::DEVICE_GPU>; // Ewald forces
	friend class Force_LCAO; // Ewald forces

    H_Ewald_pw();
    ~H_Ewald_pw();

	// the Ewald energy
    static double ewald_energy;

	// compute the Ewald energy
    static void compute_ewald(const UnitCell &cell, ModulePW::PW_Basis* rho_basis);

	private:

    static void rgen(
        const ModuleBase::Vector3<double> &dtau,
        const double &rmax,
        int *irr,
        const ModuleBase::Matrix3 &at,
        const ModuleBase::Matrix3 &bg,
        ModuleBase::Vector3<double> *r,
        double *r2,
        int  &nrm
    );

	// the coefficient of ewald method
	static double alpha;
    static int mxr;

};

#endif //ewald energy
