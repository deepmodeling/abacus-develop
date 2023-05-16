#ifndef H_EWALD_PW_H
#define H_EWALD_PW_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_cell/unitcell.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_hamilt_pw/hamilt_pwdft/forces.h"
#include "module_hamilt_pw/hamilt_pwdft/stress_func.h"

class H_Ewald_pw 
{
  public:
    H_Ewald_pw();
    ~H_Ewald_pw();

    // compute the Ewald energy
    static double compute_ewald(const UnitCell& cell,
                                const ModulePW::PW_Basis* rho_basis,
                                const ModuleBase::ComplexMatrix& strucFac);

  public:
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
