#ifndef LCAO_DOMAIN_H
#define LCAO_DOMAIN_H

#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"

namespace LCAO_domain
{

    /**
     * @brief prepare gird integration
    */
	void grid_prepare(
			const Grid_Technique& gt, 
			Gint_Gamma &gint_gamma,
			Gint_k &gint_k,
			const ModulePW::PW_Basis& rhopw, 
			const ModulePW::PW_Basis_Big& bigpw);

    /**
     * @brief set the elements of force-related matrices in LCAO method 
    */
    void set_force (
            const Parallel_Orbitals &pv,
            const int& iw1_all,
            const int& iw2_all,
            const double& vx,
            const double& vy,
            const double& vz,
            const char &dtype,
            double* dsloc_x,
            double* dsloc_y,
            double* dsloc_z,
            double* dhloc_fixed_x,
            double* dhloc_fixed_y,
            double* dhloc_fixed_z);

    /**
     * @brief set the elements of stress-related matrices in LCAO method 
    */
	void set_stress (
			const Parallel_Orbitals &pv,
			const int& iw1_all,
			const int& iw2_all,
			const double& vx,
			const double& vy,
			const double& vz,
			const char &dtype,
			const ModuleBase::Vector3<double> &dtau,
			double* dsloc_11,
			double* dsloc_12,
			double* dsloc_13,
			double* dsloc_22,
			double* dsloc_23,
			double* dsloc_33,
			double* dhloc_fixed_11,
			double* dhloc_fixed_12,
			double* dhloc_fixed_13,
			double* dhloc_fixed_22,
			double* dhloc_fixed_23,
			double* dhloc_fixed_33);

}

#endif
