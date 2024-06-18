#ifndef LCAO_DOMAIN_H
#define LCAO_DOMAIN_H

#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_base/vector3.h"
#include "module_hamilt_lcao/hamilt_lcaodft/force_stress_arrays.h"

namespace LCAO_domain
{

    // can used in gamma algorithm.
	void build_Nonlocal_beta_new(
			LCAO_Matrix &lm,
			double* Hloc,
			const UnitCell& ucell,
			const LCAO_Orbitals& orb,
			const ORB_gen_tables& uot,
			Grid_Driver* GridD);


    void build_Nonlocal_mu_new(
			LCAO_Matrix &lm,
            ForceStressArrays &fsr, // mohan 2024-06-16 
			double* HlocR,
			const bool& calc_deri,
			const UnitCell& ucell,
			const LCAO_Orbitals& orb,
			const ORB_gen_tables& uot,
			Grid_Driver* GridD);

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

    /**
     * @brief set each element without derivatives 
    */
	void single_overlap(
            LCAO_Matrix& lm,
			const LCAO_Orbitals& orb,
			const ORB_gen_tables& uot,
			const Parallel_Orbitals& pv,
			const UnitCell& ucell,
			const int nspin,
			const bool cal_stress,
			const int iw1_all,
            const int iw2_all,
			const int m1, 
			const int m2,
			const char &dtype,
			const int T1,
			const int L1,
			const int N1,
			const int T2,
			const int L2,
			const int N2,
			const ModuleBase::Vector3<double> &dtau,
			const ModuleBase::Vector3<double> &tau1,
			const ModuleBase::Vector3<double> &tau2,
			const int npol,
			const int jj,
			const int jj0,
			const int kk,
			const int kk0,
			int& nnr,  // output value
			int& total_nnr, // output value
			double *olm,    // output value
			double *HSloc); // output value

    /**
     * @brief set each element of T matrices 
    */
	void single_derivative(
            LCAO_Matrix& lm,
            ForceStressArrays& fsr,
			const LCAO_Orbitals& orb,
			const ORB_gen_tables& uot,
			const Parallel_Orbitals& pv,
			const UnitCell& ucell,
			const int nspin,
			const bool cal_stress,
            const int iw1_all,
            const int iw2_all,
			const int m1, 
			const int m2,
			const char &dtype,
			const int T1,
			const int L1,
			const int N1,
			const int T2,
			const int L2,
			const int N2,
			const ModuleBase::Vector3<double> &dtau,
			const ModuleBase::Vector3<double> &tau1,
			const ModuleBase::Vector3<double> &tau2,
			const int npol,
			const int jj,
			const int jj0,
			const int kk,
			const int kk0,
			int& nnr,  // output value
			int& total_nnr, // output value
			double *olm); // output value 

    /**
     * @brief set the elements of S and T matrices 
    */
    void build_ST_new(
			LCAO_Matrix& lm,
			ForceStressArrays& fsr,
			const char& dtype,
			const bool& cal_deri,
			const UnitCell& ucell,
			const LCAO_Orbitals& orb,
			const Parallel_Orbitals& pv,
			const ORB_gen_tables& uot,
			Grid_Driver* GridD,
			double* SHlocR,
			bool cal_syns = false,
			double dmax = 0.0);
}

#endif
