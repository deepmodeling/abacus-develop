#ifndef LCAO_gen_fixedH_H
#define LCAO_gen_fixedH_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_base/vector3.h"

class LCAO_gen_fixedH
{

  public:
    LCAO_Matrix* LM;

    LCAO_gen_fixedH();
    ~LCAO_gen_fixedH();

	void single_overlap(
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

	void single_derivative(
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

    void build_ST_new(const char& dtype,
                      const bool& cal_deri,
                      const UnitCell& ucell,
                      const LCAO_Orbitals& orb,
                      const Parallel_Orbitals& pv,
                      const ORB_gen_tables& uot,
                      Grid_Driver* GridD,
                      double* SHlocR,
                      bool cal_syns = false,
                      double dmax = 0.0);
	// cal_syns : calculate asynchronous overlap matrix for Hefei-NAMD

};

#endif
