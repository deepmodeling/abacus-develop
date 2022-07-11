#ifndef VNL_IN_PW_H
#define VNL_IN_PW_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_base/complexarray.h"
#include "../module_base/complexmatrix.h"
#include "VL_in_pw.h"
#ifdef __LCAO
#include "../module_orbital/ORB_gen_tables.h"
#endif
#include "../src_lcao/wavefunc_in_pw.h"
#include "../module_cell/unitcell_pseudo.h"

//==========================================================
// Calculate the non-local pseudopotential in reciprocal
// space using plane wave as basis set.
//==========================================================
class pseudopot_cell_vnl: public pseudopot_cell_vl
{

public:

	pseudopot_cell_vnl();
	~pseudopot_cell_vnl();

	friend class Stress_Func;
	friend class Forces;
	friend class Epsilon0_vasp;
	friend class Potential;
	friend class Hamilt_PW;
	friend class HamiltPW;
	friend class WF_atomic;
	friend class wavefunc;
	friend class Stochastic_hchi;

	void init(const int ntype, const bool allocate_vkb=1);

    double cell_factor; //LiuXh add 20180619

	int nkb; // total number of beta functions considering all atoms

	int lmaxkb; // max angular momentum for non-local projectors

	void init_vnl(UnitCell_pseudo &cell);



	void getvnl(const int &ik, ModuleBase::ComplexMatrix& vkb_in)const;

	void getvnl_alpha(const int &ik);

	void init_vnl_alpha(void);

//===============================================================
// MEMBER VARIABLES :
// NAME : nqx(number of interpolation points)
// NAME : nqxq(size of interpolation table)
// NAME : nhm(max number of different beta functions per atom)
// NAME : lmaxq
// NAME : dq(space between points in the pseudopotential tab)
//===============================================================
// private:
	int calculate_nqx(const double &ecutwfc,const double &dq);

	int nhm;

	int lmaxq;

	ModuleBase::matrix indv;		// indes linking  atomic beta's to beta's in the solid
	ModuleBase::matrix nhtol;      	// correspondence n <-> angular momentum l
	ModuleBase::matrix nhtolm;     	// correspondence n <-> combined lm index for (l,m)
	ModuleBase::matrix nhtoj;		// new added

	ModuleBase::realArray dvan;		//(:,:,:),  the D functions of the solid
	ModuleBase::ComplexArray dvan_so;	//(:,:,:),  spin-orbit case,  added by zhengdy-soc

	ModuleBase::realArray tab;		//(:,:,:), interpolation table for PPs
	ModuleBase::realArray tab_alpha;
	ModuleBase::realArray tab_at;	//(:,:,:), interpolation table for atomic wfc

	ModuleBase::realArray deeq;		//(:,:,:,:), the integral of V_eff and Q_{nm}
#ifdef __CUDA
	double *d_deeq;
#endif
	ModuleBase::ComplexArray deeq_nc;	//(:,:,:,:), the spin-orbit case
	ModuleBase::realArray becsum;	//(:,:,:,:), \sum_i  f(i) <psi(i)/beta_1><beta_m/psi(i)> //used in charge


	mutable ModuleBase::ComplexMatrix vkb;	// all beta functions in reciprocal space
	std::complex<double> ***vkb1_alpha;
	std::complex<double> ***vkb_alpha;

	// other variables
	std::complex<double> Cal_C(int alpha, int lu, int mu, int L, int M);

	double CG(int l1, int m1, int l2, int m2, int L, int M);

	void print_vnl(std::ofstream &ofs);
	#ifdef __LCAO
	ORB_gaunt_table MGT;
	#endif
};

#endif // VNL_IN_PW
