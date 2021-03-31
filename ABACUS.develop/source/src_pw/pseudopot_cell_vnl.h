//==========================================================
// AUTHOR : Lixin He,mohan
// DATE : 2008-11-08
//==========================================================
#ifndef PSEUDOPOT_CELL_VNL_H
#define PSEUDOPOT_CELL_VNL_H

#include "tools.h"
#include "pseudopot_cell_vl.h"
#include "../src_lcao/ORB_gen_tables.h"
#include "wavefunc_in_pw.h"
#include "unitcell_pseudo.h"

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
	friend class WF_atomic;
	friend class wavefunc;
	friend class Stochastic_hchi;

	void init(const int ntype, const bool allocate_vkb=1);

    double cell_factor; //LiuXh add 20180619

// NAME : nkb(total number of beta functions, with struct.fact.)
	int nkb; // be called in hm.hpw.init

// NAME : lmaxkb(max angular momentum,(see pseudo_h))
	int lmaxkb;

	void init_vnl(UnitCell_pseudo &cell);

private:

	void getvnl(const int &ik);

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

	int calculate_nqx(const double &ecutwfc,const double &dq);

	int nhm;
	int lmaxq;

	matrix indv;		// indes linking  atomic beta's to beta's in the solid
	matrix nhtol;      	// correspondence n <-> angular momentum l
	matrix nhtolm;     	// correspondence n <-> combined lm index for (l,m)
	matrix nhtoj;		// new added

	realArray dvan;		//(:,:,:),  the D functions of the solid
	ComplexArray dvan_so;	//(:,:,:),  spin-orbit case,  added by zhengdy-soc

	realArray tab;		//(:,:,:), interpolation table for PPs
	realArray tab_alpha;
	realArray tab_at;	//(:,:,:), interpolation table for atomic wfc

	realArray deeq;		//(:,:,:,:), the integral of V_eff and Q_{nm}
	ComplexArray deeq_nc;	//(:,:,:,:), the spin-orbit case
	realArray becsum;	//(:,:,:,:), \sum_i  f(i) <psi(i)/beta_1><beta_m/psi(i)> //used in charge

	ComplexMatrix vkb;	// all beta functions in reciprocal space
	complex<double> ***vkb1_alpha;
	complex<double> ***vkb_alpha;
	

	complex<double> Cal_C(int alpha, int lu, int mu, int L, int M);
	double CG(int l1, int m1, int l2, int m2, int L, int M);
	void print_vnl(ofstream &ofs);
	
	ORB_gaunt_table MGT;
};

#endif // PSEUDOPOT_CELL_VNL_H
