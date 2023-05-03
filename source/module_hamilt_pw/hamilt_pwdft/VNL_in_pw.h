#ifndef VNL_IN_PW_H
#define VNL_IN_PW_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/complexarray.h"
#include "module_base/complexmatrix.h"
#include "VL_in_pw.h"
#ifdef __LCAO
#include "module_basis/module_ao/ORB_gen_tables.h"
#endif
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_hamilt_lcao/hamilt_lcaodft/wavefunc_in_pw.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_pw/hamilt_pwdft/forces.h"
#include "module_hamilt_pw/hamilt_pwdft/stress_func.h"
#include "module_psi/psi.h"

//==========================================================
// Calculate the non-local pseudopotential in reciprocal
// space using plane wave as basis set.
//==========================================================
class pseudopot_cell_vnl: public pseudopot_cell_vl
{

public:

	pseudopot_cell_vnl();
	~pseudopot_cell_vnl();

	friend class Stress_Func<double>;
	friend class Forces<double>;
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

	void init_vnl(UnitCell &cell, ModulePW::PW_Basis_K* wfc_basis);

    template <typename FPTYPE, typename Device>
    void getvnl(Device * ctx, const int &ik, ModulePW::PW_Basis_K *wfc_basis, std::complex<FPTYPE>* vkb_in)const;

    void getvnl(const int &ik, ModulePW::PW_Basis_K *wfc_basis, ModuleBase::ComplexMatrix& vkb_in)const;

	// void getvnl_alpha(const int &ik, ModulePW::PW_Basis_K *wfc_basis);

	void init_vnl_alpha(void);

	void initgradq_vnl(const UnitCell &cell);

	void getgradq_vnl(const int ik);

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

	ModuleBase::realArray tab;		//(:,:,:), interpolation table for PPs: 4pi/sqrt(V) * \int betar(r)jl(qr)rdr
	ModuleBase::realArray tab_alpha;
	ModuleBase::realArray tab_at;	//(:,:,:), interpolation table for atomic wfc
	ModuleBase::realArray tab_dq;   //4pi/sqrt(V) * \int betar(r)*djl(qr)/d(qr)*r^2 dr

	ModuleBase::realArray deeq;		//(:,:,:,:), the integral of V_eff and Q_{nm}
	bool multi_proj = false;
    float *s_deeq = nullptr;
    double *d_deeq = nullptr;
	ModuleBase::ComplexArray deeq_nc;	//(:,:,:,:), the spin-orbit case
    std::complex<float> *c_deeq_nc = nullptr; // GPU array of deeq_nc
    std::complex<double> *z_deeq_nc = nullptr; // GPU array of deeq_nc
	ModuleBase::realArray becsum;	//(:,:,:,:), \sum_i  f(i) <psi(i)/beta_1><beta_m/psi(i)> //used in charge


	mutable ModuleBase::ComplexMatrix vkb;	// all beta functions in reciprocal space
	mutable ModuleBase::ComplexArray gradvkb; // gradient of beta functions
	std::complex<double> ***vkb1_alpha;
	std::complex<double> ***vkb_alpha;

	// other variables
	std::complex<double> Cal_C(int alpha, int lu, int mu, int L, int M);

	double CG(int l1, int m1, int l2, int m2, int L, int M);

	void print_vnl(std::ofstream &ofs);

	//calculate the effective coefficient matrix for non-local pseudopotential projectors
	void cal_effective_D();
	#ifdef __LCAO
	ORB_gaunt_table MGT;
    #endif

    template <typename FPTYPE> FPTYPE * get_nhtol_data() const;
    template <typename FPTYPE> FPTYPE * get_nhtolm_data() const;
    template <typename FPTYPE> FPTYPE * get_indv_data() const;
    template <typename FPTYPE> FPTYPE * get_tab_data() const;
    template <typename FPTYPE> FPTYPE * get_deeq_data() const;
    template <typename FPTYPE> std::complex<FPTYPE> * get_vkb_data() const;
    template <typename FPTYPE> std::complex<FPTYPE> * get_deeq_nc_data() const;

private:
    float * s_nhtol = nullptr, * s_nhtolm = nullptr, * s_indv = nullptr, * s_tab = nullptr;
    std::complex<float> * c_vkb = nullptr;

    double * d_nhtol = nullptr, * d_nhtolm = nullptr, * d_indv = nullptr, * d_tab = nullptr;
    std::complex<double> * z_vkb = nullptr;

	ModulePW::PW_Basis_K* wfcpw = nullptr;
};

#endif // VNL_IN_PW
