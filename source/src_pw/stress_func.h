#ifndef STRESS_FUNC_H
#define STRESS_FUNC_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include "../module_base/realarray.h"
#include "../module_base/vector3.h"

//-------------------------------------------------------------------
// mohan reconstruction note: 2021-02-07
// the stress code needs reconstructions (by Daye Zheng)
// 1) add explanations for each function, each variable, for 
// main procedures, make the code readable
// 2) divide the stress class into several files, each file
// deals with only one part of the stress, it is convenient for 
// the next-step reconstruction, for example, we want to make 
// pw as an external variable instead of a global variable
// 3) for PW and LCAO, keeps only one copy of the code, for example
// the ewald term needs only one copy, it will reduce the 
// time to maintain both copies of stress codes.
// 4) remain openning interfaces for others to contribute stress
// codes, for example, molecular dynamics will have an ionic stress
// term, +U? exx? may have other stress terms.
// 5) delete useless comments and tests, if you have a useless code,
// please explicitly explain why you want to keep the test
// 6) format should be beautiful! code should be readable like a 
// note (let readers be comfortable)
//-------------------------------------------------------------------

//----------------------------------------------------------------
// compute the stress terms in terms of the plane wave basis set
// the stress terms include: 
// 1) the stress from the electron kinetic energy 
// 2) the stress from the local pseudopotentials
// 3) the stress from the non-local pseudopotentials
// 4) the stress from the Hartree term
// 5) the stress from the non-linear core correction (if any)
// 6) the strees from the exchange-correlation functional term
// 7) the stress from the ewald term (ion-ion intraction under 
//		periodic boundary conditions). 
// 8) the stress from ionic contributions (for molecular dynamics)
//----------------------------------------------------------------

class Stress_Func
{
	public: 

	Stress_Func(){};
	~Stress_Func(){};

//stress functions
// 1) the stress from the electron kinetic energy
	void stress_kin(ModuleBase::matrix& sigma);  //electron kinetic part in PW basis

// 2) the stress from the Hartree term
	void stress_har(ModuleBase::matrix& sigma, const bool is_pw);  //hartree part in PW or LCAO basis

// 3) the stress from the ewald term (ion-ion intraction under 
//		periodic boundary conditions). 
	void stress_ewa(ModuleBase::matrix& sigma, const bool is_pw);     //ewald part in PW or LCAO basis

// 4) the stress from the local pseudopotentials
	void stress_loc(ModuleBase::matrix& sigma, const bool is_pw);  //local pseudopotential part in PW or LCAO
	
	void dvloc_of_g (const int& msh,
			const double* rab,
			const double* r,
			const double* vloc_at,
			const double& zp,
			double*  dvloc);	//used in local pseudopotential stress

	void dvloc_coul (const double& zp, double* dvloc); //used in local pseudopotential stress

// 5) the stress from the non-linear core correction (if any)
	void stress_cc(ModuleBase::matrix& sigma, const bool is_pw); 			//nonlinear core correction stress in PW or LCAO basis

	void deriv_drhoc (
			const bool &numeric,
			const int mesh,
			const double *r,
			const double *rab,
			const double *rhoc,
			double *drhocg);	//used in nonlinear core correction stress

// 6) the stress from the exchange-correlation functional term
	void stress_gga(ModuleBase::matrix& sigma);			//gga part in both PW and LCAO basis
	void stress_mgga(ModuleBase::matrix& sigma);			//gga part in PW basis

// 7) the stress from the non-local pseudopotentials
	void stress_nl(ModuleBase::matrix& sigma);			//nonlocal part in PW basis


	void get_dvnl1(
			ModuleBase::ComplexMatrix &vkb,
			const int ik,
			const int ipol);	//used in nonlocal part in PW basis
	void dylmr2 (
			const int nylm,
			const int ngy,
			ModuleBase::Vector3<double> *gk,
			ModuleBase::matrix &dylm,
			const int ipol);	//used in get_dvnl1()
	void get_dvnl2(
			ModuleBase::ComplexMatrix &vkb,
			const int ik);		//used in nonlocal part in PW basis
	double Polynomial_Interpolation_nl(
			const ModuleBase::realArray &table,
			const int &dim1,
			const int &dim2,
			const double &table_interval,
			const double &x);	//used in get_dvnl2()

	//functions for stress print
	void print_stress(const std::string &name, const ModuleBase::matrix& f, const bool screen, bool ry)const;

	void printstress_total (const ModuleBase::matrix& scs, bool ry);
	
	static double stress_invalid_threshold_ev;

};

#endif








