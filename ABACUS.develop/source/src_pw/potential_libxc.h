//==========================================================
// AUTHOR : Peize Lin
// DATE :   2017-09-14
// UPDATE : 2021-02-28
//==========================================================

#ifdef USE_LIBXC

#ifndef POTENTIAL_LIBXC_H
#define POTENTIAL_LIBXC_H

#include "module_base/matrix.h"
#include "src_global/vector3.h"
#include "src_global/global_variable.h"
#include "src_global/global_function.h"
#include <xc.h>
#include <vector>
#include <tuple>

class Potential_Libxc
{
	public:

	//------------------------------------------------
	// evaluate the exchange-correlation (XC) energy
	// by using the input charge density rho_in and rho_core_in
	//------------------------------------------------
	// [etxc, vtxc, v] = v_xc(...)
	static std::tuple<double,double,matrix> v_xc(
		const double * const * const rho_in,
		const double * const rho_core_in);
		
	private:

	//-------------------------------------------
	// return the type of XC functional by 
	// calling init_func()
	//-------------------------------------------
	static std::vector<XC(func_type)> init_func();

	//------------------------------------------------
	// evaluate three quantities: rho, sigma, and gdr
	// according to the input types of XC functionals
	//------------------------------------------------
	// [rho, sigma, gdr] = cal_input(...)
	static std::tuple< 
		std::vector<double>, 
		std::vector<double>, 
		std::vector<std::vector<Vector3<double>>> > 
	cal_input(
		const std::vector<XC(func_type)> &funcs, 
		const double * const * const rho_in,
		const double * const rho_core_in );

	//----------------------------
	// decide the value of spin
	//----------------------------
	static int nspin0() // may need updates from SOC
	{
		if     (NSPIN==1 || (NSPIN==4 && (!DOMAG && !DOMAG_Z)))		return 1;
		else if(NSPIN==2 || (NSPIN==4 && ( DOMAG ||  DOMAG_Z)))		return 2;
		else throw runtime_error(TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
	}
};

#endif

#endif	// ifdef USE_LIBXC
