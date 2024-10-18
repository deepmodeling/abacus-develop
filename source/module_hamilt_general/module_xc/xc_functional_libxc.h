#ifndef XC_FUNCTIONAL_LIBXC_H
#define XC_FUNCTIONAL_LIBXC_H

#ifdef USE_LIBXC

#include "module_base/matrix.h"
#include "module_base/vector3.h"

#include <xc.h>

#include <tuple>
#include <vector>

	class Charge;

namespace XC_Functional_Libxc
{
//-------------------
//  xc_functional_libxc.cpp
//-------------------

	// sets functional type, which allows combination of LIBXC keyword connected by "+"
	//		for example, "XC_LDA_X+XC_LDA_C_PZ"
	extern std::pair<int,std::vector<int>> set_xc_type_libxc(const std::string xc_func_in);

	// converts func_id into corresponding xc_func_type vector
	extern std::vector<xc_func_type> init_func(const std::vector<int> &func_id, const int xc_polarized);

	extern void finish_func(std::vector<xc_func_type> &funcs);


//-------------------
//  xc_functional_libxc_vxc.cpp
//-------------------

	extern std::tuple<double,double,ModuleBase::matrix> v_xc_libxc(
		const std::vector<int> &func_id,
		const int &nrxx, // number of real-space grid
		const double &omega, // volume of cell
		const double tpiba,
		const Charge* const chr); // charge density

	// for mGGA functional
	extern std::tuple<double,double,ModuleBase::matrix,ModuleBase::matrix> v_xc_meta(
		const std::vector<int> &func_id,
		const int &nrxx, // number of real-space grid
		const double &omega, // volume of cell
		const double tpiba,
		const Charge* const chr);


//-------------------
//  xc_functional_libxc_tools.cpp
//-------------------

	// converting rho (abacus=>libxc)
	extern std::vector<double> convert_rho(
		const int nspin,
		const std::size_t nrxx,
		const Charge* const chr);

	// converting rho (abacus=>libxc)
	extern std::tuple<std::vector<double>, std::vector<double>> convert_rho_amag_nspin4(
		const int nspin,
		const std::size_t nrxx,
		const Charge* const chr);

	// calculating grho
	extern std::vector<std::vector<ModuleBase::Vector3<double>>> cal_gdr(
		const int nspin,
		const std::size_t nrxx,
		const std::vector<double> &rho,
		const double tpiba,
		const Charge* const chr);

	// converting grho (abacus=>libxc)
	extern std::vector<double> convert_sigma(
		const std::vector<std::vector<ModuleBase::Vector3<double>>> &gdr);

	// sgn for threshold mask
	extern std::vector<double> cal_sgn(
		const double rho_threshold,
		const double grho_threshold,
		const xc_func_type &func,
		const int nspin,
		const std::size_t nrxx,
		const std::vector<double> &rho,
		const std::vector<double> &sigma);

	// converting etxc from exc (libxc=>abacus)
	extern double convert_etxc(
		const int nspin,
		const std::size_t nrxx,
		const std::vector<double> &sgn,
		const std::vector<double> &rho,
		std::vector<double> exc);

	// converting vtxc and v from vrho and vsigma (libxc=>abacus)
	extern double convert_vtxc_v(
		const xc_func_type &func,
		const int nspin,
		const std::size_t nrxx,
		const std::vector<double> &sgn,
		const std::vector<double> &rho,
		const std::vector<std::vector<ModuleBase::Vector3<double>>> &gdr,
		const std::vector<double> &vrho,
		const std::vector<double> &vsigma,
		const double tpiba,
		const Charge* const chr,
		ModuleBase::matrix &v);

	// dh for gga v
	extern std::vector<std::vector<double>> cal_dh(
		const int nspin,
		const std::size_t nrxx,
		const std::vector<double> &sgn,
		const std::vector<std::vector<ModuleBase::Vector3<double>>> &gdr,
		const std::vector<double> &vsigma,
		const double tpiba,
		const Charge* const chr);

	// convert v for NSPIN=4
	extern ModuleBase::matrix convert_v_nspin4(
		const std::size_t nrxx,
		const Charge* const chr,
		const std::vector<double> &amag,
		const ModuleBase::matrix &v);


//-------------------
//  xc_functional_libxc_wrapper_xc.cpp
//-------------------

	extern void xc_spin_libxc(
		const std::vector<int> &func_id,
		const double &rhoup, const double &rhodw,
		double &exc, double &vxcup, double &vxcdw);


//-------------------
//  xc_functional_libxc_wrapper_gcxc.cpp
//-------------------

	// the entire GGA functional, for nspin=1 case
	extern void gcxc_libxc(
		const std::vector<int> &func_id,
		const double &rho, const double &grho,
		double &sxc, double &v1xc, double &v2xc);

	// the entire GGA functional, for nspin=2 case
	extern void gcxc_spin_libxc(
		const std::vector<int> &func_id,
		const double rhoup, const double rhodw,
		const ModuleBase::Vector3<double> gdr1, const ModuleBase::Vector3<double> gdr2,
		double &sxc, double &v1xcup, double &v1xcdw, double &v2xcup, double &v2xcdw, double &v2xcud);


//-------------------
//  xc_functional_libxc_wrapper_tauxc.cpp
//-------------------

	// wrapper for the mGGA functionals
	extern void tau_xc(
		const std::vector<int> &func_id,
		const double &rho, const double &grho, const double &atau, double &sxc,
		double &v1xc, double &v2xc, double &v3xc);

	extern void tau_xc_spin(
		const std::vector<int> &func_id,
		double rhoup, double rhodw,
		ModuleBase::Vector3<double> gdr1, ModuleBase::Vector3<double> gdr2,
		double tauup, double taudw,
		double &sxc, double &v1xcup, double &v1xcdw, double &v2xcup, double &v2xcdw, double &v2xcud,
		double &v3xcup, double &v3xcdw);
} // namespace XC_Functional_Libxc

#endif // USE_LIBXC

#endif // XC_FUNCTIONAL_LIBXC_H