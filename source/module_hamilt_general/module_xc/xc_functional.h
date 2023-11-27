//==========================================================
// AUTHOR : mohan
// DATE : 2009-04-04
//==========================================================
#ifndef XC_FUNCTIONAL_H
#define XC_FUNCTIONAL_H

#ifdef USE_LIBXC
#include <xc.h>
#else
#include "xc_funcs.h"
#endif	// ifdef USE_LIBXC

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/vector3.h"
#include "module_base/matrix.h"
#include "exx_info.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_cell/unitcell.h"
class XC_Functional
{
	public:

	XC_Functional();
	~XC_Functional();

//-------------------
// subroutines, grouped according to the file they are in:
//-------------------

//-------------------
//  xc_functional_vxc.cpp
//-------------------

// This file contains interface to the xc_functional class
// it includes 3 subroutines:
// 1. v_xc : which takes rho as input, and [etxc, vtxc, v_xc] as output
// 2. v_xc_libxc : which does the same thing as v_xc, but calling libxc
// NOTE : it is only used for nspin = 1 and 2, the nspin = 4 case is treated in v_xc
// 3. v_xc_meta : which takes rho and tau as input, and v_xc as output

	// compute the exchange-correlation energy 
	// [etxc, vtxc, v] = v_xc(...)
    static std::tuple<double,double,ModuleBase::matrix> v_xc(
		const int &nrxx, // number of real-space grid
		const Charge* const chr,
		const UnitCell *ucell); // charge density

	// using libxc
    static std::tuple<double,double,ModuleBase::matrix> v_xc_libxc(
		const int &nrxx, // number of real-space grid
		const double &omega, // volume of cell
		const double tpiba,
		const Charge* const chr); // charge density

	// for mGGA functional
	static std::tuple<double,double,ModuleBase::matrix,ModuleBase::matrix> v_xc_meta(
		const int &nrxx, // number of real-space grid
		const double &omega, // volume of cell
		const double tpiba,
		const Charge* const chr);

//-------------------
//  xc_functional.cpp
//-------------------

// This file contains subroutines for setting the functional
// it includes 4 subroutines:
// 1. get_func_type : which returns the type of functional (func_type):
//		0 = none; 1 = lda; 2 = gga; 3 = mgga; 4 = hybrid lda/gga; 5 = hybrid mgga
// 2. set_xc_type : sets the value of:
//		func_id, which is the LIBXC id of functional
//		func_type, which is as specified in get_func_type
//		use_libxc, whether to use LIBXC. The rule is to NOT use it for functionals that we already have.
// 3. set_xc_type_libxc : sets functional type, which allows combination of LIBXC keyword connected by "+"
//		for example, "XC_LDA_X+XC_LDA_C_PZ"
// 4. init_func : which converts func_id into corresponding xc_func_type vector

	static int get_func_type();
	static void set_xc_type(const std::string xc_func_in);
	static void get_hybrid_alpha(const double alpha_in);
#ifdef USE_LIBXC
	static void set_xc_type_libxc(const std::string xc_func_in);
	static std::vector<xc_func_type> init_func(const int xc_polarized);
	static void finish_func(std::vector<xc_func_type> &funcs);	
#endif

	private:

	static std::vector<int> func_id; // libxc id of functional
	static int func_type; //0:none, 1:lda, 2:gga, 3:mgga, 4:hybrid lda/gga, 5:hybrid mgga
	static bool use_libxc;

	//exx_hybrid_alpha for mixing exx in hybrid functional:
	static double hybrid_alpha;

//-------------------
//  xc_functional_wrapper_xc.cpp
//-------------------

// This file contains wrapper for the LDA functionals
// it includes 3 subroutines:
// 1. xc, which is the wrapper of LDA part
// (i.e. LDA functional and LDA part of GGA functional)
// 2. xc_spin, which is the spin polarized counterpart of xc
// 3. xc_spin_libxc, which is the wrapper for LDA functional, spin polarized

// NOTE : In our own realization of GGA functional, the LDA part
// and gradient correction are calculated separately.
// The LDA part is provided in xc, while the gradient correction is 
// provided in gradcorr through gcxc/gcx_spin+gcc_spin.
// While in LIBXC, the entire GGA functional is provided.
// As a result, xc/xc_spin and xc_spin_libxc are different for GGA,
// the former gives nonzero result, while the latter returns 0.
// Furthermore, the reason for not having xc_libxc is that something like
// xc_libxc which evaluates functional for individual grid points
// is not efficient. For nspin = 1 and 2, v_xc_libxc evaluates potential
// on the entire grid. I'm having xc_spin_libxc because v_xc_libxc
// does not support nspin = 4.

	public:

	// LDA
	static void xc(const double &rho, double &exc, double &vxc);

	// LSDA
	static void xc_spin(const double &rho, const double &zeta,
			double &exc, double &vxcup, double &vxcdw);
	static void xc_spin_libxc(const double &rhoup, const double &rhodw,
			double &exc, double &vxcup, double &vxcdw);

//-------------------
//  xc_functional_wrapper_gcxc.cpp
//-------------------

// This file contains wrapper for the GGA functionals
// it includes 4 subroutines:
// 1. gcxc, which is the wrapper for gradient correction part
// 2. gcx_spin, spin polarized, exchange only
// 3. gcc_spin, spin polarized, correlation only
// 4. gcxc_libxc, the entire GGA functional, LIBXC, for nspin=1 case
// 5. gcxc_spin_libxc, the entire GGA functional, LIBXC, for nspin=2 case

// The difference between our realization (gcxc/gcx_spin/gcc_spin) and
// LIBXC, and the reason for not having gcxc_libxc is explained
// in the NOTE in the comment for xc_functional_wrapper_wc.cpp part

	// GGA
	static void gcxc(const double &rho, const double &grho,
			double &sxc, double &v1xc, double &v2xc);
	static void gcxc_libxc(const double &rho, const double &grho,
			double &sxc, double &v1xc, double &v2xc);

	// spin polarized GGA
	static void gcx_spin(double rhoup, double rhodw, double grhoup2, double grhodw2,
            double &sx, double &v1xup, double &v1xdw, double &v2xup,
            double &v2xdw);
	static void gcc_spin(double rho, double &zeta, double grho, double &sc,
            double &v1cup, double &v1cdw, double &v2c);

	static void gcxc_spin_libxc(double rhoup, double rhodw, 
		ModuleBase::Vector3<double> gdr1, ModuleBase::Vector3<double> gdr2,
        double &sxc, double &v1xcup, double &v1xcdw, double &v2xcup, double &v2xcdw, double &v2xcud);

//-------------------
//  xc_functional_wrapper_tauxc.cpp
//-------------------

// This file contains wrapper for the mGGA functionals
// it includes 1 subroutine:
// 1. tau_xc
// 2. tau_xc_spin

// NOTE : mGGA is realized through LIBXC

#ifdef USE_LIBXC
	// mGGA
	static void tau_xc(const double &rho, const double &grho, const double &atau, double &sxc,
          double &v1xc, double &v2xc, double &v3xc);

	static void tau_xc_spin(double rhoup, double rhodw, 
			ModuleBase::Vector3<double> gdr1, ModuleBase::Vector3<double> gdr2,
			double tauup, double taudw,
			double &sxc, double &v1xcup, double &v1xcdw, double &v2xcup, double &v2xcdw, double &v2xcud,
			double &v3xcup, double &v3xcdw);
#endif 

//-------------------
//  xc_functional_gradcorr.cpp
//-------------------

// This file contains subroutines realted to gradient calculations
// it contains 5 subroutines:
// 1. gradcorr, which calculates gradient correction
// 2. grad_wfc, which calculates gradient of wavefunction
//		it is used in stress_func_mgga.cpp
// 3. grad_rho, which calculates gradient of density
// 4. grad_dot, which calculates divergence of something
// 5. noncolin_rho, which diagonalizes the spin density matrix
//  and gives the spin up and spin down components of the charge.

    static void gradcorr(double& etxc,
                         double& vtxc,
                         ModuleBase::matrix& v,
                         const Charge* const chr,
                         ModulePW::PW_Basis* rhopw,
                         const UnitCell* ucell,
                         std::vector<double>& stress_gga,
                         const bool is_stress = 0);
    static void grad_wfc(const std::complex<double>* rhog,
                         const int ik,
                         std::complex<double>* grad,
                         const ModulePW::PW_Basis_K* wfc_basis,
                         const double tpiba);
    static void grad_rho(const std::complex<double>* rhog,
                         ModuleBase::Vector3<double>* gdr,
                         const ModulePW::PW_Basis* rho_basis,
                         const double tpiba);
    static void grad_dot(const ModuleBase::Vector3<double>* h,
                         double* dh,
                         ModulePW::PW_Basis* rho_basis,
                         const double tpiba);
    static void noncolin_rho(double* rhoout1,
                             double* rhoout2,
                             double* seg,
                             const double* const* const rho,
                             const int nrxx,
                             const double* ux_,
                             const bool lsign_);

    //-------------------
    //  xc_funct_exch_lda.cpp
    //-------------------

    // This file contains realization of LDA exchange functionals
    // Spin unpolarized ones:
    //  1. slater: ordinary Slater exchange with alpha=2/3
    //  2. slater1: Slater exchange with alpha=1
    //  3. slater_rxc : Slater exchange with alpha=2/3 and Relativistic exchange
    // And their spin polarized counterparts:
    //  1. slater_spin
    //  2. slater1_spin
    //  3. slater_rxc_spin

	// For LDA exchange energy
	static void slater(const double &rs, double &ex, double &vx);
	static void slater1(const double &rs, double &ex, double &vx);
	static void slater_rxc(const double &rs, double &ex, double &vx);

	// For LSDA exchange energy
	static void slater_spin( const double &rho, const double &zeta,
		double &ex, double &vxup, double &vxdw);
	static void slater1_spin( const double &rho, const double &zeta,
		double &ex, double &vxup, double &vxdw);
	static void slater_rxc_spin( const double &rho, const double &z,
		double &ex, double &vxup, double &vxdw);

//-------------------
//  xc_funct_corr_lda.cpp
//-------------------

// This file contains realization of LDA correlation functionals
// Spin unpolarized ones:
//  1. pw : Perdew-Wang LDA correlation
//  2. pz : Perdew-Zunger LDA correlation
//  3. lyp : Lee-Yang-Parr correlation
//  4. vwn : Vosko-Wilk-Nusair LDA correlation
//  5. wigner : Wigner
//  6. hl : Hedin-Lunqvist
//  7. gl : Gunnarson-Lunqvist
// And some of their spin polarized counterparts:
//  1. pw_spin
//  2. pz_spin, which calls pz_polarized

	// For LDA correlation energy
	static void pw(const double &rs, const int &iflag, double &ec, double &vc);
	static void pz(const double &rs, const int &iflag, double &ec, double &vc);
	static void lyp(const double &rs, double &ec, double &vc);
	static void vwn(const double &rs, double &ec, double &vc);
	static void wigner(const double &rs, double &ec, double &vc);
	static void hl(const double &rs, double &ec, double &vc);
	static void gl(const double &rs, double &ec, double &vc);

	// For LSDA correlation energy
	static void pw_spin( const double &rs, const double &zeta,
        double &ec, double &vcup, double &vcdw);
	static void pz_spin( const double &rs, const double &zeta,
    	double &ec, double &vcup, double &vcdw);
	static void pz_polarized( const double &rs, double &ec, double &vc);

//-------------------
//  xc_funct_exch_gga.cpp
//-------------------

// This file contains realizations of gradient correction to exchange part
// Spin unpolarized ones:
//  1. becke88 : Becke88 exchange
//  2. ggax : PW91 exchange
//  3. pbex : PBE exchange (and revPBE)
//  4. optx : OPTX, Handy et al.
//  5. wcx : Wu-Cohen exchange
// And some of their spin polarized counterparts:
//  1. becke88_spin

	static void becke88(const double &rho, const double &grho, double &sx, double &v1x, double &v2x);
	static void ggax(const double &rho, const double &grho, double &sx, double &v1x, double &v2x);
	static void pbex(const double &rho, const double &grho, const int &iflag,
		double &sx, double &v1x, double &v2x);
	static void optx(const double rho, const double grho, double &sx, double &v1x, double &v2x);
	static void wcx(const double &rho,const double &grho, double &sx, double &v1x, double &v2x);

	static void becke88_spin(double rho, double grho, double &sx, double &v1x,
		double &v2x);

//-------------------
//  xc_funct_corr_gga.cpp
//-------------------

// This file contains realizations of gradient correction to correlation part
// Spin unpolarized ones:
//  1. perdew86 : P86
//  2. ggac : PW91
//  3. pbec
//  4. glyp
// And some of their spin polarized counterparts:
//  1. perdew86_spin
//  2. ggac_spin
//  3. pbec_spin

	static void perdew86(const double rho, const double grho, double &sc, double &v1c, double &v2c);
	static void ggac(const double &rho,const double &grho, double &sc, double &v1c, double &v2c);
	static void pbec(const double &rho, const double &grho, const int &flag,
		double &sc, double &v1c, double &v2c);
	static void glyp(const double &rho, const double &grho, double &sc, double &v1c, double &v2c);

	static void perdew86_spin(double rho, double zeta, double grho, double &sc,
		double &v1cup, double &v1cdw, double &v2c);
	//static void ggac_spin(double rho, double zeta, double grho, double &sc,
	//	double &v1cup, double &v1cdw, double &v2c);
	static void pbec_spin(double rho, double zeta, double grho, const int &flag, double &sc,
		double &v1cup, double &v1cdw, double &v2c);

//-------------------
//  xc_funct_hcth.cpp
//-------------------
// This file contains realizations of the HCTH GGA functional
// hcth calls pwcorr

	static void hcth(const double rho, const double grho, double &sx, double &v1x, double &v2x);
	static void pwcorr(const double r, const double c[], double &g, double &dg);

};

#endif //XC_FUNCTION_H
