// This file contains wrapper for the mGGA functionals
// it includes 1 subroutine:
// 1. tau_xc

#ifdef USE_LIBXC

#include "xc_functional_libxc.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include <array>

//tau_xc and tau_xc_spin: interface for calling xc_mgga_exc_vxc from LIBXC
//XC_POLARIZED, XC_UNPOLARIZED: internal flags used in LIBXC, denote the polarized(nspin=1) or unpolarized(nspin=2) calculations, definition can be found in xc.h from LIBXC
void XC_Functional_Libxc::tau_xc(
            const std::vector<int> &func_id,
    const double &rho, const double &grho, const double &atau, double &sxc,
          double &v1xc, double &v2xc, double &v3xc)
{
    double s, v1, v2, v3;
	double lapl_rho, vlapl_rho;
    lapl_rho = grho;
    std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func(func_id, XC_UNPOLARIZED);

    sxc = 0.0; v1xc = 0.0; v2xc = 0.0; v3xc = 0.0;

    for(xc_func_type &func : funcs)
    {
        xc_mgga_exc_vxc(&func,1,&rho,&grho,&lapl_rho,&atau,&s,&v1,&v2,&vlapl_rho,&v3);
#ifdef __EXX
        if (func.info->number == XC_MGGA_X_SCAN && XC_Functional::get_func_type() == 5)
        {
            s *= (1.0 - GlobalC::exx_info.info_global.hybrid_alpha);
            v1 *= (1.0 - GlobalC::exx_info.info_global.hybrid_alpha);
            v2 *= (1.0 - GlobalC::exx_info.info_global.hybrid_alpha);
            v3 *= (1.0 - GlobalC::exx_info.info_global.hybrid_alpha);
        }
#endif
        sxc += s * rho;
        v2xc += v2 * 2.0;
        v1xc += v1;
        v3xc += v3;
    }
    XC_Functional_Libxc::finish_func(funcs);

	return;
}


void XC_Functional_Libxc::tau_xc_spin(
        const std::vector<int> &func_id,
        double rhoup, double rhodw,
        ModuleBase::Vector3<double> gdr1, ModuleBase::Vector3<double> gdr2,
        double tauup, double taudw,
        double &sxc, double &v1xcup, double &v1xcdw, double &v2xcup, double &v2xcdw, double &v2xcud,
        double &v3xcup, double &v3xcdw)
{
    sxc = v1xcup = v1xcdw = 0.0;
    v2xcup = v2xcdw = v2xcud = 0.0;
    v3xcup = v3xcdw = 0.0;

    const std::array<double,2> rho = {rhoup, rhodw};
    const std::array<double,3> grho = {gdr1.norm2(), gdr1 * gdr2, gdr2.norm2()};
    const std::array<double,2> tau = {tauup, taudw};

	std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func(func_id, XC_POLARIZED);

    for(xc_func_type &func : funcs)
    {
        if( func.info->family == XC_FAMILY_MGGA || func.info->family == XC_FAMILY_HYB_MGGA)
        {
            constexpr double rho_threshold = 1E-6;
            constexpr double grho_threshold = 1E-10;
            std::array<double,2> sgn = {1.0, 1.0};
            if(func.info->kind==XC_CORRELATION)
            {
                if ( rho[0]<rho_threshold || sqrt(std::abs(grho[0]))<grho_threshold )
                    sgn[0] = 0.0;
                if ( rho[1]<rho_threshold || sqrt(std::abs(grho[2]))<grho_threshold )
                    sgn[1] = 0.0;
            }

            double s = 0.0;
            std::array<double,2> v1xc, v3xc, lapl, vlapl;
            std::array<double,3> v2xc;
            // call Libxc function: xc_mgga_exc_vxc
            xc_mgga_exc_vxc( &func, 1, rho.data(), grho.data(), lapl.data(), tau.data(), &s, v1xc.data(), v2xc.data(), vlapl.data(), v3xc.data());

            sxc += s * (rho[0] * sgn[0] + rho[1] * sgn[1]);
            v1xcup += v1xc[0] * sgn[0];
            v1xcdw += v1xc[1] * sgn[1];
            v2xcup += 2.0 * v2xc[0] * sgn[0];
            v2xcud += v2xc[1] * sgn[0] * sgn[1];
            v2xcdw += 2.0 * v2xc[2] * sgn[1];
            v3xcup += v3xc[0] * sgn[0];
            v3xcdw += v3xc[1] * sgn[1];
        }
    }

    XC_Functional_Libxc::finish_func(funcs);
}

#endif