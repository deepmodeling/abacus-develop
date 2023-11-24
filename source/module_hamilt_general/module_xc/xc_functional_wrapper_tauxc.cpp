// This file contains wrapper for the mGGA functionals
// it includes 1 subroutine:
// 1. tau_xc

#include "xc_functional.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

//tau_xc and tau_xc_spin: interface for calling xc_mgga_exc_vxc from LIBXC
//XC_POLARIZED, XC_UNPOLARIZED: internal flags used in LIBXC, denote the polarized(nspin=1) or unpolarized(nspin=2) calculations, definition can be found in xc.h from LIBXC
#ifdef USE_LIBXC
void XC_Functional::tau_xc(const double &rho, const double &grho, const double &atau, double &sxc,
          double &v1xc, double &v2xc, double &v3xc)
{
    double s, v1, v2, v3;
	double lapl_rho, vlapl_rho;
    lapl_rho = grho;
    std::vector<xc_func_type> funcs = init_func(XC_UNPOLARIZED);
    
    sxc = 0.0; v1xc = 0.0; v2xc = 0.0; v3xc = 0.0;

    for(xc_func_type &func : funcs)
    {
        xc_mgga_exc_vxc(&func,1,&rho,&grho,&lapl_rho,&atau,&s,&v1,&v2,&vlapl_rho,&v3);
#ifdef __EXX
        if (func.info->number == XC_MGGA_X_SCAN && get_func_type() == 5)
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
    finish_func(funcs);

	return;
}


void XC_Functional::tau_xc_spin(double rhoup, double rhodw, 
        ModuleBase::Vector3<double> gdr1, ModuleBase::Vector3<double> gdr2,
        double tauup, double taudw,
        double &sxc, double &v1xcup, double &v1xcdw, double &v2xcup, double &v2xcdw, double &v2xcud,
        double &v3xcup, double &v3xcdw)
{

	std::vector<xc_func_type> funcs = init_func(XC_POLARIZED);

    double *rho, *grho, *tau, *v1xc, *v2xc, *v3xc, *sgn, s;
    sxc = v1xcup = v1xcdw = 0.0;
    v2xcup = v2xcdw = v2xcud = 0.0;
    v3xcup = v3xcdw = 0.0;
    rho = new double[2];
    grho= new double[3];
    tau = new double[2];
    v1xc= new double[2];
    v2xc= new double[3];
    v3xc= new double[2];
    sgn = new double[2];
    
    rho[0] = rhoup;
    rho[1] = rhodw;
    grho[0] = gdr1.norm2();
    grho[1] = gdr1 * gdr2;
    grho[2] = gdr2.norm2();
    tau[0] = tauup;
    tau[1] = taudw;

    double lapl[2] = {0.0,0.0};
    double vlapl[2];

    const double rho_threshold = 1E-6;
    const double grho_threshold = 1E-10;

    for(xc_func_type &func : funcs)
    {
        if( func.info->family == XC_FAMILY_MGGA || func.info->family == XC_FAMILY_HYB_MGGA)
        {
            sgn[0] = sgn[1] = 1.0;
            // call Libxc function: xc_mgga_exc_vxc
            xc_mgga_exc_vxc( &func, 1, rho, grho, lapl, tau, &s, v1xc, v2xc, vlapl, v3xc);
            if(func.info->kind==XC_CORRELATION)
            {
                if ( rho[0]<rho_threshold || sqrt(std::abs(grho[0]))<grho_threshold )
                    sgn[0] = 0.0;
                if ( rho[1]<rho_threshold || sqrt(std::abs(grho[2]))<grho_threshold )
                    sgn[1] = 0.0;
            }
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

    delete[] grho;
    delete[] rho;
    delete[] tau;
    delete[] v1xc;
    delete[] v2xc;
    delete[] v3xc;
    delete[] sgn;
    finish_func(funcs);

	return;
}

#endif