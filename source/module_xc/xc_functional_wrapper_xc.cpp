// This file contains wrapper for the LDA functionals
// it includes 3 subroutines:
// 1. xc, which is the wrapper of LDA part
// (i.e. LDA functional and LDA part of GGA functional)
// 2. xc_spin, which is the spin polarized counterpart of xc
// 3. xc_spin_libxc, which is the wrapper for LDA functional, spin polarized

#include "xc_functional.h"
#include <stdexcept>

void XC_Functional::xc(const double &rho, double &exc, double &vxc)
{

	double third = 1.0 / 3.0;
	double pi34 = 0.6203504908994e0 ; // pi34=(3/4pi)^(1/3)
	double rs;
    double e,v;
    
    exc = vxc = 0.00;
	
	rs = pi34 / std::pow(rho, third);

    for(int id : func_id)
    {
        switch( id )
        {
            // Exchange functionals containing slater exchange
            case XC_LDA_X: case XC_GGA_X_PBE: case XC_GGA_X_PBE_R: 
            case XC_GGA_X_WC: case XC_GGA_X_B88: case XC_GGA_X_PW91:
            //  SLA,PBX,rPBX,PBXsol,WC,B88,PW91_X
                XC_Functional::slater(rs, e, v);break;

            // Exchange functionals containing attenuated slater exchange
            case XC_HYB_GGA_XC_PBEH:
            //  PBE0
                XC_Functional::slater(rs, e, v);
                e *= 0.75; v*= 0.75;
                break;

            // Correlation functionals containing PW correlation
            case XC_GGA_C_PBE: case XC_GGA_C_PW91: case XC_LDA_C_PW:
            //   PBC,PW91,PWLDA
                XC_Functional::pw(rs, 0, e, v);break;

            // Correlation functionals containing PZ correlation
            case XC_LDA_C_PZ: case XC_GGA_C_P86:
            //  PZ,P86
                XC_Functional::pz(rs, 0, e, v);break;

            // Correlation functionals containing LYP correlation
            case XC_GGA_C_LYP:
            //  BLYP
                XC_Functional::lyp(rs, e, v);break;
                            
            default:
                e = v = 0.0;
        }
        exc += e;
        vxc += v;
    }
	return;
}

void XC_Functional::xc_spin(const double &rho, const double &zeta,
		double &exc, double &vxcup, double &vxcdw)
{
	static const double small = 1.e-10;
    double e, vup, vdw;
    exc = vxcup = vxcdw = 0.0;

	static const double third = 1.0 / 3.0;
	static const double pi34 = 0.62035049089940; 
	const double rs = pi34 / pow(rho, third);//wigner_sitz_radius;

    for(int id : func_id)
    {
        switch( id )
        {
            // Exchange functionals containing slater exchange
            case XC_LDA_X: case XC_GGA_X_PBE: case XC_GGA_X_PBE_R: 
            case XC_GGA_X_WC: case XC_GGA_X_B88: case XC_GGA_X_PW91:
            //  SLA,PBX,rPBX,PBXsol,WC,B88,PW91_X
                XC_Functional::slater_spin(rho, zeta, e, vup, vdw);	break;
            
            // Exchange functionals containing attenuated slater exchange
            case XC_HYB_GGA_XC_PBEH:
            //  PBE0
                XC_Functional::slater_spin(rho, zeta, e, vup, vdw);
                e *= 0.75; vup *= 0.75; vdw *=0.75;
                break;

            // Correlation functionals containing PZ correlation
            case XC_LDA_C_PZ: case XC_GGA_C_P86:
            //  PZ,P86
                XC_Functional::pz_spin(rs, zeta, e, vup, vdw); break;

            // Correlation functionals containing PW correlationtests/integrate/101_PW_OU_pseudopot
            case XC_GGA_C_PBE:
            //   PBC,PBCsol
                XC_Functional::pw_spin(rs, zeta, e, vup, vdw); break;

            // Correlation functionals that already contains LDA components
            // so sets to 0 here
            case XC_GGA_X_HCTH_A: case XC_GGA_C_HCTH_A: case XC_GGA_X_OPTX:
            //HCTH_X  HTCH_C      OPTX
                e = vup = vdw =0.0;

            // Cases that are only realized in LIBXC
            default:
                throw std::domain_error("functional unfinished in "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));	break;

        }
        exc += e;
        vxcup += vup;
        vxcdw += vdw;
	}
	return;
}

void XC_Functional::xc_spin_libxc(const double &rhoup, const double &rhodw,
		double &exc, double &vxcup, double &vxcdw)
{
#ifdef USE_LIBXC
    double e, vup, vdw;
    double *rho_ud, *vxc_ud;
    exc = vxcup = vxcdw = 0.0;

    rho_ud = new double[2];
    vxc_ud = new double[2];
    rho_ud[0] = rhoup;
    rho_ud[1] = rhodw;

    std::vector<xc_func_type> funcs = init_func(XC_POLARIZED);

    for(xc_func_type &func : funcs)
    {
        if( func.info->family == XC_FAMILY_LDA)
        {
            // call Libxc function: xc_lda_exc_vxc
            xc_lda_exc_vxc( &func, 1, rho_ud, &e, vxc_ud);
        }
        exc += e;
        vxcup += vxc_ud[0];
        vxcdw += vxc_ud[1];
    }    


    delete[] rho_ud;
    delete[] vxc_ud;
#else
    ModuleBase::WARNING_QUIT("xc_spin_libxc","compile with LIBXC to use this subroutine");
#endif 
}