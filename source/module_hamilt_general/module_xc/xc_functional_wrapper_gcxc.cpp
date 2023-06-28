// This file contains wrapper for the GGA functionals
// it includes 4 subroutines:
// 1. gcxc, which is the wrapper for gradient correction part
// 2. gcx_spin, spin polarized, exchange only
// 3. gcc_spin, spin polarized, correlation only
// 4. gcxc_libxc, the entire GGA functional, LIBXC, for nspin=1 case
// 5. gcxc_spin_libxc, the entire GGA functional, LIBXC, for nspin=2 case

#include "xc_functional.h"
#include <stdexcept>
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/global_function.h"

void XC_Functional::gcxc(const double &rho, const double &grho, double &sxc,
          double &v1xc, double &v2xc)
{
    //-----------------------------------------------------------------------
    //     gradient corrections for exchange and correlation - Hartree a.u.
    //     exchange  :  Becke88
    //                  GGA (Generalized Gradient Approximation), PW91
    //                  PBE
    //                  revPBE
    //     correlation: Perdew86
    //                  GGA (PW91)
    //                  Lee-Yang-Parr
    //                  PBE
    //
    //     input:  rho, grho=|\nabla rho|^2
    //     definition:  E_x = \int E_x(rho,grho) dr
    //     output: sx = E_x(rho,grho)
    //             v1x= D(E_x)/D(rho)
    //             v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
    //             sc, v1c, v2c as above for correlation
    //
    // use funct
    // USE kinds
    // implicit none
    // real rho, grho, sx, sc, v1x, v2x, v1c, v2c;
    const double small = 1.e-6;
    const double smallg = 1.e-10;
    double s,v1,v2;
    sxc = v1xc = v2xc = 0.0;

    if (rho <= small || grho < smallg)
    {
        return;
    }

    for(int id : func_id)
    {
        switch( id )
        {
            case XC_GGA_X_B88: //B88
                XC_Functional::becke88(rho, grho, s, v1, v2);break;
            case XC_GGA_X_PW91: //PW91_X
                XC_Functional::ggax(rho, grho, s, v1, v2);break;
            case XC_GGA_X_PBE: //PBX
                XC_Functional::pbex(rho, grho, 0, s, v1, v2);break;
            case XC_GGA_X_PBE_R: //revised PBX
                XC_Functional::pbex(rho, grho, 1, s, v1, v2);break;
            case XC_GGA_X_HCTH_A: //HCTH_X
                XC_Functional::hcth(rho, grho, s, v1, v2);break; //XC together
            case XC_GGA_C_HCTH_A: //HCTH_C
                s = 0.0; v1 = 0.0; v2 = 0.0;break;
            case XC_GGA_X_OPTX: //OPTX
                XC_Functional::optx(rho, grho, s, v1, v2);break;
            case XC_GGA_X_PBE_SOL: //PBXsol
                XC_Functional::pbex(rho, grho, 2, s, v1, v2);break;
            case XC_GGA_X_WC: //Wu-Cohen
                XC_Functional::wcx (rho, grho, s, v1, v2);break;
            case XC_GGA_C_P86: //P86
                XC_Functional::perdew86(rho, grho, s, v1, v2);break;
            case XC_GGA_C_PW91: //PW91_C
                XC_Functional::ggac(rho, grho, s, v1, v2);break;
            case XC_GGA_C_PBE: //PBC
                XC_Functional::pbec(rho, grho, 0, s, v1, v2);break;
            case XC_GGA_C_PBE_SOL: //PBCsol
                XC_Functional::pbec(rho, grho, 1, s, v1, v2);break;
            case XC_GGA_C_LYP: //BLYP
                XC_Functional::glyp(rho, grho, s, v1, v2); break;
            case XC_HYB_GGA_XC_PBEH: //PBE0
                double sx, v1x, v2x, sc, v1c, v2c;
                XC_Functional::pbex(rho, grho, 0, sx, v1x, v2x);
                sx *= (1.0 - XC_Functional::hybrid_alpha); 
                v1x *= (1.0 - XC_Functional::hybrid_alpha); 
                v2x *= (1.0 - XC_Functional::hybrid_alpha);
                XC_Functional::pbec(rho, grho, 0, sc, v1c, v2c);
                s = sx + sc;
                v1 = v1x + v1c;
                v2 = v2x + v2c;
                break;
            default: //SCAN_X,SCAN_C,HSE, and so on
                throw std::domain_error("functional unfinished in "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
        }
        sxc += s;
        v1xc += v1;
        v2xc += v2;
    }

    return;
}

//-----------------------------------------------------------------------
void XC_Functional::gcx_spin(double rhoup, double rhodw, double grhoup2, double grhodw2,
              double &sx, double &v1xup, double &v1xdw, double &v2xup, double &v2xdw)
{
    //--------------------------------------------------------------------
    //     gradient corrections for exchange - Hartree a.u.
    //     Implemented:  Becke88, GGA (PW91), PBE, revPBE
    //
    // use funct
    // USE kinds
    // implicit none

    //     dummy arguments

    // real rhoup, rhodw, grhoup2, grhodw2, sx, v1xup, v1xdw,
    //	v2xup, v2xdw;
    // up and down charge
    // up and down gradient of the charge
    // exchange and correlation energies
    // derivatives of exchange wr. rho
    // derivatives of exchange wr. grho

    // parameter :
    double small = 1.e-10;
    double sxup, sxdw;
    int iflag;

    // exchange
    double rho = rhoup + rhodw;

    sx = 0.00;
    v1xup = 0.00; v2xup = 0.00;
    v1xdw = 0.00; v2xdw = 0.00;
    sxup  = 0.00; sxdw  = 0.00;

    if (rho <= small)
    {
        return;
    }

    // not the correct way to do things, will change later
    // should put exchange and correlation together
    // like the others
    //for(int id : func_id)
    //{
        int id = func_id[0];
        switch( id )
        {
            case XC_GGA_X_B88: //B88
                if (rhoup > small && sqrt(fabs(grhoup2)) > small)
                {
                    XC_Functional::becke88_spin(rhoup, grhoup2, sxup, v1xup, v2xup);
                }
                if (rhodw > small && sqrt(fabs(grhodw2)) > small)
                {
                    XC_Functional::becke88_spin(rhodw, grhodw2, sxdw, v1xdw, v2xdw);
                }
                break;
            case XC_GGA_X_PBE: //PBX
                if (rhoup > small && sqrt(fabs(grhoup2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhoup, 4.0 * grhoup2, 0, sxup, v1xup, v2xup);
                }
                if (rhodw > small && sqrt(fabs(grhodw2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhodw, 4.0 * grhodw2, 0, sxdw, v1xdw, v2xdw);
                }
                break;
            case XC_GGA_X_PBE_R: //revised PBX
                if (rhoup > small && sqrt(fabs(grhoup2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhoup, 4.0 * grhoup2, 1, sxup, v1xup, v2xup);
                }
                if (rhodw > small && sqrt(fabs(grhodw2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhodw, 4.0 * grhodw2, 1, sxdw, v1xdw, v2xdw);
                }
                break;
            case XC_HYB_GGA_XC_PBEH: //PBE0
                if (rhoup > small && sqrt(fabs(grhoup2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhoup, 4.0 * grhoup2, 0, sxup, v1xup, v2xup);
                    sxup *= (1.0 - XC_Functional::hybrid_alpha); 
                    v1xup *= (1.0 - XC_Functional::hybrid_alpha); 
                    v2xup *= (1.0 - XC_Functional::hybrid_alpha);
                }
                if (rhodw > small && sqrt(fabs(grhodw2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhodw, 4.0 * grhodw2, 0, sxdw, v1xdw, v2xdw);
        	    	sxdw *= (1.0 - XC_Functional::hybrid_alpha); 
                    v1xdw *= (1.0 - XC_Functional::hybrid_alpha); 
                    v2xdw *= (1.0 - XC_Functional::hybrid_alpha);            
                }
                break;
            case XC_GGA_X_PBE_SOL: //PBXsol
                if (rhoup > small && sqrt(fabs(grhoup2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhoup, 4.0 * grhoup2, 2, sxup, v1xup, v2xup);
                }
                if (rhodw > small && sqrt(fabs(grhodw2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhodw, 4.0 * grhodw2, 2, sxdw, v1xdw, v2xdw);
                }
                break;
            default:
                sxup = 0.0; sxdw = 0.0;
                v1xup = 0.0; v2xup = 0.0;
                v1xdw = 0.0; v2xdw = 0.0;
        }
        sx = 0.50 * (sxup + sxdw);
        v2xup = 2.0 * v2xup;
        v2xdw = 2.0 * v2xdw;       
    //}

    return;
} //end subroutine gcx_spin

//
//-----------------------------------------------------------------------
void XC_Functional::gcc_spin(double rho, double &zeta, double grho, double &sc,
              double &v1cup, double &v1cdw, double &v2c)
{
    //-------------------------------------------------------------------
    //     gradient corrections for correlations - Hartree a.u.
    //     Implemented:  Perdew86, GGA (PW91), PBE

    // use funct
    // USE kinds
    // implicit none

    //     dummy arguments

    // real(kind=DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
    // the total charge
    // the magnetization
    // the gradient of the charge squared
    // exchange and correlation energies
    // derivatives of correlation wr. rho
    // derivatives of correlation wr. grho

    // parameter :
    double small = 1.0e-10;
	double epsr = 1.0e-6;

    double x;

    sc = 0.00;
    v1cup = 0.00; v1cdw = 0.00;
    v2c = 0.00;
    if (std::abs(zeta) - 1.0 > small || rho <= small || sqrt(std::abs(grho)) <= small)
    {
        return;
    }
    else
    {
        // ... ( - 1.0 + epsr )  <  zeta  <  ( 1.0 - epsr )
        // zeta = SIGN( MIN( ABS( zeta ), ( 1.D0 - epsr ) ) , zeta )
        x = std::min(std::abs(zeta), (1.0 - epsr));
		if(zeta>0)
		{
			zeta = x;
		}
		else
		{
			zeta = -x;
		}
    } //endif

    if(func_id[0]==XC_HYB_GGA_XC_PBEH)
    {
        XC_Functional::pbec_spin(rho, zeta, grho, 1, sc, v1cup, v1cdw, v2c);
        return;
    }

    //for(int id : func_id)
    //{
        int id = func_id[1];
        switch( id )
        {          
            case XC_GGA_C_P86: //P86
                XC_Functional::perdew86_spin(rho, zeta, grho, sc, v1cup, v1cdw, v2c);break;
            case XC_GGA_C_PW91: //PW91_C
                ModuleBase::WARNING_QUIT("xc_wrapper_gcxc","there seems to be something wrong with ggac_spin, better use libxc version instead");break;
                //XC_Functional::ggac_spin(rho, zeta, grho, sc, v1cup, v1cdw, v2c);
            case XC_GGA_C_PBE: //PBC
                XC_Functional::pbec_spin(rho, zeta, grho, 1, sc, v1cup, v1cdw, v2c);break;
            case XC_GGA_C_PBE_SOL: //PBCsol
                XC_Functional::pbec_spin(rho, zeta, grho, 2, sc, v1cup, v1cdw, v2c);break;
        }

    //}
    return;
} //end subroutine gcc_spin

void XC_Functional::gcxc_libxc(const double &rho, const double &grho, double &sxc,
          double &v1xc, double &v2xc)
{
#ifdef USE_LIBXC

    const double small = 1.e-6;
    const double smallg = 1.e-10;
    double s,v1,v2;
    sxc = v1xc = v2xc = 0.0;

    if (rho <= small || grho < smallg)
    {
        return;
    }

    std::vector<xc_func_type> funcs = init_func(XC_UNPOLARIZED);
    

    for(xc_func_type &func : funcs)
    {
        xc_gga_exc_vxc(&func, 1, &rho, &grho, &s, &v1, &v2);
        
        sxc += s * rho;
        v2xc += v2 * 2.0;
        v1xc += v1;
    }
    finish_func(funcs);

    return;
#else
    ModuleBase::WARNING_QUIT("gcxc_libxc","compile with LIBXC to use this subroutine");
#endif
} // end subroutine gcxc_libxc

void XC_Functional::gcxc_spin_libxc(double rhoup, double rhodw, 
        ModuleBase::Vector3<double> gdr1, ModuleBase::Vector3<double> gdr2,
        double &sxc, double &v1xcup, double &v1xcdw, double &v2xcup, double &v2xcdw, double &v2xcud)
{
#ifdef USE_LIBXC
	std::vector<xc_func_type> funcs = init_func(XC_POLARIZED);
    double *rho, *grho, *v1xc, *v2xc, *sgn, s;
    sxc = v1xcup = v1xcdw = 0.0;
    v2xcup = v2xcdw = v2xcud = 0.0;
    rho = new double[2];
    grho= new double[3];
    v1xc= new double[2];
    v2xc= new double[3];
    sgn = new double[2];
    
    rho[0] = rhoup;
    rho[1] = rhodw;
    grho[0] = gdr1.norm2();
    grho[1] = gdr1 * gdr2;
    grho[2] = gdr2.norm2();

    const double rho_threshold = 1E-6;
    const double grho_threshold = 1E-10;

    for(xc_func_type &func : funcs)
    {
        if( func.info->family == XC_FAMILY_GGA || func.info->family == XC_FAMILY_HYB_GGA)
        {
            sgn[0] = sgn[1] = 1.0;
            // call Libxc function: xc_gga_exc_vxc
            xc_gga_exc_vxc( &func, 1, rho, grho, &s, v1xc, v2xc);
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
        }
    }
    delete[] grho;
    delete[] rho;
    delete[] v1xc;
    delete[] v2xc;
    delete[] sgn;
    finish_func(funcs);

    return;
#else
    ModuleBase::WARNING_QUIT("xc_spin_libxc","compile with LIBXC to use this subroutine");
#endif
} 
