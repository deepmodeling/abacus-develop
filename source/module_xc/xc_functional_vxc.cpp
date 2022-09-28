// This file contains interface to xc_functional class:
// 1. v_xc : which takes rho as input, and v_xc as output
// 2. v_xc_libxc : which does the same thing as v_xc, but calling libxc
// NOTE : it is only used for nspin = 1 and 2, the nspin = 4 case is treated in v_xc
// 3. v_xc_meta : which takes rho and tau as input, and v_xc as output

#include "xc_functional.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"
#include "./xc_functional.h"
#include "../src_pw/global.h"

// [etxc, vtxc, v] = XC_Functional::v_xc(...)
std::tuple<double,double,ModuleBase::matrix> XC_Functional::v_xc(
	const int &nrxx, // number of real-space grid
	const int &ncxyz, // total number of charge grid
	const double &omega, // volume of cell
    const double*const*const rho_in,
	const double*const rho_core) // core charge density
{
    ModuleBase::TITLE("XC_Functional","v_xc");
    ModuleBase::timer::tick("XC_Functional","v_xc");

    if( (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 2) && use_libxc)
    {
#ifdef USE_LIBXC
        return v_xc_libxc(nrxx, ncxyz, omega, rho_in, rho_core);
#else
        ModuleBase::WARNING_QUIT("v_xc","compile with LIBXC");
#endif
    }

    //Exchange-Correlation potential Vxc(r) from n(r)
    double etxc = 0.0;
    double vtxc = 0.0;
    ModuleBase::matrix v(GlobalV::NSPIN, nrxx);

    // the square of the e charge
    // in Rydeberg unit, so * 2.0.
    double e2 = 2.0;

    double rhox = 0.0;
    double arhox = 0.0;
    double zeta = 0.0;
    double exc = 0.0;
    double vxc[2];

    int ir, is;

    double vanishing_charge = 1.0e-10;

    if (GlobalV::NSPIN == 1 || ( GlobalV::NSPIN ==4 && !GlobalV::DOMAG && !GlobalV::DOMAG_Z))
    {
        // spin-unpolarized case
        for (int ir = 0;ir < nrxx;ir++)
        {
	    // total electron charge density
            rhox = rho_in[0][ir] + rho_core[ir];
            arhox = abs(rhox);
            if (arhox > vanishing_charge)
            {
                XC_Functional::xc(arhox, exc, vxc[0]);
                v(0,ir) = e2 * vxc[0];
				// consider the total charge density
                etxc += e2 * exc * rhox;
				// only consider rho_in
                vtxc += v(0, ir) * rho_in[0][ir];
            } // endif
        } //enddo
    }
    else if(GlobalV::NSPIN ==2)
    {
        // spin-polarized case
        for (ir = 0;ir < nrxx;ir++)
        {
            rhox = rho_in[0][ir] + rho_in[1][ir] + rho_core[ir]; //HLX(05-29-06): bug fixed
            arhox = abs(rhox);

            if (arhox > vanishing_charge)
            {
                zeta = (rho_in[0][ir] - rho_in[1][ir]) / arhox; //HLX(05-29-06): bug fixed

                if (abs(zeta)  > 1.0)
                {
                    zeta = (zeta > 0.0) ? 1.0 : (-1.0);
                }

                double rhoup = arhox * (1.0+zeta) / 2.0;
                double rhodw = arhox * (1.0-zeta) / 2.0;
                XC_Functional::xc_spin(arhox, zeta, exc, vxc[0], vxc[1]);

                for (is = 0;is < GlobalV::NSPIN;is++)
                {
                    v(is, ir) = e2 * vxc[is];
                }

                etxc += e2 * exc * rhox;
                vtxc += v(0, ir) * rho_in[0][ir] + v(1, ir) * rho_in[1][ir];
            }
        }
    }
    else if(GlobalV::NSPIN == 4)//noncollinear case added by zhengdy
    {
        for( ir = 0;ir<nrxx; ir++)
        {
            double amag = sqrt( pow(rho_in[1][ir],2) + pow(rho_in[2][ir],2) + pow(rho_in[3][ir],2) );

            rhox = rho_in[0][ir] + rho_core[ir];

            arhox = abs( rhox );

            if ( arhox > vanishing_charge )
            {
                zeta = amag / arhox;

                if ( abs( zeta ) > 1.0 )
                {
                    zeta = (zeta > 0.0) ? 1.0 : (-1.0);
                }//end if

                if(use_libxc)
                {
                    double rhoup = arhox * (1.0+zeta) / 2.0;
                    double rhodw = arhox * (1.0-zeta) / 2.0;
                    XC_Functional::xc_spin_libxc(rhoup, rhodw, exc, vxc[0], vxc[1]);
                }
                else
                {
                    double rhoup = arhox * (1.0+zeta) / 2.0;
                    double rhodw = arhox * (1.0-zeta) / 2.0;
                    XC_Functional::xc_spin(arhox, zeta, exc, vxc[0], vxc[1]);
                }

                etxc += e2 * exc * rhox;

                v(0, ir) = e2*( 0.5 * ( vxc[0] + vxc[1]) );
                vtxc += v(0,ir) * rho_in[0][ir];

                double vs = 0.5 * ( vxc[0] - vxc[1] );
                if ( amag > vanishing_charge )
                {
                    for(int ipol = 1;ipol< 4;ipol++)
                    {
                        v(ipol, ir) = e2 * vs * rho_in[ipol][ir] / amag;
                        vtxc += v(ipol,ir) * rho_in[ipol][ir];
                    }//end do
                }//end if
            }//end if
        }//end do
    }//end if
    // energy terms, local-density contributions

    // add gradient corrections (if any)
    // mohan modify 2009-12-15

    // the dummy variable dum contains gradient correction to stress
    // which is not used here
    std::vector<double> dum;
    gradcorr(etxc, vtxc, v, dum);

    // parallel code : collect vtxc,etxc
    // mohan add 2008-06-01
#ifdef __MPI
    Parallel_Reduce::reduce_double_pool( etxc );
    Parallel_Reduce::reduce_double_pool( vtxc );
#endif
    etxc *= omega / ncxyz;
    vtxc *= omega / ncxyz;

    ModuleBase::timer::tick("XC_Functional","v_xc");
    return std::make_tuple(etxc, vtxc, std::move(v));
}

#ifdef USE_LIBXC

std::tuple<double,double,ModuleBase::matrix> XC_Functional::v_xc_libxc(
        const int &nrxx, // number of real-space grid
        const int &ncxyz, // total number of charge grid
        const double &omega, // volume of cell
        const double * const * const rho_in,
        const double * const rho_core_in)
{
    ModuleBase::TITLE("XC_Functional","v_xc");
    ModuleBase::timer::tick("XC_Functional","v_xc");

    const int nspin = GlobalV::NSPIN;

    double etxc = 0.0;
    double vtxc = 0.0;
    ModuleBase::matrix v(nspin,nrxx);

    //----------------------------------------------------------
    // xc_func_type is defined in Libxc package
    // to understand the usage of xc_func_type,
    // use can check on website, for example:
    // https://www.tddft.org/programs/libxc/manual/libxc-5.1.x/
    //----------------------------------------------------------

    std::vector<xc_func_type> funcs = init_func( ( (1==nspin) ? XC_UNPOLARIZED:XC_POLARIZED ) );

    bool is_gga = false;
    for( xc_func_type &func : funcs )
    {
        switch( func.info->family )
        {
            case XC_FAMILY_GGA:
            case XC_FAMILY_HYB_GGA:
                is_gga = true;
                break;
        }
    }

    // converting rho
    std::vector<double> rho;
    rho.resize(nrxx*nspin);
    for( int is=0; is!=nspin; ++is )
    {
        for( int ir=0; ir!=nrxx; ++ir )
        {
            rho[ir*nspin+is] = rho_in[is][ir] + 1.0/nspin*rho_core_in[ir];
        }
    }

    std::vector<std::vector<ModuleBase::Vector3<double>>> gdr;
    std::vector<double> sigma;
    if(is_gga)
    {
        // calculating grho
        gdr.resize( nspin );
        for( int is=0; is!=nspin; ++is )
        {
            std::vector<double> rhor(nrxx);
            for(int ir=0; ir<nrxx; ++ir)
            {
                rhor[ir] = rho[ir*nspin+is];
            }
            //------------------------------------------
            // initialize the charge density array in reciprocal space
            // bring electron charge density from real space to reciprocal space
            //------------------------------------------
            std::vector<std::complex<double>> rhog(GlobalC::rhopw->npw);
            GlobalC::CHR.set_rhog(rhor.data(), rhog.data());

            //-------------------------------------------
            // compute the gradient of charge density and
            // store the gradient in gdr[is]
            //-------------------------------------------
            gdr[is].resize(nrxx);
            XC_Functional::grad_rho(rhog.data(), gdr[is].data(), GlobalC::rhopw);
        }

        // converting grho
        sigma.resize( nrxx * ((1==nspin)?1:3) );

        if( 1==nspin )
        {
            for( int ir=0; ir!=nrxx; ++ir )
            {
                sigma[ir] = gdr[0][ir]*gdr[0][ir];
            } 
        }
        else
        {
            for( int ir=0; ir!=nrxx; ++ir )
            {
                sigma[ir*3]   = gdr[0][ir]*gdr[0][ir];
                sigma[ir*3+1] = gdr[0][ir]*gdr[1][ir];
                sigma[ir*3+2] = gdr[1][ir]*gdr[1][ir];
            }
        }
    }

    std::vector<double> exc   ( nrxx                    );
    std::vector<double> vrho  ( nrxx * nspin            );
    std::vector<double> vsigma( nrxx * ((1==nspin)?1:3) );

    for( xc_func_type &func : funcs )
    {
        // jiyy add for threshold
        const double rho_threshold = 1E-6;
        const double grho_threshold = 1E-10;
        // sgn for threshold mask
        std::vector<double> sgn( nrxx * nspin, 1.0);

        // in the case of GGA correlation for polarized case,
        // a cutoff for grho is required to ensure that libxc gives reasonable results
        if(nspin==2 && func.info->family != XC_FAMILY_LDA && func.info->kind==XC_CORRELATION)
        {
            for( size_t ir=0; ir!=nrxx; ++ir )
            {
                if ( rho[ir*2]<rho_threshold || sqrt(abs(sigma[ir*3]))<grho_threshold )
                    sgn[ir*2] = 0.0;
                if ( rho[ir*2+1]<rho_threshold || sqrt(abs(sigma[ir*3+2]))<grho_threshold )
                    sgn[ir*2+1] = 0.0;
            }
        }

        xc_func_set_dens_threshold(&func, rho_threshold);
        switch( func.info->family )
        {
            case XC_FAMILY_LDA:
                // call Libxc function: xc_lda_exc_vxc
                xc_lda_exc_vxc( &func, nrxx, rho.data(), 
                    exc.data(), vrho.data() );
                break;
            case XC_FAMILY_GGA:
            case XC_FAMILY_HYB_GGA:
                // call Libxc function: xc_gga_exc_vxc
                xc_gga_exc_vxc( &func, nrxx, rho.data(), sigma.data(),
                    exc.data(), vrho.data(), vsigma.data() );
                break;
            default:
                throw std::domain_error("func.info->family ="+ModuleBase::GlobalFunc::TO_STRING(func.info->family)
                    +" unfinished in "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
                break;
        }

        for( int is=0; is!=nspin; ++is )
        {
            for( int ir=0; ir!= nrxx; ++ir )
            {
                etxc += ModuleBase::e2 * exc[ir] * rho[ir*nspin+is] * sgn[ir*nspin+is];
            }   
        }

        for( int is=0; is!=nspin; ++is )
        {
            for( int ir=0; ir!= nrxx; ++ir )
            {
                const double v_tmp = ModuleBase::e2 * vrho[ir*nspin+is] * sgn[ir*nspin+is];
                v(is,ir) += v_tmp;
                vtxc += v_tmp * rho_in[is][ir];
            }
        }

        if(func.info->family == XC_FAMILY_GGA || func.info->family == XC_FAMILY_HYB_GGA)
        {
            std::vector<std::vector<ModuleBase::Vector3<double>>> h( nspin, std::vector<ModuleBase::Vector3<double>>(nrxx) );
            if( 1==nspin )
            {
                for( int ir=0; ir!= nrxx; ++ir )
                {
                    h[0][ir] = 2.0 * gdr[0][ir] * vsigma[ir] * 2.0 * sgn[ir];
                }
            }
            else
            {
                for( int ir=0; ir!= nrxx; ++ir )
                {
                    h[0][ir] = 2.0 * (gdr[0][ir] * vsigma[ir*3  ] * sgn[ir*2  ] * 2.0 
                                    + gdr[1][ir] * vsigma[ir*3+1] * sgn[ir*2]   * sgn[ir*2+1]);
                    h[1][ir] = 2.0 * (gdr[1][ir] * vsigma[ir*3+2] * sgn[ir*2+1] * 2.0 
                                    + gdr[0][ir] * vsigma[ir*3+1] * sgn[ir*2]   * sgn[ir*2+1]);
                }
            }

            // define two dimensional array dh [ nspin, nrxx ]
            std::vector<std::vector<double>> dh(nspin, std::vector<double>( nrxx));
            for( int is=0; is!=nspin; ++is )
            {
                XC_Functional::grad_dot( ModuleBase::GlobalFunc::VECTOR_TO_PTR(h[is]), ModuleBase::GlobalFunc::VECTOR_TO_PTR(dh[is]), GlobalC::rhopw );
            }

            for( int is=0; is!=nspin; ++is )
            {
                for( int ir=0; ir!= nrxx; ++ir )
                {
                    vtxc -= dh[is][ir] * rho[ir*nspin+is];
                    v(is,ir) -= dh[is][ir];
                }
            }
        }
    }

    //-------------------------------------------------
    // for MPI, reduce the exchange-correlation energy
    //-------------------------------------------------
#ifdef __MPI
    Parallel_Reduce::reduce_double_pool( etxc );
    Parallel_Reduce::reduce_double_pool( vtxc );
#endif

    etxc *= omega / ncxyz;
    vtxc *= omega / ncxyz;

    ModuleBase::timer::tick("XC_Functional","v_xc");
    return std::make_tuple( etxc, vtxc, std::move(v) );
}

//the interface to libxc xc_mgga_exc_vxc(xc_func,n,rho,grho,laplrho,tau,e,v1,v2,v3,v4)
//xc_func : LIBXC data type, contains information on xc functional
//n: size of array, nspin*nnr
//rho,grho,laplrho: electron density, its gradient and laplacian
//tau(kin_r): kinetic energy density
//e: energy density
//v1-v4: derivative of energy density w.r.t rho, gradient, laplacian and tau
//v1 and v2 are combined to give v; v4 goes into vofk

//XC_POLARIZED, XC_UNPOLARIZED: internal flags used in LIBXC, denote the polarized(nspin=1) or unpolarized(nspin=2) calculations, definition can be found in xc.h from LIBXC

// [etxc, vtxc, v, vofk] = XC_Functional::v_xc(...)
tuple<double,double,ModuleBase::matrix,ModuleBase::matrix> XC_Functional::v_xc_meta(
	const int &nrxx, // number of real-space grid
	const int &ncxyz, // total number of charge grid
	const double &omega, // volume of cell
	const double * const * const rho_in,
	const double * const rho_core_in,
	const double * const * const kin_r_in)
{
    ModuleBase::TITLE("XC_Functional","v_xc");
    ModuleBase::timer::tick("XC_Functional","v_xc");

    if(GlobalV::NSPIN==4)
    {
        ModuleBase::WARNING_QUIT("v_xc_meta","meta-GGA has not been implemented for nspin = 4 yet");
    }

    double e2 = 2.0;

	//output of the subroutine
    double etxc = 0.0;
    double vtxc = 0.0;
	ModuleBase::matrix v(GlobalV::NSPIN,nrxx);
	ModuleBase::matrix vofk(GlobalV::NSPIN,nrxx);

	//----------------------------------------------------------
	// xc_func_type is defined in Libxc package
	// to understand the usage of xc_func_type,
	// use can check on website, for example:
	// https://www.tddft.org/programs/libxc/manual/libxc-5.1.x/
	//----------------------------------------------------------
	
    const int nspin = GlobalV::NSPIN;
    std::vector<xc_func_type> funcs = init_func( ( (1==nspin) ? XC_UNPOLARIZED:XC_POLARIZED ) );

    // converting rho
    std::vector<double> rho;
    rho.resize(nrxx*nspin);
    for( int is=0; is!=nspin; ++is )
    {
        for( int ir=0; ir!=nrxx; ++ir )
        {
            rho[ir*nspin+is] = rho_in[is][ir] + 1.0/nspin*rho_core_in[ir];
        }
    }

    std::vector<std::vector<ModuleBase::Vector3<double>>> gdr;
    std::vector<double> sigma;

    // calculating grho
    gdr.resize( nspin );
    for( int is=0; is!=nspin; ++is )
    {
        std::vector<double> rhor(nrxx);
        for(int ir=0; ir<nrxx; ++ir)
        {
            rhor[ir] = rho[ir*nspin+is];
        }
        //------------------------------------------
        // initialize the charge density array in reciprocal space
        // bring electron charge density from real space to reciprocal space
        //------------------------------------------
        std::vector<std::complex<double>> rhog(GlobalC::rhopw->npw);
        GlobalC::CHR.set_rhog(rhor.data(), rhog.data());

        //-------------------------------------------
        // compute the gradient of charge density and
        // store the gradient in gdr[is]
        //-------------------------------------------
        gdr[is].resize(nrxx);
        XC_Functional::grad_rho(rhog.data(), gdr[is].data(), GlobalC::rhopw);
    }

    // converting grho
    sigma.resize( nrxx * ((1==nspin)?1:3) );

    if( 1==nspin )
    {
        for( int ir=0; ir!=nrxx; ++ir )
        {
            sigma[ir] = gdr[0][ir]*gdr[0][ir];
        } 
    }
    else
    {
        for( int ir=0; ir!=nrxx; ++ir )
        {
            sigma[ir*3]   = gdr[0][ir]*gdr[0][ir];
            sigma[ir*3+1] = gdr[0][ir]*gdr[1][ir];
            sigma[ir*3+2] = gdr[1][ir]*gdr[1][ir];
        }
    }

    //converting kin_r
    std::vector<double> kin_r;
    kin_r.resize(nrxx*nspin);
    for( int is=0; is!=nspin; ++is )
    {
        for( int ir=0; ir!=nrxx; ++ir )
        {
            kin_r[ir*nspin+is] = kin_r_in[is][ir] / 2.0;
        }
    }

    std::vector<double> exc    ( nrxx                    );
    std::vector<double> vrho   ( nrxx * nspin            );
    std::vector<double> vsigma ( nrxx * ((1==nspin)?1:3) );
	std::vector<double> vtau   ( nrxx * nspin            );
    std::vector<double> vlapl  ( nrxx * nspin            );

    const double rho_th  = 1e-8;
    const double grho_th = 1e-12;
    const double tau_th  = 1e-8;
    // sgn for threshold mask
    std::vector<double> sgn( nrxx * nspin, 1.0);

    if(nspin == 1)
    {
        for( size_t ir=0; ir!=nrxx; ++ir )
        {
            if ( rho[ir]<rho_th || sqrt(abs(sigma[ir]))<grho_th || abs(kin_r[ir]<tau_th))
            {
                sgn[ir] = 0.0;
            }
        }
    }
    else
    {
        for( size_t ir=0; ir!=nrxx; ++ir )
        {
            if ( rho[ir*2]<rho_th || sqrt(abs(sigma[ir*3]))<grho_th || abs(kin_r[ir*2]<tau_th))
                sgn[ir*2] = 0.0;
            if ( rho[ir*2+1]<rho_th || sqrt(abs(sigma[ir*3+2]))<grho_th || abs(kin_r[ir*2+1]<tau_th))
                sgn[ir*2+1] = 0.0;
        }
    }

    for ( xc_func_type &func : funcs )
    {
        assert(func.info->family == XC_FAMILY_MGGA);
        xc_mgga_exc_vxc(&func, nrxx, rho.data(), sigma.data(), sigma.data(),
            kin_r.data(), exc.data(), vrho.data(), vsigma.data(), vlapl.data(), vtau.data());

        //process etxc
        for( int is=0; is!=nspin; ++is )
        {
            for( int ir=0; ir!= nrxx; ++ir )
            {
#ifdef __EXX
                if (func.info->number == XC_MGGA_X_SCAN && get_func_type() == 5)
                {
                    exc[ir] *= (1.0 - GlobalC::exx_global.info.hybrid_alpha);
                }
#endif
                etxc += ModuleBase::e2 * exc[ir] * rho[ir*nspin+is]  * sgn[ir*nspin+is];
            }   
        }

        //process vtxc
        for( int is=0; is!=nspin; ++is )
        {
            for( int ir=0; ir!= nrxx; ++ir )
            {
#ifdef __EXX
                if (func.info->number == XC_MGGA_X_SCAN && get_func_type() == 5)
                {
                    vrho[ir*nspin+is] *= (1.0 - GlobalC::exx_global.info.hybrid_alpha);
                }
#endif
                const double v_tmp = ModuleBase::e2 * vrho[ir*nspin+is]  * sgn[ir*nspin+is];
                v(is,ir) += v_tmp;
                vtxc += v_tmp * rho_in[is][ir];
            }
        }

        //process vsigma
        std::vector<std::vector<ModuleBase::Vector3<double>>> h( nspin, std::vector<ModuleBase::Vector3<double>>(nrxx) );
        if( 1==nspin )
        {
            for( int ir=0; ir!= nrxx; ++ir )
            {
#ifdef __EXX
                if (func.info->number == XC_MGGA_X_SCAN && get_func_type() == 5)
                {
                    vsigma[ir] *= (1.0 - GlobalC::exx_global.info.hybrid_alpha);
                }
#endif
                h[0][ir] = 2.0 * gdr[0][ir] * vsigma[ir] * 2.0 * sgn[ir];
            }
        }
        else
        {
            for( int ir=0; ir!= nrxx; ++ir )
            {
#ifdef __EXX
                if (func.info->number == XC_MGGA_X_SCAN && get_func_type() == 5)
                {
                    vsigma[ir*3] *= (1.0 - GlobalC::exx_global.info.hybrid_alpha);
                    vsigma[ir*3+1] *= (1.0 - GlobalC::exx_global.info.hybrid_alpha);
                    vsigma[ir*3+2] *= (1.0 - GlobalC::exx_global.info.hybrid_alpha);
                }
#endif
                h[0][ir] = 2.0 * (gdr[0][ir] * vsigma[ir*3  ] * sgn[ir*2  ] * 2.0 
                                + gdr[1][ir] * vsigma[ir*3+1] * sgn[ir*2]   * sgn[ir*2+1]);
                h[1][ir] = 2.0 * (gdr[1][ir] * vsigma[ir*3+2] * sgn[ir*2+1] * 2.0 
                                + gdr[0][ir] * vsigma[ir*3+1] * sgn[ir*2]   * sgn[ir*2+1]);
            }
        }

        // define two dimensional array dh [ nspin, nrxx ]
        std::vector<std::vector<double>> dh(nspin, std::vector<double>( nrxx));
        for( int is=0; is!=nspin; ++is )
        {
            XC_Functional::grad_dot( ModuleBase::GlobalFunc::VECTOR_TO_PTR(h[is]), ModuleBase::GlobalFunc::VECTOR_TO_PTR(dh[is]), GlobalC::rhopw );
        }

        for( int is=0; is!=nspin; ++is )
        {
            for( int ir=0; ir!= nrxx; ++ir )
            {
                vtxc -= dh[is][ir] * rho[ir*nspin+is];
                v(is,ir) -= dh[is][ir];
            }
        }

        //process vtau
        for( int is=0; is!=nspin; ++is )
        {
            for( int ir=0; ir!= nrxx; ++ir )
            {
#ifdef __EXX
                if (func.info->number == XC_MGGA_X_SCAN && get_func_type() == 5)
                {
                    vtau[ir*nspin+is] *= (1.0 - GlobalC::exx_global.info.hybrid_alpha);
                }
#endif
                vofk(is,ir) += vtau[ir*nspin+is]  * sgn[ir*nspin+is];
            }
        }
	}
		
	//-------------------------------------------------
	// for MPI, reduce the exchange-correlation energy
	//-------------------------------------------------
#ifdef __MPI
	Parallel_Reduce::reduce_double_pool( etxc );
	Parallel_Reduce::reduce_double_pool( vtxc );
#endif

    etxc *= omega / ncxyz;
    vtxc *= omega / ncxyz;

    ModuleBase::timer::tick("XC_Functional","v_xc_meta");
	return std::make_tuple( etxc, vtxc, move(v), move(vofk) );

}

#endif	//ifdef USE_LIBXC