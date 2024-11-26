#ifdef USE_LIBXC

#include "xc_functional.h"
#include "xc_functional_libxc.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_base/global_variable.h"
#include "module_parameter/parameter.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"

#include <xc.h>

#include <vector>

std::tuple<double,double,ModuleBase::matrix> XC_Functional_Libxc::v_xc_libxc(		// Peize Lin update for nspin==4 at 2023.01.14
        const std::vector<int> &func_id,
        const int &nrxx, // number of real-space grid
        const double &omega, // volume of cell
        const double tpiba,
        const Charge* const chr)
{
    ModuleBase::TITLE("XC_Functional_Libxc","v_xc_libxc");
    ModuleBase::timer::tick("XC_Functional_Libxc","v_xc_libxc");

    const int nspin =
        (PARAM.inp.nspin == 1 || ( PARAM.inp.nspin ==4 && !PARAM.globalv.domag && !PARAM.globalv.domag_z))
        ? 1 : 2;

    //----------------------------------------------------------
    // xc_func_type is defined in Libxc package
    // to understand the usage of xc_func_type,
    // use can check on website, for example:
    // https://www.tddft.org/programs/libxc/manual/libxc-5.1.x/
    //----------------------------------------------------------

    std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func( func_id, (1==nspin) ? XC_UNPOLARIZED:XC_POLARIZED );

    const bool is_gga = [&funcs]()
    {
        for( xc_func_type &func : funcs )
        {
            switch( func.info->family )
            {
                case XC_FAMILY_GGA:
                case XC_FAMILY_HYB_GGA:
                    return true;
            }
        }
        return false;
    }();

    // converting rho
    std::vector<double> rho;
    std::vector<double> amag;
    if(1==nspin || 2==PARAM.inp.nspin)
    {
        rho = XC_Functional_Libxc::convert_rho(nspin, nrxx, chr);
    }
    else
    {
        std::tuple<std::vector<double>,std::vector<double>> rho_amag = XC_Functional_Libxc::convert_rho_amag_nspin4(nspin, nrxx, chr);
        rho = std::get<0>(std::move(rho_amag));
        amag = std::get<1>(std::move(rho_amag));
    }

    std::vector<std::vector<ModuleBase::Vector3<double>>> gdr;
    std::vector<double> sigma;
    if(is_gga)
    {
        gdr = XC_Functional_Libxc::cal_gdr(nspin, nrxx, rho, tpiba, chr);
        sigma = XC_Functional_Libxc::convert_sigma(gdr);
    }

    double etxc = 0.0;
    double vtxc = 0.0;
    ModuleBase::matrix v(nspin,nrxx);

    for( xc_func_type &func : funcs )
    {
        // jiyy add for threshold
        constexpr double rho_threshold = 1E-6;
        constexpr double grho_threshold = 1E-10;

        xc_func_set_dens_threshold(&func, rho_threshold);

        // sgn for threshold mask
        const std::vector<double> sgn = XC_Functional_Libxc::cal_sgn(rho_threshold, grho_threshold, func, nspin, nrxx, rho, sigma);

        std::vector<double> exc   ( nrxx                    );
        std::vector<double> vrho  ( nrxx * nspin            );
        std::vector<double> vsigma( nrxx * ((1==nspin)?1:3) );
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
                throw std::domain_error("func.info->family ="+std::to_string(func.info->family)
                    +" unfinished in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
                break;
        }

        etxc += XC_Functional_Libxc::convert_etxc(nspin, nrxx, sgn, rho, exc);
        const std::pair<double,ModuleBase::matrix> vtxc_v = XC_Functional_Libxc::convert_vtxc_v(
            func, nspin, nrxx,
            sgn, rho, gdr,
            vrho, vsigma,
            tpiba, chr);
        vtxc += std::get<0>(vtxc_v);
        v += std::get<1>(vtxc_v);
    } // end for( xc_func_type &func : funcs )

    if(4==PARAM.inp.nspin)
    {
        v = XC_Functional_Libxc::convert_v_nspin4(nrxx, chr, amag, v);
    }

    //-------------------------------------------------
    // for MPI, reduce the exchange-correlation energy
    //-------------------------------------------------
    #ifdef __MPI
    Parallel_Reduce::reduce_pool(etxc);
    Parallel_Reduce::reduce_pool(vtxc);
    #endif

    etxc *= omega / chr->rhopw->nxyz;
    vtxc *= omega / chr->rhopw->nxyz;

    XC_Functional_Libxc::finish_func(funcs);

    ModuleBase::timer::tick("XC_Functional_Libxc","v_xc_libxc");
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
std::tuple<double,double,ModuleBase::matrix,ModuleBase::matrix> XC_Functional_Libxc::v_xc_meta(
    const std::vector<int> &func_id,
    const int &nrxx, // number of real-space grid
    const double &omega, // volume of cell
    const double tpiba,
    const Charge* const chr)
{
    ModuleBase::TITLE("XC_Functional_Libxc","v_xc_meta");
    ModuleBase::timer::tick("XC_Functional_Libxc","v_xc_meta");

    double e2 = 2.0;

    //output of the subroutine
    double etxc = 0.0;
    double vtxc = 0.0;
    ModuleBase::matrix v(PARAM.inp.nspin,nrxx);
    ModuleBase::matrix vofk(PARAM.inp.nspin,nrxx);

    //----------------------------------------------------------
    // xc_func_type is defined in Libxc package
    // to understand the usage of xc_func_type,
    // use can check on website, for example:
    // https://www.tddft.org/programs/libxc/manual/libxc-5.1.x/
    //----------------------------------------------------------

    const int nspin = PARAM.inp.nspin;
    std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func( func_id, ( (1==nspin) ? XC_UNPOLARIZED:XC_POLARIZED ) );

    const std::vector<double> rho = XC_Functional_Libxc::convert_rho(nspin, nrxx, chr);
    const std::vector<std::vector<ModuleBase::Vector3<double>>> gdr
        = XC_Functional_Libxc::cal_gdr(nspin, nrxx, rho, tpiba, chr);
    const std::vector<double> sigma = XC_Functional_Libxc::convert_sigma(gdr);

    //converting kin_r
    std::vector<double> kin_r;
    kin_r.resize(nrxx*nspin);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
    for( int is=0; is<nspin; ++is )
    {
        for( int ir=0; ir<nrxx; ++ir )
        {
            kin_r[ir*nspin+is] = chr->kin_r[is][ir] / 2.0;
        }
    }

    std::vector<double> exc    ( nrxx                    );
    std::vector<double> vrho   ( nrxx * nspin            );
    std::vector<double> vsigma ( nrxx * ((1==nspin)?1:3) );
    std::vector<double> vtau   ( nrxx * nspin            );
    std::vector<double> vlapl  ( nrxx * nspin            );

    constexpr double rho_th  = 1e-8;
    constexpr double grho_th = 1e-12;
    constexpr double tau_th  = 1e-8;
    // sgn for threshold mask
    std::vector<double> sgn( nrxx * nspin);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int i = 0; i < nrxx * nspin; ++i)
    {
        sgn[i] = 1.0;
    }

    if(nspin == 1)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for( int ir=0; ir<nrxx; ++ir )
        {
            if ( rho[ir]<rho_th || sqrt(std::abs(sigma[ir]))<grho_th || std::abs(kin_r[ir])<tau_th)
            {
                sgn[ir] = 0.0;
            }
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for( int ir=0; ir<nrxx; ++ir )
        {
            if ( rho[ir*2]<rho_th || sqrt(std::abs(sigma[ir*3]))<grho_th || std::abs(kin_r[ir*2])<tau_th)
                sgn[ir*2] = 0.0;
            if ( rho[ir*2+1]<rho_th || sqrt(std::abs(sigma[ir*3+2]))<grho_th || std::abs(kin_r[ir*2+1])<tau_th)
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
#ifdef _OPENMP
#pragma omp parallel for reduction(+:etxc) schedule(static, 256)
#endif
            for( int ir=0; ir< nrxx; ++ir )
            {
#ifdef __EXX
                if (func.info->number == XC_MGGA_X_SCAN && XC_Functional::get_func_type() == 5)
                {
                    exc[ir] *= (1.0 - XC_Functional::get_hybrid_alpha());
                }
#endif
                etxc += ModuleBase::e2 * exc[ir] * rho[ir*nspin+is]  * sgn[ir*nspin+is];
            }
        }

        //process vtxc
#ifdef _OPENMP
#pragma omp parallel for collapse(2) reduction(+:vtxc) schedule(static, 256)
#endif
        for( int is=0; is<nspin; ++is )
        {
            for( int ir=0; ir< nrxx; ++ir )
            {
#ifdef __EXX
                if (func.info->number == XC_MGGA_X_SCAN && XC_Functional::get_func_type() == 5)
                {
                    vrho[ir*nspin+is] *= (1.0 - XC_Functional::get_hybrid_alpha());
                }
#endif
                const double v_tmp = ModuleBase::e2 * vrho[ir*nspin+is]  * sgn[ir*nspin+is];
                v(is,ir) += v_tmp;
                vtxc += v_tmp * chr->rho[is][ir];
            }
        }

        //process vsigma
        std::vector<std::vector<ModuleBase::Vector3<double>>> h(
            nspin,
            std::vector<ModuleBase::Vector3<double>>(nrxx) );
        if( 1==nspin )
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for( int ir=0; ir< nrxx; ++ir )
            {
#ifdef __EXX
                if (func.info->number == XC_MGGA_X_SCAN && XC_Functional::get_func_type() == 5)
                {
                    vsigma[ir] *= (1.0 - XC_Functional::get_hybrid_alpha());
                }
#endif
                h[0][ir] = 2.0 * gdr[0][ir] * vsigma[ir] * 2.0 * sgn[ir];
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 64)
#endif
            for( int ir=0; ir< nrxx; ++ir )
            {
#ifdef __EXX
                if (func.info->number == XC_MGGA_X_SCAN && XC_Functional::get_func_type() == 5)
                {
                    vsigma[ir*3]   *= (1.0 - XC_Functional::get_hybrid_alpha());
                    vsigma[ir*3+1] *= (1.0 - XC_Functional::get_hybrid_alpha());
                    vsigma[ir*3+2] *= (1.0 - XC_Functional::get_hybrid_alpha());
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
            XC_Functional::grad_dot( h[is].data(),
                dh[is].data(), chr->rhopw,
                tpiba);
        }

        double rvtxc = 0.0;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) reduction(+:rvtxc) schedule(static, 256)
#endif
        for( int is=0; is<nspin; ++is )
        {
            for( int ir=0; ir< nrxx; ++ir )
            {
                rvtxc += dh[is][ir] * rho[ir*nspin+is];
                v(is,ir) -= dh[is][ir];
            }
        }
        vtxc -= rvtxc;

        //process vtau
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
        for( int is=0; is<nspin; ++is )
        {
            for( int ir=0; ir< nrxx; ++ir )
            {
#ifdef __EXX
                if (func.info->number == XC_MGGA_X_SCAN && XC_Functional::get_func_type() == 5)
                {
                    vtau[ir*nspin+is] *= (1.0 - XC_Functional::get_hybrid_alpha());
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
    Parallel_Reduce::reduce_pool(etxc);
    Parallel_Reduce::reduce_pool(vtxc);
#endif

    etxc *= omega / chr->rhopw->nxyz;
    vtxc *= omega / chr->rhopw->nxyz;

    XC_Functional_Libxc::finish_func(funcs);

    ModuleBase::timer::tick("XC_Functional_Libxc","v_xc_meta");
    return std::make_tuple( etxc, vtxc, std::move(v), std::move(vofk) );
}

#endif