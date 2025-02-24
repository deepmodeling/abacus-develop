#ifdef USE_LIBXC

#include "xc_functional.h"
#include "xc_functional_libxc.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_base/global_variable.h"
#include "module_parameter/parameter.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"

#include "NCLibxc/NCLibxc.h"

#include <xc.h>

#include <vector>
// ...existing code...
#include <iomanip>

using namespace NCXC;


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
    ModuleBase::matrix v_nspin4(4,nrxx);

    std::vector<double> torque(3 * nrxx, 0.0);

    
    bool use_v_nspin4 = false;
    //double etxc_col = 0.0; //for collinear tests
    //double vtxc_col = 0.0; //for collinear tests
   // double rho0_0 = 0.0; //for collinear tests
   // double rho1_0 = 0.0; //for collinear tests
   // double v1_col0 = 0.0; //for collinear tests
    //double v2_col0 = 0.0; //for collinear tests


    for( xc_func_type &func : funcs )
    {
        // jiyy add for threshold
        constexpr double rho_threshold = 0.0;
        constexpr double grho_threshold = 0.0;

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
        if(PARAM.inp.multicolin)//  added by Xiaoyu Zhang, Peking University, 2024.10.09.  multicollinear method 
        {
            use_v_nspin4 = true;
            etxc=0;
            vtxc=0;
            
            std::complex<double> twoi(0.0, 2.0);
            std::complex<double> two(2.0, 0.0);
            double e2=2.0;


            NCLibxc::print_NCLibxc();
            if(!is_gga)
            { //LDA 
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024) reduction(+:etxc) reduction(+:vtxc)
#endif
                for(int ir = 0;ir<nrxx; ++ir){
                    double exc = 0.0;
                    for(int ipol=0;ipol<4;++ipol){
                        v_nspin4(ipol, ir) = 0;
                    }
                    std::vector<double> n = {chr->rho[0][ir] + chr->rho_core[ir]};
                    std::vector<double> mx = {chr->rho[1][ir]};
                    std::vector<double> my = {chr->rho[2][ir]};
                    std::vector<double> mz = {chr->rho[3][ir]};
                    double amag = sqrt( pow(chr->rho[1][ir],2) + pow(chr->rho[2][ir],2) + pow(chr->rho[3][ir],2) );
                     if (n[0] - amag <= 0.0) { //ensure the rhoup and rhodn to libxc are positive
                        continue;
                    }
                    std::vector<double> E_MC;
                    std::vector<Matrix2x2> V_MC;
                   for(const int &id : func_id){
                        //auto [E_MC, V_MC] = NCLibxc::lda_mc(id, n, mx, my, mz);
                       std::tie(E_MC, V_MC) = NCLibxc::lda_mc(id, n, mx, my, mz);
                        exc = e2*E_MC[0];
                        auto v00 = V_MC[0][0][0];
                        auto v01 = V_MC[0][0][1];
                        auto v10 = V_MC[0][1][0];
                        auto v11 = V_MC[0][1][1];
                        v_nspin4(0, ir) += std::real(e2*(v00+v11)/two);
                        v_nspin4(1, ir) += std::real(e2*(v01+v10)/two);
                        v_nspin4(2, ir) += std::real(e2*(v10-v01)/twoi);
                        v_nspin4(3, ir) += std::real(e2*(v00-v11)/two);
                        etxc += exc * n[0];
                        vtxc += std::real(e2*(v00+v11)/two) * n[0] + std::real(e2*(v01+v10)/two) * mx[0] + std::real(e2*(v10-v01)/twoi) * my[0] + std::real(e2*(v00-v11)/two) * mz[0];// vtxc is used the calculation of the total energy(Ts more specifically), because abacus doesn't directly programme the kinetic operator and instead uses the sum of occupied orbital energy
                    }
                }   
            }
            else
            { // gga
                std::vector<double> n(nrxx);
                std::vector<double> mx(nrxx);
                std::vector<double> my(nrxx);
                std::vector<double> mz(nrxx);
                std::vector<double> gradx_n(nrxx);
                std::vector<double> grady_n(nrxx);
                std::vector<double> gradz_n(nrxx);
                std::vector<double> gradx_mx(nrxx);
                std::vector<double> grady_mx(nrxx);
                std::vector<double> gradz_mx(nrxx);
                std::vector<double> gradx_my(nrxx);
                std::vector<double> grady_my(nrxx);
                std::vector<double> gradz_my(nrxx);
                std::vector<double> gradx_mz(nrxx);
                std::vector<double> grady_mz(nrxx);
                std::vector<double> gradz_mz(nrxx);
                std::vector<double> grad2xx_n(nrxx);
                std::vector<double> grad2yy_n(nrxx);
                std::vector<double> grad2zz_n(nrxx);
                std::vector<double> grad2xx_mx(nrxx);
                std::vector<double> grad2yy_mx(nrxx);
                std::vector<double> grad2zz_mx(nrxx);
                std::vector<double> grad2xx_my(nrxx);
                std::vector<double> grad2yy_my(nrxx);
                std::vector<double> grad2zz_my(nrxx);
                std::vector<double> grad2xx_mz(nrxx);
                std::vector<double> grad2yy_mz(nrxx);
                std::vector<double> grad2zz_mz(nrxx);
                std::vector<double> grad2xy_n(nrxx);
                std::vector<double> grad2yz_n(nrxx);
                std::vector<double> grad2xz_n(nrxx);
                std::vector<double> grad2xy_mx(nrxx);
                std::vector<double> grad2yz_mx(nrxx);
                std::vector<double> grad2xz_mx(nrxx);
                std::vector<double> grad2xy_my(nrxx);
                std::vector<double> grad2yz_my(nrxx);
                std::vector<double> grad2xz_my(nrxx);
                std::vector<double> grad2xy_mz(nrxx);
                std::vector<double> grad2yz_mz(nrxx);
                std::vector<double> grad2xz_mz(nrxx);

                std::vector<double> grad3xxx_n(nrxx);
                std::vector<double> grad3xxy_n(nrxx);
                std::vector<double> grad3xxz_n(nrxx);
                std::vector<double> grad3xyy_n(nrxx);
                std::vector<double> grad3xyz_n(nrxx);
                std::vector<double> grad3xzz_n(nrxx);
                std::vector<double> grad3yyy_n(nrxx);
                std::vector<double> grad3yyz_n(nrxx);
                std::vector<double> grad3yzz_n(nrxx);
                std::vector<double> grad3zzz_n(nrxx);

                std::vector<double> grad3xxx_mx(nrxx);
                std::vector<double> grad3xxy_mx(nrxx);
                std::vector<double> grad3xxz_mx(nrxx);
                std::vector<double> grad3xyy_mx(nrxx);
                std::vector<double> grad3xyz_mx(nrxx);
                std::vector<double> grad3xzz_mx(nrxx);
                std::vector<double> grad3yyy_mx(nrxx);
                std::vector<double> grad3yyz_mx(nrxx);
                std::vector<double> grad3yzz_mx(nrxx);
                std::vector<double> grad3zzz_mx(nrxx);

                std::vector<double> grad3xxx_my(nrxx);
                std::vector<double> grad3xxy_my(nrxx);
                std::vector<double> grad3xxz_my(nrxx);
                std::vector<double> grad3xyy_my(nrxx);
                std::vector<double> grad3xyz_my(nrxx);
                std::vector<double> grad3xzz_my(nrxx);
                std::vector<double> grad3yyy_my(nrxx);
                std::vector<double> grad3yyz_my(nrxx);
                std::vector<double> grad3yzz_my(nrxx);
                std::vector<double> grad3zzz_my(nrxx);

                std::vector<double> grad3xxx_mz(nrxx);
                std::vector<double> grad3xxy_mz(nrxx);
                std::vector<double> grad3xxz_mz(nrxx);
                std::vector<double> grad3xyy_mz(nrxx);
                std::vector<double> grad3xyz_mz(nrxx);
                std::vector<double> grad3xzz_mz(nrxx);
                std::vector<double> grad3yyy_mz(nrxx);
                std::vector<double> grad3yyz_mz(nrxx);
                std::vector<double> grad3yzz_mz(nrxx);
                std::vector<double> grad3zzz_mz(nrxx);

                std::vector<double> grad4xxxx_n(nrxx);
                std::vector<double> grad4xxxy_n(nrxx);
                std::vector<double> grad4xxxz_n(nrxx);
                std::vector<double> grad4xxyy_n(nrxx);
                std::vector<double> grad4xxyz_n(nrxx);
                std::vector<double> grad4xxzz_n(nrxx);
                std::vector<double> grad4xyyy_n(nrxx);
                std::vector<double> grad4xyyz_n(nrxx);
                std::vector<double> grad4xyzz_n(nrxx);
                std::vector<double> grad4xzzz_n(nrxx);
                std::vector<double> grad4yyyy_n(nrxx);
                std::vector<double> grad4yyyz_n(nrxx);
                std::vector<double> grad4yyzz_n(nrxx);
                std::vector<double> grad4yzzz_n(nrxx);
                std::vector<double> grad4zzzz_n(nrxx);

                std::vector<double> grad4xxxx_mx(nrxx);
                std::vector<double> grad4xxxy_mx(nrxx);
                std::vector<double> grad4xxxz_mx(nrxx);
                std::vector<double> grad4xxyy_mx(nrxx);
                std::vector<double> grad4xxyz_mx(nrxx);
                std::vector<double> grad4xxzz_mx(nrxx);
                std::vector<double> grad4xyyy_mx(nrxx);
                std::vector<double> grad4xyyz_mx(nrxx);
                std::vector<double> grad4xyzz_mx(nrxx);
                std::vector<double> grad4xzzz_mx(nrxx);
                std::vector<double> grad4yyyy_mx(nrxx);
                std::vector<double> grad4yyyz_mx(nrxx);
                std::vector<double> grad4yyzz_mx(nrxx);
                std::vector<double> grad4yzzz_mx(nrxx);
                std::vector<double> grad4zzzz_mx(nrxx);

                std::vector<double> grad4xxxx_my(nrxx);
                std::vector<double> grad4xxxy_my(nrxx);
                std::vector<double> grad4xxxz_my(nrxx);
                std::vector<double> grad4xxyy_my(nrxx);
                std::vector<double> grad4xxyz_my(nrxx);
                std::vector<double> grad4xxzz_my(nrxx);
                std::vector<double> grad4xyyy_my(nrxx);
                std::vector<double> grad4xyyz_my(nrxx);
                std::vector<double> grad4xyzz_my(nrxx);
                std::vector<double> grad4xzzz_my(nrxx);
                std::vector<double> grad4yyyy_my(nrxx);
                std::vector<double> grad4yyyz_my(nrxx);
                std::vector<double> grad4yyzz_my(nrxx);
                std::vector<double> grad4yzzz_my(nrxx);
                std::vector<double> grad4zzzz_my(nrxx);

                std::vector<double> grad4xxxx_mz(nrxx);
                std::vector<double> grad4xxxy_mz(nrxx);
                std::vector<double> grad4xxxz_mz(nrxx);
                std::vector<double> grad4xxyy_mz(nrxx);
                std::vector<double> grad4xxyz_mz(nrxx);
                std::vector<double> grad4xxzz_mz(nrxx);
                std::vector<double> grad4xyyy_mz(nrxx);
                std::vector<double> grad4xyyz_mz(nrxx);
                std::vector<double> grad4xyzz_mz(nrxx);
                std::vector<double> grad4xzzz_mz(nrxx);
                std::vector<double> grad4yyyy_mz(nrxx);
                std::vector<double> grad4yyyz_mz(nrxx);
                std::vector<double> grad4yyzz_mz(nrxx);
                std::vector<double> grad4yzzz_mz(nrxx);
                std::vector<double> grad4zzzz_mz(nrxx);

                std::vector<std::vector<ModuleBase::Vector3<double>>> gdr;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    n[ir] = chr->rho[0][ir] + chr->rho_core[ir];
                    mx[ir] = chr->rho[1][ir];
                    my[ir] = chr->rho[2][ir];
                    mz[ir] = chr->rho[3][ir];
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, n, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    gradx_n[ir] = gdr[0][ir].x;
                    grady_n[ir] = gdr[0][ir].y;
                    gradz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, mx, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    gradx_mx[ir] = gdr[0][ir].x;
                    grady_mx[ir] = gdr[0][ir].y;
                    gradz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, my, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    gradx_my[ir] = gdr[0][ir].x;
                    grady_my[ir] = gdr[0][ir].y;
                    gradz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, mz, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    gradx_mz[ir] = gdr[0][ir].x;
                    grady_mz[ir] = gdr[0][ir].y;
                    gradz_mz[ir] = gdr[0][ir].z;
                }

                // first derivatives of n -> second derivatives
                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, gradx_n, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad2xx_n[ir] = gdr[0][ir].x;
                    grad2xy_n[ir] = gdr[0][ir].y;
                    grad2xz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grady_n, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    // x-component unused for grad2yy_n, so we only take y,z here
                    grad2yy_n[ir] = gdr[0][ir].y;
                    grad2yz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, gradz_n, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    // x,y components unused for grad2zz_n
                    grad2zz_n[ir] = gdr[0][ir].z;
                }

                // first derivatives of mx -> second derivatives
                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, gradx_mx, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad2xx_mx[ir] = gdr[0][ir].x;
                    grad2xy_mx[ir] = gdr[0][ir].y;
                    grad2xz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grady_mx, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad2yy_mx[ir] = gdr[0][ir].y;
                    grad2yz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, gradz_mx, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad2zz_mx[ir] = gdr[0][ir].z;
                }

                // first derivatives of my -> second derivatives
                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, gradx_my, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad2xx_my[ir] = gdr[0][ir].x;
                    grad2xy_my[ir] = gdr[0][ir].y;
                    grad2xz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grady_my, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad2yy_my[ir] = gdr[0][ir].y;
                    grad2yz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, gradz_my, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad2zz_my[ir] = gdr[0][ir].z;
                }

                // first derivatives of mz -> second derivatives
                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, gradx_mz, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad2xx_mz[ir] = gdr[0][ir].x;
                    grad2xy_mz[ir] = gdr[0][ir].y;
                    grad2xz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grady_mz, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad2yy_mz[ir] = gdr[0][ir].y;
                    grad2yz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, gradz_mz, tpiba, chr);
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024)
                #endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad2zz_mz[ir] = gdr[0][ir].z;
                }
                

                // === Second to Third Derivatives of n ===
                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2xx_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3xxx_n[ir] = gdr[0][ir].x;
                    grad3xxy_n[ir] = gdr[0][ir].y;
                    grad3xxz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2xy_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3xyy_n[ir] = gdr[0][ir].y;
                    grad3xyz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2xz_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3xzz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2yy_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3yyy_n[ir] = gdr[0][ir].y;
                    grad3yyz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2yz_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3yzz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2zz_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3zzz_n[ir] = gdr[0][ir].z;
                }

                // === Second to Third Derivatives of mx ===
                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2xx_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3xxx_mx[ir] = gdr[0][ir].x;
                    grad3xxy_mx[ir] = gdr[0][ir].y;
                    grad3xxz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2xy_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3xyy_mx[ir] = gdr[0][ir].y;
                    grad3xyz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2xz_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3xzz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2yy_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3yyy_mx[ir] = gdr[0][ir].y;
                    grad3yyz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2yz_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3yzz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2zz_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3zzz_mx[ir] = gdr[0][ir].z;
                }

                // === Second to Third Derivatives of my ===
                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2xx_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3xxx_my[ir] = gdr[0][ir].x;
                    grad3xxy_my[ir] = gdr[0][ir].y;
                    grad3xxz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2xy_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3xyy_my[ir] = gdr[0][ir].y;
                    grad3xyz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2xz_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3xzz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2yy_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3yyy_my[ir] = gdr[0][ir].y;
                    grad3yyz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2yz_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3yzz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2zz_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3zzz_my[ir] = gdr[0][ir].z;
                }

                // === Second to Third Derivatives of mz ===
                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2xx_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3xxx_mz[ir] = gdr[0][ir].x;
                    grad3xxy_mz[ir] = gdr[0][ir].y;
                    grad3xxz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2xy_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3xyy_mz[ir] = gdr[0][ir].y;
                    grad3xyz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2xz_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3xzz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2yy_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3yyy_mz[ir] = gdr[0][ir].y;
                    grad3yyz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2yz_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3yzz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad2zz_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad3zzz_mz[ir] = gdr[0][ir].z;
                }

                // === Third to Fourth Derivatives of n ===
                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xxx_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xxxx_n[ir] = gdr[0][ir].x;
                    grad4xxxy_n[ir] = gdr[0][ir].y;
                    grad4xxxz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xxy_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xxyy_n[ir] = gdr[0][ir].y;
                    grad4xxyz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xxz_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xxzz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xyy_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xyyy_n[ir] = gdr[0][ir].y;
                    grad4xyyz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xyz_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xyzz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xzz_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xzzz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3yyy_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4yyyy_n[ir] = gdr[0][ir].y;
                    grad4yyyz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3yyz_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4yyzz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3yzz_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4yzzz_n[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3zzz_n, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4zzzz_n[ir] = gdr[0][ir].z;
                }

                // === Third to Fourth Derivatives of mx ===
                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xxx_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xxxx_mx[ir] = gdr[0][ir].x;
                    grad4xxxy_mx[ir] = gdr[0][ir].y;
                    grad4xxxz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xxy_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xxyy_mx[ir] = gdr[0][ir].y;
                    grad4xxyz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xxz_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xxzz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xyy_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xyyy_mx[ir] = gdr[0][ir].y;
                    grad4xyyz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xyz_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xyzz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xzz_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xzzz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3yyy_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4yyyy_mx[ir] = gdr[0][ir].y;
                    grad4yyyz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3yyz_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4yyzz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3yzz_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4yzzz_mx[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3zzz_mx, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4zzzz_mx[ir] = gdr[0][ir].z;
                }

                // === Third to Fourth Derivatives of my ===
                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xxx_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xxxx_my[ir] = gdr[0][ir].x;
                    grad4xxxy_my[ir] = gdr[0][ir].y;
                    grad4xxxz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xxy_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xxyy_my[ir] = gdr[0][ir].y;
                    grad4xxyz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xxz_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xxzz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xyy_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xyyy_my[ir] = gdr[0][ir].y;
                    grad4xyyz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xyz_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xyzz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xzz_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xzzz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3yyy_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4yyyy_my[ir] = gdr[0][ir].y;
                    grad4yyyz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3yyz_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4yyzz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3yzz_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4yzzz_my[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3zzz_my, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4zzzz_my[ir] = gdr[0][ir].z;
                }

                // === Third to Fourth Derivatives of mz ===
                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xxx_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xxxx_mz[ir] = gdr[0][ir].x;
                    grad4xxxy_mz[ir] = gdr[0][ir].y;
                    grad4xxxz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xxy_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xxyy_mz[ir] = gdr[0][ir].y;
                    grad4xxyz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xxz_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xxzz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xyy_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xyyy_mz[ir] = gdr[0][ir].y;
                    grad4xyyz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xyz_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xyzz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3xzz_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4xzzz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3yyy_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4yyyy_mz[ir] = gdr[0][ir].y;
                    grad4yyyz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3yyz_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4yyzz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3yzz_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4yzzz_mz[ir] = gdr[0][ir].z;
                }

                gdr = XC_Functional_Libxc::cal_gdr(1, nrxx, grad3zzz_mz, tpiba, chr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (int ir = 0; ir < nrxx; ++ir) {
                    grad4zzzz_mz[ir] = gdr[0][ir].z;
                }



                std::vector<double> E_MC;
                std::vector<Matrix2x2> V_MC;
                
                for (const int& id : func_id) {
                    std::tie(E_MC, V_MC) = NCLibxc::gga_mc(
    id, 
    n, mx, my, mz,
    gradx_n, grady_n, gradz_n,
    gradx_mx, grady_mx, gradz_mx,
    gradx_my, grady_my, gradz_my,
    gradx_mz, grady_mz, gradz_mz,
    grad2xx_n, grad2yy_n, grad2zz_n,
    grad2xy_n, grad2yz_n, grad2xz_n,
    grad2xx_mx, grad2yy_mx, grad2zz_mx,
    grad2xy_mx, grad2yz_mx, grad2xz_mx,
    grad2xx_my, grad2yy_my, grad2zz_my,
    grad2xy_my, grad2yz_my, grad2xz_my,
    grad2xx_mz, grad2yy_mz, grad2zz_mz,
    grad2xy_mz, grad2yz_mz, grad2xz_mz,
    grad3xxx_n, grad3xxy_n, grad3xxz_n,
    grad3xyy_n, grad3xyz_n, grad3xzz_n,
    grad3yyy_n, grad3yyz_n, grad3yzz_n,
    grad3zzz_n,
    grad3xxx_mx, grad3xxy_mx, grad3xxz_mx,
    grad3xyy_mx, grad3xyz_mx, grad3xzz_mx,
    grad3yyy_mx, grad3yyz_mx, grad3yzz_mx,
    grad3zzz_mx,
    grad3xxx_my, grad3xxy_my, grad3xxz_my,
    grad3xyy_my, grad3xyz_my, grad3xzz_my,
    grad3yyy_my, grad3yyz_my, grad3yzz_my,
    grad3zzz_my,
    grad3xxx_mz, grad3xxy_mz, grad3xxz_mz,
    grad3xyy_mz, grad3xyz_mz, grad3xzz_mz,
    grad3yyy_mz, grad3yyz_mz, grad3yzz_mz,
    grad3zzz_mz,
    grad4xxxx_n, grad4xxxy_n, grad4xxxz_n,
    grad4xxyy_n, grad4xxyz_n, grad4xxzz_n,
    grad4xyyy_n, grad4xyyz_n, grad4xyzz_n,
    grad4xzzz_n, grad4yyyy_n, grad4yyyz_n,
    grad4yyzz_n, grad4yzzz_n, grad4zzzz_n,
    grad4xxxx_mx, grad4xxxy_mx, grad4xxxz_mx,
    grad4xxyy_mx, grad4xxyz_mx, grad4xxzz_mx,
    grad4xyyy_mx, grad4xyyz_mx, grad4xyzz_mx,
    grad4xzzz_mx, grad4yyyy_mx, grad4yyyz_mx,
    grad4yyzz_mx, grad4yzzz_mx, grad4zzzz_mx,
    grad4xxxx_my, grad4xxxy_my, grad4xxxz_my,
    grad4xxyy_my, grad4xxyz_my, grad4xxzz_my,
    grad4xyyy_my, grad4xyyz_my, grad4xyzz_my,
    grad4xzzz_my, grad4yyyy_my, grad4yyyz_my,
    grad4yyzz_my, grad4yzzz_my, grad4zzzz_my,
    grad4xxxx_mz, grad4xxxy_mz, grad4xxxz_mz,
    grad4xxyy_mz, grad4xxyz_mz, grad4xxzz_mz,
    grad4xyyy_mz, grad4xyyz_mz, grad4xyzz_mz,
    grad4xzzz_mz, grad4yyyy_mz, grad4yyyz_mz,
    grad4yyzz_mz, grad4yzzz_mz, grad4zzzz_mz
);

             if (PARAM.inp.xc_torque){
                std::vector<double> torque_tmp = NCLibxc::gga_torque(
                    mx,my,mz,V_MC
                );
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static, 1024) 
                #endif
                    for (int ir = 0; ir < 3 * nrxx; ++ir) {
                        torque[ir] += torque_tmp[ir];
                    }
            }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024) reduction(+:etxc) reduction(+:vtxc)
#endif
                    for (int ir = 0; ir < nrxx; ++ir) {
                        double exc = e2 * E_MC[ir];
                        auto v00 = V_MC[ir][0][0];
                        auto v01 = V_MC[ir][0][1];
                        auto v10 = V_MC[ir][1][0];
                        auto v11 = V_MC[ir][1][1];
                        v_nspin4(0, ir) += std::real(e2 * (v00 + v11) / two);
                        v_nspin4(1, ir) += std::real(e2 * (v01 + v10) / two);
                        v_nspin4(2, ir) += std::real(e2 * (v10 - v01) / twoi);
                        v_nspin4(3, ir) += std::real(e2 * (v00 - v11) / two);
                        etxc += exc * n[ir];
                        vtxc += std::real(e2 * (v00 + v11) / two) * n[ir] + std::real(e2 * (v01 + v10) / two) * mx[ir] + std::real(e2 * (v10 - v01) / twoi) * my[ir] + std::real(e2 * (v00 - v11) / two) * mz[ir];
                    }
                }
/*
                //collinear limit tests
                std::vector<double> e(nrxx);
                // Prepare inputs for postlibxc_gga
                size_t size = n.size();
                std::vector<double> rho0(size), rho1(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (size_t i = 0; i < size; ++i) {
                    rho0[i] = (n[i] + mz[i]) / 2.0;
                    rho1[i] = (n[i] - mz[i]) / 2.0;
                }

                // Gradient calculations
                std::vector<double> gradx_rho0(size), grady_rho0(size), gradz_rho0(size);
                std::vector<double> gradx_rho1(size), grady_rho1(size), gradz_rho1(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (size_t i = 0; i < size; ++i) {
                    gradx_rho0[i] = (gradx_n[i] + gradx_mz[i]) / 2.0;
                    grady_rho0[i] = (grady_n[i] + grady_mz[i]) / 2.0;
                    gradz_rho0[i] = (gradz_n[i] + gradz_mz[i]) / 2.0;

                    gradx_rho1[i] = (gradx_n[i] - gradx_mz[i]) / 2.0;
                    grady_rho1[i] = (grady_n[i] - grady_mz[i]) / 2.0;
                    gradz_rho1[i] = (gradz_n[i] - gradz_mz[i]) / 2.0;
                }

                // Second derivatives
                std::vector<double> grad2xx_rho0(size), grad2yy_rho0(size), grad2zz_rho0(size);
                std::vector<double> grad2xy_rho0(size), grad2yz_rho0(size), grad2xz_rho0(size);
                std::vector<double> grad2xx_rho1(size), grad2yy_rho1(size), grad2zz_rho1(size);
                std::vector<double> grad2xy_rho1(size), grad2yz_rho1(size), grad2xz_rho1(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                for (size_t i = 0; i < size; ++i) {
                    grad2xx_rho0[i] = (grad2xx_n[i] + grad2xx_mz[i]) / 2.0;
                    grad2yy_rho0[i] = (grad2yy_n[i] + grad2yy_mz[i]) / 2.0;
                    grad2zz_rho0[i] = (grad2zz_n[i] + grad2zz_mz[i]) / 2.0;
                    grad2xy_rho0[i] = (grad2xy_n[i] + grad2xy_mz[i]) / 2.0;
                    grad2yz_rho0[i] = (grad2yz_n[i] + grad2yz_mz[i]) / 2.0;
                    grad2xz_rho0[i] = (grad2xz_n[i] + grad2xz_mz[i]) / 2.0;

                    grad2xx_rho1[i] = (grad2xx_n[i] - grad2xx_mz[i]) / 2.0;
                    grad2yy_rho1[i] = (grad2yy_n[i] - grad2yy_mz[i]) / 2.0;
                    grad2zz_rho1[i] = (grad2zz_n[i] - grad2zz_mz[i]) / 2.0;
                    grad2xy_rho1[i] = (grad2xy_n[i] - grad2xy_mz[i]) / 2.0;
                    grad2yz_rho1[i] = (grad2yz_n[i] - grad2yz_mz[i]) / 2.0;
                    grad2xz_rho1[i] = (grad2xz_n[i] - grad2xz_mz[i]) / 2.0;
                }

                // Third derivatives
                std::vector<double> grad3xxx_rho0(size), grad3xxy_rho0(size), grad3xxz_rho0(size);
                std::vector<double> grad3xyy_rho0(size), grad3xyz_rho0(size), grad3xzz_rho0(size);
                std::vector<double> grad3yyy_rho0(size), grad3yyz_rho0(size), grad3yzz_rho0(size);
                std::vector<double> grad3zzz_rho0(size);

                std::vector<double> grad3xxx_rho1(size), grad3xxy_rho1(size), grad3xxz_rho1(size);
                std::vector<double> grad3xyy_rho1(size), grad3xyz_rho1(size), grad3xzz_rho1(size);
                std::vector<double> grad3yyy_rho1(size), grad3yyz_rho1(size), grad3yzz_rho1(size);
                std::vector<double> grad3zzz_rho1(size);

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif

                                for (size_t i = 0; i < size; ++i) {
                                    // rho0 third derivatives
                                    grad3xxx_rho0[i] = (grad3xxx_n[i] + grad3xxx_mz[i]) / 2.0;
                                    grad3xxy_rho0[i] = (grad3xxy_n[i] + grad3xxy_mz[i]) / 2.0;
                                    grad3xxz_rho0[i] = (grad3xxz_n[i] + grad3xxz_mz[i]) / 2.0;
                                    grad3xyy_rho0[i] = (grad3xyy_n[i] + grad3xyy_mz[i]) / 2.0;
                                    grad3xyz_rho0[i] = (grad3xyz_n[i] + grad3xyz_mz[i]) / 2.0;
                                    grad3xzz_rho0[i] = (grad3xzz_n[i] + grad3xzz_mz[i]) / 2.0;
                                    grad3yyy_rho0[i] = (grad3yyy_n[i] + grad3yyy_mz[i]) / 2.0;
                                    grad3yyz_rho0[i] = (grad3yyz_n[i] + grad3yyz_mz[i]) / 2.0;
                                    grad3yzz_rho0[i] = (grad3yzz_n[i] + grad3yzz_mz[i]) / 2.0;
                                    grad3zzz_rho0[i] = (grad3zzz_n[i] + grad3zzz_mz[i]) / 2.0;

                                    // rho1 third derivatives
                                    grad3xxx_rho1[i] = (grad3xxx_n[i] - grad3xxx_mz[i]) / 2.0;
                                    grad3xxy_rho1[i] = (grad3xxy_n[i] - grad3xxy_mz[i]) / 2.0;
                                    grad3xxz_rho1[i] = (grad3xxz_n[i] - grad3xxz_mz[i]) / 2.0;
                                    grad3xyy_rho1[i] = (grad3xyy_n[i] - grad3xyy_mz[i]) / 2.0;
                                    grad3xyz_rho1[i] = (grad3xyz_n[i] - grad3xyz_mz[i]) / 2.0;
                                    grad3xzz_rho1[i] = (grad3xzz_n[i] - grad3xzz_mz[i]) / 2.0;
                                    grad3yyy_rho1[i] = (grad3yyy_n[i] - grad3yyy_mz[i]) / 2.0;
                                    grad3yyz_rho1[i] = (grad3yyz_n[i] - grad3yyz_mz[i]) / 2.0;
                                    grad3yzz_rho1[i] = (grad3yzz_n[i] - grad3yzz_mz[i]) / 2.0;
                                    grad3zzz_rho1[i] = (grad3zzz_n[i] - grad3zzz_mz[i]) / 2.0;
                                }

                // Fourth derivatives
                std::vector<double> grad4xxxx_rho0(size), grad4xxxy_rho0(size), grad4xxxz_rho0(size);
                std::vector<double> grad4xxyy_rho0(size), grad4xxyz_rho0(size), grad4xxzz_rho0(size);
                std::vector<double> grad4xyyy_rho0(size), grad4xyyz_rho0(size), grad4xyzz_rho0(size);
                std::vector<double> grad4xzzz_rho0(size), grad4yyyy_rho0(size), grad4yyyz_rho0(size);
                std::vector<double> grad4yyzz_rho0(size), grad4yzzz_rho0(size), grad4zzzz_rho0(size);

                std::vector<double> grad4xxxx_rho1(size), grad4xxxy_rho1(size), grad4xxxz_rho1(size);
                std::vector<double> grad4xxyy_rho1(size), grad4xxyz_rho1(size), grad4xxzz_rho1(size);
                std::vector<double> grad4xyyy_rho1(size), grad4xyyz_rho1(size), grad4xyzz_rho1(size);
                std::vector<double> grad4xzzz_rho1(size), grad4yyyy_rho1(size), grad4yyyz_rho1(size);
                std::vector<double> grad4yyzz_rho1(size), grad4yzzz_rho1(size), grad4zzzz_rho1(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
                // Compute fourth derivatives for rho0 and rho1
for (size_t i = 0; i < size; ++i) {
    // rho0 fourth derivatives
    grad4xxxx_rho0[i] = (grad4xxxx_n[i] + grad4xxxx_mz[i]) / 2.0;
    grad4xxxy_rho0[i] = (grad4xxxy_n[i] + grad4xxxy_mz[i]) / 2.0;
    grad4xxxz_rho0[i] = (grad4xxxz_n[i] + grad4xxxz_mz[i]) / 2.0;
    grad4xxyy_rho0[i] = (grad4xxyy_n[i] + grad4xxyy_mz[i]) / 2.0;
    grad4xxyz_rho0[i] = (grad4xxyz_n[i] + grad4xxyz_mz[i]) / 2.0;
    grad4xxzz_rho0[i] = (grad4xxzz_n[i] + grad4xxzz_mz[i]) / 2.0;
    grad4xyyy_rho0[i] = (grad4xyyy_n[i] + grad4xyyy_mz[i]) / 2.0;
    grad4xyyz_rho0[i] = (grad4xyyz_n[i] + grad4xyyz_mz[i]) / 2.0;
    grad4xyzz_rho0[i] = (grad4xyzz_n[i] + grad4xyzz_mz[i]) / 2.0;
    grad4xzzz_rho0[i] = (grad4xzzz_n[i] + grad4xzzz_mz[i]) / 2.0;
    grad4yyyy_rho0[i] = (grad4yyyy_n[i] + grad4yyyy_mz[i]) / 2.0;
    grad4yyyz_rho0[i] = (grad4yyyz_n[i] + grad4yyyz_mz[i]) / 2.0;
    grad4yyzz_rho0[i] = (grad4yyzz_n[i] + grad4yyzz_mz[i]) / 2.0;
    grad4yzzz_rho0[i] = (grad4yzzz_n[i] + grad4yzzz_mz[i]) / 2.0;
    grad4zzzz_rho0[i] = (grad4zzzz_n[i] + grad4zzzz_mz[i]) / 2.0;

    // rho1 fourth derivatives
    grad4xxxx_rho1[i] = (grad4xxxx_n[i] - grad4xxxx_mz[i]) / 2.0;
    grad4xxxy_rho1[i] = (grad4xxxy_n[i] - grad4xxxy_mz[i]) / 2.0;
    grad4xxxz_rho1[i] = (grad4xxxz_n[i] - grad4xxxz_mz[i]) / 2.0;
    grad4xxyy_rho1[i] = (grad4xxyy_n[i] - grad4xxyy_mz[i]) / 2.0;
    grad4xxyz_rho1[i] = (grad4xxyz_n[i] - grad4xxyz_mz[i]) / 2.0;
    grad4xxzz_rho1[i] = (grad4xxzz_n[i] - grad4xxzz_mz[i]) / 2.0;
    grad4xyyy_rho1[i] = (grad4xyyy_n[i] - grad4xyyy_mz[i]) / 2.0;
    grad4xyyz_rho1[i] = (grad4xyyz_n[i] - grad4xyyz_mz[i]) / 2.0;
    grad4xyzz_rho1[i] = (grad4xyzz_n[i] - grad4xyzz_mz[i]) / 2.0;
    grad4xzzz_rho1[i] = (grad4xzzz_n[i] - grad4xzzz_mz[i]) / 2.0;
    grad4yyyy_rho1[i] = (grad4yyyy_n[i] - grad4yyyy_mz[i]) / 2.0;
    grad4yyyz_rho1[i] = (grad4yyyz_n[i] - grad4yyyz_mz[i]) / 2.0;
    grad4yyzz_rho1[i] = (grad4yyzz_n[i] - grad4yyzz_mz[i]) / 2.0;
    grad4yzzz_rho1[i] = (grad4yzzz_n[i] - grad4yzzz_mz[i]) / 2.0;
    grad4zzzz_rho1[i] = (grad4zzzz_n[i] - grad4zzzz_mz[i]) / 2.0;
}

                // Call postlibxc_gga
                std::vector<double> e_col, v1_col, v2_col, f1_col, f2_col, f3_col;
                

                for (const int &id : func_id) {
                    NCLibxc::postlibxc_gga(
                   id,
    rho0,
    rho1,
    gradx_rho0,
    grady_rho0,
    gradz_rho0,
    gradx_rho1,
    grady_rho1,
    gradz_rho1,
    grad2xx_rho0,
    grad2yy_rho0,
    grad2zz_rho0,
    grad2xy_rho0,
    grad2yz_rho0,
    grad2xz_rho0,
    grad2xx_rho1,
    grad2yy_rho1,
    grad2zz_rho1,
    grad2xy_rho1,
    grad2yz_rho1,
    grad2xz_rho1,
    grad3xxx_rho0,
    grad3xxy_rho0,
    grad3xxz_rho0,
    grad3xyy_rho0,
    grad3xyz_rho0,
    grad3xzz_rho0,
    grad3yyy_rho0,
    grad3yyz_rho0,
    grad3yzz_rho0,
    grad3zzz_rho0,
    grad3xxx_rho1,
    grad3xxy_rho1,
    grad3xxz_rho1,
    grad3xyy_rho1,
    grad3xyz_rho1,
    grad3xzz_rho1,
    grad3yyy_rho1,
    grad3yyz_rho1,
    grad3yzz_rho1,
    grad3zzz_rho1,
    grad4xxxx_rho0,
    grad4xxxy_rho0,
    grad4xxxz_rho0,
    grad4xxyy_rho0,
    grad4xxyz_rho0,
    grad4xxzz_rho0,
    grad4xyyy_rho0,
    grad4xyyz_rho0,
    grad4xyzz_rho0,
    grad4xzzz_rho0,
    grad4yyyy_rho0,
    grad4yyyz_rho0,
    grad4yyzz_rho0,
    grad4yzzz_rho0,
    grad4zzzz_rho0,
    grad4xxxx_rho1,
    grad4xxxy_rho1,
    grad4xxxz_rho1,
    grad4xxyy_rho1,
    grad4xxyz_rho1,
    grad4xxzz_rho1,
    grad4xyyy_rho1,
    grad4xyyz_rho1,
    grad4xyzz_rho1,
    grad4xzzz_rho1,
    grad4yyyy_rho1,
    grad4yyyz_rho1,
    grad4yyzz_rho1,
    grad4yzzz_rho1,
    grad4zzzz_rho1,
                    e_col, v1_col, v2_col, f1_col, f2_col, f3_col
                    );
                    rho0_0 = rho0[0];
                    rho1_0 = rho1[0];
                    v1_col0 += v1_col[0];
                    v2_col0 += v2_col[0];
                    for (int ir = 0; ir < nrxx; ++ir) {
                        etxc_col += n[ir] * e_col[ir]*e2;
                        vtxc_col+= e2*(rho0[ir] * v1_col[ir] + rho1[ir] * v2_col[ir]);
                    }
                }*/
                
            }
        }
    }

    //-------------------------------------------------
    // for MPI, reduce the exchange-correlation energy
    //-------------------------------------------------
    #ifdef __MPI
    Parallel_Reduce::reduce_pool(etxc);
    Parallel_Reduce::reduce_pool(vtxc);
    //Parallel_Reduce::reduce_pool(etxc_col);
    //Parallel_Reduce::reduce_pool(vtxc_col);
    #endif
    if(PARAM.inp.xc_torque){
        NCLibxc::print_torque(torque);
    }
    
/*
    std::cout << std::fixed << std::setprecision(8) << "etxc: " << etxc << std::endl;
    std::cout << std::fixed << std::setprecision(8) << "etxc_col: " << etxc_col << std::endl;
    std::cout << std::fixed << std::setprecision(8) << "vtxc: " << vtxc << std::endl;
    std::cout << std::fixed << std::setprecision(8) << "vtxc_col: " << vtxc_col << std::endl; 

    // Debug prints for first point values
    std::cout << std::fixed << std::setprecision(8) << "n[0]: " << chr->rho[0][0] << std::endl;
    std::cout << std::fixed << std::setprecision(8) << "mx[0]: " << chr->rho[1][0] << std::endl;
    std::cout << std::fixed << std::setprecision(8) << "my[0]: " << chr->rho[2][0] << std::endl;
    std::cout << std::fixed << std::setprecision(8) << "mz[0]: " << chr->rho[3][0] << std::endl;
    std::cout << std::fixed << std::setprecision(8) << "rho0[0]: " << rho0_0 << std::endl;
    std::cout << std::fixed << std::setprecision(8) << "rho1[0]: " << rho1_0 << std::endl;
    std::cout << std::fixed << std::setprecision(8) << "v1_col[0]:"  << v1_col0 << std::endl;
    std::cout << std::fixed << std::setprecision(8) << "v2_col[0]:"  << v2_col0 << std::endl;
    std::cout << std::fixed << std::setprecision(8) << "nspin1[0]: " << v_nspin4(0,0) << std::endl;
    std::cout << std::fixed << std::setprecision(8) << "nspin2[0]: " << v_nspin4(1,0) << std::endl;
    std::cout << std::fixed << std::setprecision(8) << "nspin3[0]: " << v_nspin4(2,0) << std::endl;
    std::cout << std::fixed << std::setprecision(8) << "nspin4[0]: " << v_nspin4(3,0) << std::endl;
*/
    etxc *= omega / chr->rhopw->nxyz;
    vtxc *= omega / chr->rhopw->nxyz;

    XC_Functional_Libxc::finish_func(funcs);

    ModuleBase::timer::tick("XC_Functional_Libxc","v_xc_libxc");
    if(use_v_nspin4)
    {
        return std::make_tuple(etxc, vtxc, std::move(v_nspin4));
    }
    else
    {
        return std::make_tuple(etxc, vtxc, std::move(v));
    }
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
            if ( rho[ir*2]<rho_th || sqrt(std::abs(sigma[ir*3]))<grho_th || std::abs(kin_r[ir*2])<tau_th) {
                sgn[ir*2] = 0.0;
}
            if ( rho[ir*2+1]<rho_th || sqrt(std::abs(sigma[ir*3+2]))<grho_th || std::abs(kin_r[ir*2+1])<tau_th) {
                sgn[ir*2+1] = 0.0;
}
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
