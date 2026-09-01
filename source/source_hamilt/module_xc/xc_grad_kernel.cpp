// Second stage of the gradient-correction pipeline: evaluate the
// GGA/meta-GGA kernel on every real-space grid point (own implementation
// or libxc) and accumulate vtxcgc, etxcgc, stress_gga and v.
// The stage functions and their shared parameter/buffer structs are
// declared in xc_grad_internal.h; the public entry point
// XC_Functional::gradcorr lives in xc_grad.cpp.

#include "xc_functional.h"
#include "xc_grad_internal.h"
#include "source_base/constants.h"

#ifdef __LIBXC
#include "libxc_abacus.h"
#endif

void gradcorr_xc_kernel(const GradCorrParams& params,
                        double& vtxcgc,
                        double& etxcgc,
                        std::vector<double>& stress_gga,
                        ModuleBase::matrix& v,
                        GradCorrBuffers& buf)
{
    const Charge* const chr = params.chr;
    ModulePW::PW_Basis* rhopw = params.rhopw;
    const int nspin = params.nspin;
    const int nspin0 = params.nspin0;
    const double fac = params.fac;
    const bool is_stress = params.is_stress;
    const bool igcc_is_lyp = params.igcc_is_lyp;
    const bool domag = params.domag;
    const bool domag_z = params.domag_z;
    const double hybrid_alpha_in = params.hybrid_alpha;
    const double hse_omega_in = params.hse_omega;
    const bool use_libxc = params.use_libxc;
    const int func_type = params.func_type;
    const std::vector<int>& func_id = params.func_id;

    const double epsr = 1.0e-6;

#ifdef _OPENMP
#pragma omp parallel
    {
        std::vector<double> local_stress_gga;
        double local_vtxcgc = 0.0;
        double local_etxcgc = 0.0;

        if(is_stress)
        {
            local_stress_gga.resize(9);
            for(int i=0;i<9;i++)
            {
                local_stress_gga[i] = 0.0;
            }
        }
#else
    std::vector<double> &local_stress_gga = stress_gga;
    double &local_vtxcgc = vtxcgc;
    double &local_etxcgc = etxcgc;
#endif

        double grho2a = 0.0;
        double grho2b = 0.0;
        double sxc = 0.0;
        double v1xc = 0.0;
        double v2xc = 0.0;

        if(nspin0==1)
        {
            double segno = 0.0;
#ifdef _OPENMP
#pragma omp for
#endif
            for(int ir=0; ir<rhopw->nrxx; ir++)
            {
                const double arho = std::abs( buf.rhotmp1[ir] );
                if(!is_stress)
                {
                    buf.h1[ir].x = 0.0;
                    buf.h1[ir].y = 0.0;
                    buf.h1[ir].z = 0.0;
                }

                if(arho > epsr)
                {
                    grho2a = buf.gdr1[ir].norm2();

                    //normally values in rhotmp can either be >= 0 or < 0.
                    if( buf.rhotmp1[ir] >= 0.0 )
                    {
                        segno = 1.0;
                    }
                    else
                    {
                        segno = -1.0;
                    }
                    if (use_libxc && is_stress)
                    {
#ifdef __LIBXC
                        if(func_type == 3 || func_type == 5)
                        {
                            double v3xc = 0.0;
                            double vlaplxc = 0.0;
                            double atau = chr->kin_r[0][ir]/2.0;
                            double lapl_val = (!buf.lapl1.empty()) ? buf.lapl1[ir] : 0.0;
                            XC_Functional_Libxc::tau_xc( func_id, arho, grho2a, lapl_val, atau, sxc, v1xc, v2xc, v3xc, vlaplxc, hybrid_alpha_in, hse_omega_in);
                            if(!buf.vlapl_arr1.empty()) buf.vlapl_arr1[ir] = vlaplxc;
                        }
                        else
                        {
                            XC_Functional_Libxc::gcxc_libxc( func_id, arho, grho2a, sxc, v1xc, v2xc, hybrid_alpha_in, hse_omega_in);
                        }
#endif
                    }
                    else
                    {
                        XC_Functional::gcxc( arho, grho2a, sxc, v1xc, v2xc);
                    }
                    if(is_stress)
                    {
                        double tt[3];
                        tt[0] = buf.gdr1[ir].x;
                        tt[1] = buf.gdr1[ir].y;
                        tt[2] = buf.gdr1[ir].z;
                        for(int l = 0;l< 3;l++)
                        {
                            for(int m = 0;m< l+1;m++)
                            {
                                int ind = l*3 + m;
                                local_stress_gga[ind] += tt[l] * tt[m] * ModuleBase::e2 * v2xc;
                            }
                        }
                    }
                    else
                    {
                        // first term of the gradient correction:
                        // D(rho*Exc)/D(rho)
                        v(0, ir) += ModuleBase::e2 * v1xc;

                        // h contains
                        // D(rho*Exc) / D(|grad rho|) * (grad rho) / |grad rho|
                        buf.h1[ir] = ModuleBase::e2 * v2xc * buf.gdr1[ir];

                        local_vtxcgc += ModuleBase::e2* v1xc * ( buf.rhotmp1[ir] - chr->rho_core[ir] );
                        local_etxcgc += ModuleBase::e2* sxc  * segno;
                    }
                }
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp for
#endif
            for(int ir=0; ir<rhopw->nrxx; ir++)
            {
                if(use_libxc)
                {
#ifdef __LIBXC
                    double sxc = 0.0;
                    double v1xcup = 0.0;
                    double v1xcdw = 0.0;
                    double v2xcup = 0.0;
                    double v2xcdw = 0.0;
                    double v2xcud = 0.0;
                    if(func_type == 3 || func_type == 5)
                    {
                        double v3xcup = 0.0;
                        double v3xcdw = 0.0;
                        double vlaplxcup = 0.0;
                        double vlaplxcdw = 0.0;
                        double atau1 = chr->kin_r[0][ir]/2.0;
                        double atau2 = chr->kin_r[1][ir]/2.0;
                        double laplup_val = (!buf.lapl1.empty()) ? buf.lapl1[ir] : 0.0;
                        double lapldw_val = (!buf.lapl2.empty()) ? buf.lapl2[ir] : 0.0;
                        XC_Functional_Libxc::tau_xc_spin(
                            func_id,
                            buf.rhotmp1[ir], buf.rhotmp2[ir], buf.gdr1[ir], buf.gdr2[ir],
                            laplup_val, lapldw_val, atau1, atau2, sxc, v1xcup, v1xcdw, v2xcup, v2xcdw, v2xcud, v3xcup, v3xcdw, vlaplxcup, vlaplxcdw, hybrid_alpha_in, hse_omega_in);
                        if(!buf.vlapl_arr1.empty()) buf.vlapl_arr1[ir] = vlaplxcup;
                        if(!buf.vlapl_arr2.empty()) buf.vlapl_arr2[ir] = vlaplxcdw;
                    }
                    else
                    {
                        XC_Functional_Libxc::gcxc_spin_libxc(
                            func_id,
                            buf.rhotmp1[ir], buf.rhotmp2[ir], buf.gdr1[ir], buf.gdr2[ir],
                            sxc, v1xcup, v1xcdw, v2xcup, v2xcdw, v2xcud,
                            hybrid_alpha_in, hse_omega_in);
                    }
                    if(is_stress)
                    {
                        double tt1[3],tt2[3];
                        {
                            tt1[0] = buf.gdr1[ir].x;
                            tt1[1] = buf.gdr1[ir].y;
                            tt1[2] = buf.gdr1[ir].z;
                            tt2[0] = buf.gdr2[ir].x;
                            tt2[1] = buf.gdr2[ir].y;
                            tt2[2] = buf.gdr2[ir].z;
                        }
                        for(int l = 0;l< 3;l++)
                        {
                            for(int m = 0;m< l+1;m++)
                            {
                                int ind = l*3 + m;
                                local_stress_gga [ind] += ( tt1[l] * tt1[m] * v2xcup +
                                    tt2[l] * tt2[m] * v2xcdw +
                                    (tt1[l] * tt2[m] + tt2[l] * tt1[m] ) * v2xcud ) * ModuleBase::e2;
                            }
                        }
                    }
                    else
                    {
                        // first term of the gradient correction : D(rho*Exc)/D(rho)
                        v(0,ir) += ModuleBase::e2 * v1xcup;
                        v(1,ir) += ModuleBase::e2 * v1xcdw;

                        // h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
                        buf.h1[ir] += ModuleBase::e2 * ( v2xcup * buf.gdr1[ir] + v2xcud * buf.gdr2[ir] );
                        buf.h2[ir] += ModuleBase::e2 * ( v2xcdw * buf.gdr2[ir] + v2xcud * buf.gdr1[ir] );

                        local_vtxcgc = local_vtxcgc + ModuleBase::e2 * v1xcup * ( buf.rhotmp1[ir] - chr->rho_core[ir] * fac );
                        local_vtxcgc = local_vtxcgc + ModuleBase::e2 * v1xcdw * ( buf.rhotmp2[ir] - chr->rho_core[ir] * fac );
                        local_etxcgc = local_etxcgc + ModuleBase::e2 * sxc;
                    }
#endif
                }
                else
                {
                    double v1cup = 0.0;
                    double v1cdw = 0.0;
                    double v2cup = 0.0;
                    double v2cdw = 0.0;
                    double v1xup = 0.0;
                    double v1xdw = 0.0;
                    double v2xup = 0.0;
                    double v2xdw = 0.0;
                    double v2cud = 0.0;
                    double v2c = 0.0;
                    double sx = 0.0;
                    double sc = 0.0;
                    double rh = buf.rhotmp1[ir] + buf.rhotmp2[ir];
                    grho2a = buf.gdr1[ir].norm2();
                    grho2b = buf.gdr2[ir].norm2();
                    XC_Functional::gcx_spin(buf.rhotmp1[ir], buf.rhotmp2[ir], grho2a, grho2b,
                        sx, v1xup, v1xdw, v2xup, v2xdw);

                    if(rh > epsr)
                    {
                        if(igcc_is_lyp)
                        {
                            ModuleBase::WARNING_QUIT("XC_Functional","igcc_is_lyp is not available now.");
                        }
                        else
                        {
                            double zeta = ( buf.rhotmp1[ir] - buf.rhotmp2[ir] ) / rh;
                            if(nspin==4&&(domag||domag_z))
                            {
                                zeta = fabs(zeta) * buf.neg[ir];
                            }
                            const double grh2 = (buf.gdr1[ir]+buf.gdr2[ir]).norm2();
                            XC_Functional::gcc_spin(rh, zeta, grh2, sc, v1cup, v1cdw, v2c);
                            v2cup = v2c;
                            v2cdw = v2c;
                            v2cud = v2c;
                        }
                    }
                    else
                    {
                        sc = 0.0;
                        v1cup = 0.0;
                        v1cdw = 0.0;
                        v2c = 0.0;
                        v2cup = 0.0;
                        v2cdw = 0.0;
                        v2cud = 0.0;
                    }

                    if(is_stress)
                    {
                        double tt1[3],tt2[3];
                        {
                            tt1[0] = buf.gdr1[ir].x;
                            tt1[1] = buf.gdr1[ir].y;
                            tt1[2] = buf.gdr1[ir].z;
                            tt2[0] = buf.gdr2[ir].x;
                            tt2[1] = buf.gdr2[ir].y;
                            tt2[2] = buf.gdr2[ir].z;
                        }
                        for(int l = 0;l< 3;l++)
                        {
                            for(int m = 0;m< l+1;m++)
                            {
                                int ind = l*3 + m;
                                // exchange
                                local_stress_gga [ind] += tt1[l] * tt1[m] * ModuleBase::e2 * v2xup +
                                    tt2[l] * tt2[m] * ModuleBase::e2 * v2xdw;
                                // correlation
                                local_stress_gga [ind] += ( tt1[l] * tt1[m] * v2cup +
                                    tt2[l] * tt2[m] * v2cdw +
                                    (tt1[l] * tt2[m] + tt2[l] * tt1[m] ) * v2cud ) * ModuleBase::e2;
                            }
                        }
                    }
                    else
                    {
                        // first term of the gradient correction : D(rho*Exc)/D(rho)
                        v(0,ir) = v(0,ir) + ModuleBase::e2 * ( v1xup + v1cup );
                        v(1,ir) = v(1,ir) + ModuleBase::e2 * ( v1xdw + v1cdw );

                        // h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
                        buf.h1[ir] = ModuleBase::e2 * ( ( v2xup + v2cup ) * buf.gdr1[ir] + v2cud * buf.gdr2[ir] );
                        buf.h2[ir] = ModuleBase::e2 * ( ( v2xdw + v2cdw ) * buf.gdr2[ir] + v2cud * buf.gdr1[ir] );

                        local_vtxcgc = local_vtxcgc + ModuleBase::e2 * ( v1xup + v1cup ) * ( buf.rhotmp1[ir] - chr->rho_core[ir] * fac );
                        local_vtxcgc = local_vtxcgc + ModuleBase::e2 * ( v1xdw + v1cdw ) * ( buf.rhotmp2[ir] - chr->rho_core[ir] * fac );
                        local_etxcgc = local_etxcgc + ModuleBase::e2 * ( sx + sc );
                    }
                }
            }
        }
#ifdef _OPENMP
    #pragma omp critical(xc_grad_reduce)
    {
        if(is_stress)
        {
            for(int l = 0;l< 3;l++)
            {
                for(int m = 0;m< l+1;m++)
                {
                    int ind = l*3 + m;
                    stress_gga [ind] += local_stress_gga [ind];
                }
            }
        }
        else
        {
            vtxcgc += local_vtxcgc;
            etxcgc += local_etxcgc;
        }
    }
}
#endif
}
