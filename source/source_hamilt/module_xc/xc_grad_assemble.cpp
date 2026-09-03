// Final stage of the gradient-correction pipeline: assemble the kernel
// contribution into vtxc/etxc/stress_gga/v, subtract the core density,
// and restore the nspin = 4 non-collinear case via spin rotation.
// The stage functions and their shared parameter/buffer structs are
// declared in xc_grad_internal.h; the public entry point
// XC_Functional::gradcorr lives in xc_grad.cpp.

#include "xc_functional.h"
#include "xc_grad_internal.h"
#include "source_base/constants.h"

void gradcorr_assemble_vxc(const GradCorrParams& params,
                           double vtxcgc,
                           double etxcgc,
                           double& vtxc,
                           double& etxc,
                           std::vector<double>& stress_gga,
                           ModuleBase::matrix& v,
                           GradCorrBuffers& buf)
{
    const Charge* const chr = params.chr;
    ModulePW::PW_Basis* rhopw = params.rhopw;
    const UnitCell* ucell = params.ucell;
    const int nspin = params.nspin;
    const int nspin0 = params.nspin0;
    const double fac = params.fac;
    const bool is_stress = params.is_stress;
    const bool domag = params.domag;
    const bool domag_z = params.domag_z;
    const bool use_libxc = params.use_libxc;

    // Add Laplacian stress contribution from meta-GGA functionals
    if(is_stress && use_libxc && !buf.vlapl_arr1.empty())
    {
        const int ng = rhopw->npw;
        const int nrxx = rhopw->nrxx;
        const double tpiba2 = ucell->tpiba * ucell->tpiba;
        std::vector<std::complex<double>> vlapl_g(rhopw->nmaxgr);

        for(int is = 0; is < nspin0; is++)
        {
            double* vlapl_ptr = (is == 0) ? buf.vlapl_arr1.data() : buf.vlapl_arr2.data();
            double* rho_ptr = (is == 0) ? buf.rhotmp1.data() : buf.rhotmp2.data();
            if(vlapl_ptr == nullptr || rho_ptr == nullptr) continue;

            std::vector<std::complex<double>> rho_g(rhopw->nmaxgr);
            for(int ir = 0; ir < nrxx; ir++)
                rho_g[ir] = std::complex<double>(rho_ptr[ir], 0.0);
            rhopw->real2recip(rho_g.data(), rho_g.data());

            for(int ir = 0; ir < nrxx; ir++)
                vlapl_g[ir] = std::complex<double>(vlapl_ptr[ir], 0.0);
            rhopw->real2recip(vlapl_g.data(), vlapl_g.data());

            for(int l = 0; l < 3; l++)
            {
                for(int m = 0; m <= l; m++)
                {
                    double sum = 0.0;
                    for(int ig = 0; ig < ng; ig++)
                    {
                        double g_prod = rhopw->gcar[ig][l] * rhopw->gcar[ig][m] * tpiba2;
                        sum += g_prod * (rho_g[ig].real() * vlapl_g[ig].real()
                                       + rho_g[ig].imag() * vlapl_g[ig].imag());
                    }
                    stress_gga[l*3+m] -= 2.0 * static_cast<double>(rhopw->nxyz) * ModuleBase::e2 * sum;
                }
            }
        }
    }

    if(!is_stress)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir=0; ir<rhopw->nrxx; ir++)
        {
            buf.rhotmp1[ir] -= fac * chr->rho_core[ir];
        }
        if(nspin0==2)
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for(int ir=0; ir<rhopw->nrxx; ir++)
            {
                buf.rhotmp2[ir] -= fac * chr->rho_core[ir];
            }
        }

        // second term of the gradient correction :
        // sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )

        // dh is in real sapce.
        std::vector<double> dh(rhopw->nrxx);

        for(int is=0; is<nspin0; is++)
        {
            if(is==0)
            {
                XC_Functional::grad_dot(buf.h1.data(), dh.data(), rhopw, ucell->tpiba);
            }
            if(is==1)
            {
                XC_Functional::grad_dot(buf.h2.data(), dh.data(), rhopw, ucell->tpiba);
            }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for(int ir=0; ir<rhopw->nrxx; ir++)
            {
                v(is, ir) -= dh[ir];
            }

            double sum = 0.0;
            if(is==0)
            {
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum) schedule(static, 256)
#endif
                for(int ir=0; ir<rhopw->nrxx; ir++)
                {
                    sum += dh[ir] * buf.rhotmp1[ir];
                }
            }
            else if(is==1)
            {
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum) schedule(static, 256)
#endif
                for(int ir=0; ir<rhopw->nrxx; ir++)
                {
                    sum += dh[ir] * buf.rhotmp2[ir];
                }
            }
            vtxcgc -= sum;
        }

        vtxc += vtxcgc;
        etxc += etxcgc;

        // Add Laplacian contribution from meta-GGA functionals
        // v_xc += nabla^2(vlapl) where vlapl = d(rho*eps_xc)/d(nabla^2 rho)
        if(use_libxc && !buf.vlapl_arr1.empty())
        {
            const int ng = rhopw->npw;
            const int nrxx = rhopw->nrxx;
            const double tpiba2 = ucell->tpiba * ucell->tpiba;
            std::vector<std::complex<double>> vlapl_g(rhopw->nmaxgr);

            for(int is = 0; is < nspin0; is++)
            {
                double* vlapl_ptr = (is == 0) ? buf.vlapl_arr1.data() : buf.vlapl_arr2.data();
                if(vlapl_ptr == nullptr) continue;

                for(int ir = 0; ir < nrxx; ir++)
                    vlapl_g[ir] = std::complex<double>(vlapl_ptr[ir], 0.0);
                rhopw->real2recip(vlapl_g.data(), vlapl_g.data());
                for(int ig = 0; ig < ng; ig++)
                {
                    double g2 = 0.0;
                    for(int i = 0; i < 3; i++)
                        g2 += rhopw->gcar[ig][i] * rhopw->gcar[ig][i];
                    vlapl_g[ig] *= -g2 * tpiba2;
                }
                rhopw->recip2real(vlapl_g.data(), vlapl_g.data());
                for(int ir = 0; ir < nrxx; ir++)
                {
                    v(is, ir) += ModuleBase::e2 * vlapl_g[ir].real();
                }
            }
        }

        if(nspin == 4 && (domag||domag_z))
        {
            const int nrxx = rhopw->nrxx;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
            for(int is=0;is<nspin;is++)
            {
                for(int ir=0;ir<nrxx;ir++)
                {
                    if(is<nspin0)
                    {
                        buf.vgg[is * nrxx + ir] = v(is,ir);
                    }
                    v(is,ir) = buf.vsave[is * nrxx + ir];
                }
            }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for(int ir=0;ir<nrxx;ir++)
            {
                v(0,ir) += 0.5 * (buf.vgg[ir] + buf.vgg[nrxx + ir]);
                double amag = sqrt(pow(chr->rho[1][ir],2)+pow(chr->rho[2][ir],2)+pow(chr->rho[3][ir],2));
                if(amag>1e-12)
                {
                    for(int i=1;i<4;i++)
                    {
                        v(i,ir)+= buf.neg[ir] * 0.5 *(buf.vgg[ir]-buf.vgg[nrxx + ir])*chr->rho[i][ir]/amag;
                    }
                }
            }
        }
    }
}
