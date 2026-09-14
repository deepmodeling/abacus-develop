// First stage of the gradient-correction pipeline: build the real- and
// reciprocal-space density buffers (rhotmp/rhogsum) and their gradients
// (gdr) for every spin channel, including the non-collinear spin case.
// The stage functions and their shared parameter/buffer structs are
// declared in xc_grad_internal.h; the public entry point
// XC_Functional::gradcorr lives in xc_grad.cpp.

#include "xc_functional.h"
#include "xc_grad_internal.h"

void gradcorr_prepare_rho(
    const GradCorrParams& params,
    ModuleBase::matrix& v,
    GradCorrBuffers& buf)
{
    const Charge* const chr = params.chr;
    ModulePW::PW_Basis* rhopw = params.rhopw;
    const UnitCell* ucell = params.ucell;
    const int nspin = params.nspin;
    const int nspin0 = params.nspin0;
    const double fac = params.fac;
    const bool need_laplacian = params.need_laplacian;
    const bool is_stress = params.is_stress;
    const bool domag = params.domag;
    const bool domag_z = params.domag_z;

    // doing FFT to get rho in G space: rhog1
    rhopw->real2recip(chr->rho[0], chr->rhog[0]);
    if(nspin==2)
    {
        rhopw->real2recip(chr->rho[1], chr->rhog[1]);
    }
    rhopw->real2recip(chr->rho_core, chr->rhog_core);

    // sum up (rho_core+rho) for each spin in real space
    // and reciprocal space.

    // for spin unpolarized case,
    // calculate the gradient of (rho_core+rho) in reciprocal space.
    buf.rhotmp1.resize(rhopw->nrxx);
    buf.rhogsum1.resize(rhopw->npw);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int ir=0; ir<rhopw->nrxx; ir++)
    {
        buf.rhotmp1[ir] = chr->rho[0][ir] + fac * chr->rho_core[ir];
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int ig=0; ig<rhopw->npw; ig++)
    {
        buf.rhogsum1[ig] = chr->rhog[0][ig] + fac * chr->rhog_core[ig];
    }

    buf.gdr1.resize(rhopw->nrxx);
    if(!is_stress)
    {
        buf.h1.resize(rhopw->nrxx);
    }

    XC_Functional::grad_rho(buf.rhogsum1.data(), buf.gdr1.data(), rhopw, ucell->tpiba);

    if(need_laplacian)
    {
        buf.lapl1.resize(rhopw->nrxx);
        XC_Functional::laplacian_rho(buf.rhogsum1.data(), buf.lapl1.data(), rhopw, ucell->tpiba);
        if(params.use_libxc)
        {
            buf.vlapl_arr1.resize(rhopw->nrxx, 0.0);
        }
    }

    // for spin polarized case;
    // calculate the gradient of (rho_core+rho) in reciprocal space.
    if(nspin==2)
    {
        buf.rhotmp2.resize(rhopw->nrxx);
        buf.rhogsum2.resize(rhopw->npw);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir=0; ir<rhopw->nrxx; ir++)
        {
            buf.rhotmp2[ir] = chr->rho[1][ir] + fac * chr->rho_core[ir];
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig=0; ig<rhopw->npw; ig++)
        {
            buf.rhogsum2[ig] = chr->rhog[1][ig] + fac * chr->rhog_core[ig];
        }

        buf.gdr2.resize(rhopw->nrxx);
        if(!is_stress)
        {
            buf.h2.resize(rhopw->nrxx);
        }

        XC_Functional::grad_rho(buf.rhogsum2.data(), buf.gdr2.data(), rhopw, ucell->tpiba);

        if(need_laplacian)
        {
            buf.lapl2.resize(rhopw->nrxx);
            XC_Functional::laplacian_rho(buf.rhogsum2.data(), buf.lapl2.data(), rhopw, ucell->tpiba);
            if(params.use_libxc)
            {
                buf.vlapl_arr2.resize(rhopw->nrxx, 0.0);
            }
        }
    }

    if(nspin == 4&&(domag||domag_z))
    {
        buf.rhotmp2.resize(rhopw->nrxx);
        buf.rhogsum2.resize(rhopw->npw);
        buf.neg.resize(rhopw->nrxx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir=0; ir<rhopw->nrxx; ir++)
        {
            buf.rhotmp1[ir] = 0.0;
            buf.rhotmp2[ir] = 0.0;
            buf.neg[ir] = 0.0;
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig=0; ig<rhopw->npw; ig++)
        {
            buf.rhogsum1[ig] = 0.0;
            buf.rhogsum2[ig] = 0.0;
        }
        if(!is_stress)
        {
            buf.vsave.assign(nspin * rhopw->nrxx, 0.0);
            buf.vgg.assign(nspin0 * rhopw->nrxx, 0.0);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
            for(int is = 0;is<nspin;is++)
            {
                for(int ir =0;ir<rhopw->nrxx;ir++)
                {
                    buf.vsave[is * rhopw->nrxx + ir] = v(is,ir);
                    v(is,ir) = 0;
                }
            }
        }
        XC_Functional::noncolin_rho(buf.rhotmp1.data(), buf.rhotmp2.data(), buf.neg.data(), chr->rho, rhopw->nrxx, ucell->magnet.ux_, ucell->magnet.lsign_);
        rhopw->real2recip(buf.rhotmp1.data(), buf.rhogsum1.data());
        rhopw->real2recip(buf.rhotmp2.data(), buf.rhogsum2.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir=0; ir<rhopw->nrxx; ir++)
        {
            buf.rhotmp2[ir] += fac * chr->rho_core[ir];
            buf.rhotmp1[ir] += fac * chr->rho_core[ir];
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig=0; ig<rhopw->npw; ig++)
        {
            buf.rhogsum2[ig] += fac * chr->rhog_core[ig];
            buf.rhogsum1[ig] += fac * chr->rhog_core[ig];
        }

        buf.gdr2.resize(rhopw->nrxx);
        if(!is_stress)
        {
            buf.h2.resize(rhopw->nrxx);
        }

        XC_Functional::grad_rho(buf.rhogsum1.data(), buf.gdr1.data(), rhopw, ucell->tpiba);
        XC_Functional::grad_rho(buf.rhogsum2.data(), buf.gdr2.data(), rhopw, ucell->tpiba);

        if(need_laplacian)
        {
            buf.lapl1.resize(rhopw->nrxx);
            XC_Functional::laplacian_rho(buf.rhogsum1.data(), buf.lapl1.data(), rhopw, ucell->tpiba);
            buf.lapl2.resize(rhopw->nrxx);
            XC_Functional::laplacian_rho(buf.rhogsum2.data(), buf.lapl2.data(), rhopw, ucell->tpiba);
            if(params.use_libxc)
            {
                buf.vlapl_arr1.resize(rhopw->nrxx, 0.0);
                buf.vlapl_arr2.resize(rhopw->nrxx, 0.0);
            }
        }
    }
}
