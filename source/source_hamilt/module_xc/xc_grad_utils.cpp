// Auxiliary real/reciprocal-space utilities shared by the gradient
// correction and meta-GGA stress:
// 1. grad_rho      : gradient of the density
// 2. grad_dot      : divergence of a vector field
// 3. laplacian_rho : Laplacian of the density
// 4. noncolin_rho  : diagonalizes the spin density matrix and gives the
//                    spin up and spin down components of the charge.

#include "xc_functional.h"
#include "source_base/constants.h"

void XC_Functional::grad_rho(
    const std::complex<double>* rhog,
    ModuleBase::Vector3<double>* gdr,
    const ModulePW::PW_Basis* rho_basis,
    const double tpiba)
{
    std::vector<std::complex<double>> gdrtmp(rho_basis->nmaxgr);

    // the formula is : rho(r)^prime = int iG * rho(G)e^{iGr} dG
    for(int i = 0 ; i < 3 ; ++i)
    {
        // calculate the charge density gradient in reciprocal space.
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig=0; ig<rho_basis->npw; ig++)
        {
            gdrtmp[ig] = ModuleBase::IMAG_UNIT * rhog[ig] * rho_basis->gcar[ig][i];
        }

        // bring the gdr from G --> R
        rho_basis->recip2real(gdrtmp.data(), gdrtmp.data());

        // remember to multily 2pi/a0, which belongs to G vectors.
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir=0; ir<rho_basis->nrxx; ir++)
        {
            gdr[ir][i] = gdrtmp[ir].real() * tpiba;
        }
    }

    return;
}


void XC_Functional::grad_dot(
    const ModuleBase::Vector3<double>* h,
    double* dh,
    const ModulePW::PW_Basis* rho_basis,
    const double tpiba)
{
    std::vector<std::complex<double>> aux(rho_basis->nmaxgr);
    std::vector<std::complex<double>> gaux(rho_basis->npw);

    for(int i = 0 ; i < 3 ; ++i)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir = 0; ir < rho_basis->nrxx; ++ir)
        {
            aux[ir] = std::complex<double>( h[ir][i], 0.0);
        }

        // bring to G space.
        rho_basis->real2recip(aux.data(), aux.data());
        if (i == 0)
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for(int ig = 0; ig < rho_basis->npw; ++ig)
            {
                gaux[ig] =  ModuleBase::IMAG_UNIT * aux[ig] * rho_basis->gcar[ig][i];
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for(int ig = 0; ig < rho_basis->npw; ++ig)
            {
                gaux[ig] +=  ModuleBase::IMAG_UNIT * aux[ig] * rho_basis->gcar[ig][i];
            }
        }
    }

    // bring back to R space
    rho_basis->recip2real(gaux.data(), aux.data());

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int ir=0; ir<rho_basis->nrxx; ir++)
    {
        dh[ir] = aux[ir].real() * tpiba;
    }

    return;
}

void XC_Functional::laplacian_rho(
    const std::complex<double>* rhog,
    double* lapl,
    const ModulePW::PW_Basis* rho_basis,
    const double tpiba)
{
    std::vector<std::complex<double>> lapl_tmp(rho_basis->nmaxgr);

    for(int ig=0; ig<rho_basis->npw; ig++)
    {
        double g2 = 0.0;
        for(int i=0; i<3; i++)
        {
            g2 += rho_basis->gcar[ig][i] * rho_basis->gcar[ig][i];
        }
        lapl_tmp[ig] = -rhog[ig] * g2;
    }
    rho_basis->recip2real(lapl_tmp.data(), lapl_tmp.data());
    for(int ir=0; ir<rho_basis->nrxx; ir++)
    {
        lapl[ir] = lapl_tmp[ir].real() * tpiba * tpiba;
    }
}

void XC_Functional::noncolin_rho(
    double *rhoout1,
    double *rhoout2,
    double *neg,
    const double*const*const rho,
    const int nrxx,
    const double* ux_,
    const bool lsign_)
{
    //this function diagonalizes the spin density matrix and gives as output the
    //spin up and spin down components of the charge.
    //If lsign is true up and dw are with respect to the fixed quantization axis
    //ux, otherwise rho + |m| is always rhoup and rho-|m| is always rhodw.
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int ir = 0;ir<nrxx;ir++)
    {
        neg[ir] = 1.0;
    }
    if(lsign_)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir = 0;ir<nrxx;ir++)
        {
            if(rho[1][ir]*ux_[0] + rho[2][ir]*ux_[1] + rho[3][ir]*ux_[2]>0)
            {
                neg[ir] = 1.0;
            }
            else
            {
                neg[ir] = -1.0;
            }
        }
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int ir = 0;ir<nrxx;ir++)
    {
        double amag = sqrt(pow(rho[1][ir],2)+pow(rho[2][ir],2)+pow(rho[3][ir],2));
        rhoout1[ir] = 0.5 * (rho[0][ir] + neg[ir] * amag);
        rhoout2[ir] = 0.5 * (rho[0][ir] - neg[ir] * amag);
    }
    return;
}
