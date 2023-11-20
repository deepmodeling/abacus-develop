#include "./kedf_vw.h"

#include <iostream>

#include "module_base/parallel_reduce.h"

void KEDF_vW::set_para(double dV, double vw_weight)
{
    this->dV_ = dV;
    this->vw_weight_ = vw_weight;
}

//
// EvW = -1/2 \int{sqrt(rho) nabla^2 sqrt(rho)}
//
double KEDF_vW::get_energy(double **pphi, ModulePW::PW_Basis *pw_rho)
{
    // since pphi may contain minus element, we define tempPhi = std::abs(phi), which is true sqrt(rho)
    double **tempPhi = new double *[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        tempPhi[is] = new double[pw_rho->nrxx];
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            tempPhi[is][ir] = std::abs(pphi[is][ir]);
        }
    }

    double **LapPhi = new double *[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
        LapPhi[is] = new double[pw_rho->nrxx];
    this->laplacian_phi(tempPhi, LapPhi, pw_rho);

    double energy = 0.; // in Ry
    if (GlobalV::NSPIN == 1)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            energy += tempPhi[0][ir] * LapPhi[0][ir];
        }
        energy *= this->dV_ * 0.5 * this->vw_weight_ * 2.; // vw_weight * 2 to convert Hartree to Ry
    }
    else if (GlobalV::NSPIN == 2)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ir = 0; ir < pw_rho->nrxx; ++ir)
            {
                energy += 2 * tempPhi[is][ir] * LapPhi[is][ir];
            }
        }
        energy *= 0.5 * this->dV_ * 0.5 * this->vw_weight_ * 2.; // vw_weight * 2 to convert Hartree to Ry
    }
    this->vw_energy = energy;
    Parallel_Reduce::reduce_all(this->vw_energy);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] tempPhi[is];
        delete[] LapPhi[is];
    }
    delete[] tempPhi;
    delete[] LapPhi;

    return energy;
}

double KEDF_vW::get_energy_density(double **pphi, int is, int ir, ModulePW::PW_Basis *pw_rho)
{
    // since pphi may contain minus element, we define tempPhi = std::abs(phi), which is true sqrt(rho)
    double **tempPhi = new double *[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        tempPhi[is] = new double[pw_rho->nrxx];
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            tempPhi[is][ir] = std::abs(pphi[is][ir]);
        }
    }

    double **LapPhi = new double *[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
        LapPhi[is] = new double[pw_rho->nrxx];
    this->laplacian_phi(tempPhi, LapPhi, pw_rho);

    double energyDen = 0.;                                                     // in Ry
    energyDen = 0.5 * tempPhi[is][ir] * LapPhi[is][ir] * this->vw_weight_ * 2.; // vw_weight * 2 to convert Hartree to Ry

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] tempPhi[is];
        delete[] LapPhi[is];
    }
    delete[] tempPhi;
    delete[] LapPhi;
    return energyDen;
}

//
// Vvw(r)=-1/2 * [\nabla^2 sqrt(rho(r))]/sqrt(rho(r))
// Here we calculate Vvw * 2 * phi = - sign(phi) * \nabla^2 sqrt(rho(r))
//
void KEDF_vW::vw_potential(const double *const *pphi, ModulePW::PW_Basis *pw_rho, ModuleBase::matrix &rpotential)
{
    ModuleBase::timer::tick("KEDF_vW", "vw_potential");

    // since pphi may contain minus element, we define tempPhi = std::abs(phi), which is true sqrt(rho)
    double **tempPhi = new double *[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        tempPhi[is] = new double[pw_rho->nrxx];
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            tempPhi[is][ir] = std::abs(pphi[is][ir]);
        }
    }

    // calculate the minus \nabla^2 sqrt(rho)
    double **LapPhi = new double *[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
        LapPhi[is] = new double[pw_rho->nrxx];
    this->laplacian_phi(tempPhi, LapPhi, pw_rho);

    // calculate potential
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            if (pphi[is][ir] >= 0)
            {
                rpotential(is, ir) += LapPhi[is][ir] * this->vw_weight_ * 2.; // vw_weight * 2 to convert Hartree to Ry
            }
            else
            {
                rpotential(is, ir) += -LapPhi[is][ir] * this->vw_weight_ * 2.; // vw_weight * 2 to convert Hartree to Ry
            }
        }
    }

    // calculate energy
    double energy = 0.; // in Ry
    if (GlobalV::NSPIN == 1)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            energy += tempPhi[0][ir] * LapPhi[0][ir];
        }
        energy *= this->dV_ * 0.5 * this->vw_weight_ * 2.; // vw_weight * 2 to convert Hartree to Ry
    }
    else if (GlobalV::NSPIN == 2)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ir = 0; ir < pw_rho->nrxx; ++ir)
            {
                energy += 2 * tempPhi[is][ir] * LapPhi[is][ir];
            }
        }
        energy *= 0.5 * this->dV_ * 0.5 * this->vw_weight_ * 2.; // vw_weight * 2 to convert Hartree to Ry
    }
    this->vw_energy = energy;
    Parallel_Reduce::reduce_all(this->vw_energy);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] tempPhi[is];
        delete[] LapPhi[is];
    }
    delete[] tempPhi;
    delete[] LapPhi;

    ModuleBase::timer::tick("KEDF_vW", "vw_potential");
}

void KEDF_vW::get_stress(const double *const *pphi, ModulePW::PW_Basis *pw_rho, double inpt_vWenergy)
{
    // since pphi may contain minus element, we define tempPhi = std::abs(phi), which is true sqrt(rho)
    double **tempPhi = new double *[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        tempPhi[is] = new double[pw_rho->nrxx];
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            tempPhi[is][ir] = std::abs(pphi[is][ir]);
        }
    }

    std::complex<double> **recipPhi = new std::complex<double> *[GlobalV::NSPIN];
    std::complex<double> **ggrecipPhi = new std::complex<double> *[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        recipPhi[is] = new std::complex<double>[pw_rho->npw];
        ggrecipPhi[is] = new std::complex<double>[pw_rho->npw];

        pw_rho->real2recip(tempPhi[is], recipPhi[is]);
    }

    double *ggPhi = new double[pw_rho->nrxx];

    for (int alpha = 0; alpha < 3; ++alpha)
    {
        for (int beta = alpha; beta < 3; ++beta)
        {
            this->stress(alpha, beta) = 0;
            for (int is = 0; is < GlobalV::NSPIN; ++is)
            {
                for (int ik = 0; ik < pw_rho->npw; ++ik)
                {
                    ggrecipPhi[is][ik]
                        = -recipPhi[is][ik] * pw_rho->gcar[ik][alpha] * pw_rho->gcar[ik][beta] * pw_rho->tpiba2;
                }
                pw_rho->recip2real(ggrecipPhi[is], ggPhi);
                for (int ir = 0; ir < pw_rho->nrxx; ++ir)
                {
                    this->stress(alpha, beta) += tempPhi[is][ir] * ggPhi[ir];
                }
            }
            Parallel_Reduce::reduce_all(this->stress(alpha, beta));
            this->stress(alpha, beta)
                *= -1. * this->vw_weight_ * 2. / pw_rho->nxyz; // vw_weight * 2 to convert Hartree to Ry
        }
    }
    for (int alpha = 1; alpha < 3; ++alpha)
    {
        for (int beta = 0; beta < alpha; ++beta)
        {
            this->stress(alpha, beta) = this->stress(beta, alpha);
        }
    }
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] tempPhi[is];
        delete[] recipPhi[is];
        delete[] ggrecipPhi[is];
    }
    delete[] tempPhi;
    delete[] recipPhi;
    delete[] ggrecipPhi;
    delete[] ggPhi;
}

// get minus Laplacian phi
void KEDF_vW::laplacian_phi(const double *const *pphi, double **rLapPhi, ModulePW::PW_Basis *pw_rho)
{
    std::complex<double> **recipPhi = new std::complex<double> *[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        recipPhi[is] = new std::complex<double>[pw_rho->npw];

        pw_rho->real2recip(pphi[is], recipPhi[is]);
        for (int ik = 0; ik < pw_rho->npw; ++ik)
        {
            recipPhi[is][ik] *= pw_rho->gg[ik] * pw_rho->tpiba2;
        }
        pw_rho->recip2real(recipPhi[is], rLapPhi[is]);
    }

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] recipPhi[is];
    }
    delete[] recipPhi;
}