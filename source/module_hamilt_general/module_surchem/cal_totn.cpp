#include "surchem.h"

void surchem::cal_totn(const UnitCell& cell,
                       const ModulePW::PW_Basis* rho_basis,
                       const complex<double>* Porter_g,
                       complex<double>* N,
                       complex<double>* TOTN,
                       const double* vlocal)
{
    // vloc to N
    complex<double> *vloc_g = new complex<double>[rho_basis->npw];
    ModuleBase::GlobalFunc::ZEROS(vloc_g, rho_basis->npw);

    rho_basis->real2recip(vlocal, vloc_g);  // now n is vloc in Recispace
    for (int ig = 0; ig < rho_basis->npw; ig++) {
        if(ig==rho_basis->ig_gge0)
        {
            N[ig] = Porter_g[ig];
            continue;
        }
        const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI /
                            (cell.tpiba2 * rho_basis->gg[ig]);

        N[ig] = -vloc_g[ig] / fac;
    }

    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        TOTN[ig] = N[ig] - Porter_g[ig];
    }

    delete[] vloc_g;
    return;
}

void surchem::induced_charge(const UnitCell& cell, const ModulePW::PW_Basis* rho_basis, double* induced_rho)
{
    std::complex<double> *delta_phig = new complex<double>[rho_basis->npw];
    std::complex<double> *induced_rhog = new complex<double>[rho_basis->npw];
    ModuleBase::GlobalFunc::ZEROS(induced_rhog, rho_basis->npw);
    rho_basis->real2recip(delta_phi, delta_phig);
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {   
        if(rho_basis->ig_gge0 == ig)
        {
            continue;
        }
        else
        {
            const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI /(cell.tpiba2 * rho_basis->gg[ig]);
            induced_rhog[ig] = -delta_phig[ig] / fac;
        }
    }

    rho_basis->recip2real(induced_rhog, induced_rho);

    delete[] delta_phig;
    delete[] induced_rhog;
    return;
}