#include "surchem.h"

namespace GlobalC
{
surchem solvent_model;
}

surchem::surchem()
{
    TOTN_real = nullptr;
    delta_phi = nullptr;
    epspot = nullptr;
    comp_real = nullptr;
    phi_comp_R = nullptr;
    Vcav = ModuleBase::matrix();
    Vel = ModuleBase::matrix();
    qs = 0;
}

void surchem::allocate(const int &nrxx, const int &nspin)
{
    assert(nrxx >= 0);
    assert(nspin > 0);

    delete[] TOTN_real;
    delete[] delta_phi;
    delete[] epspot;
    delete[] comp_real;
    delete[] phi_comp_R;
    if (nrxx > 0)
    {
        TOTN_real = new double[nrxx];
        delta_phi = new double[nrxx];
        epspot = new double[nrxx];
        comp_real = new double[nrxx];
        phi_comp_R = new double[nrxx];
    }
    else
    {
        TOTN_real = nullptr;
        delta_phi = nullptr;
        epspot = nullptr;
        comp_real = nullptr;
        phi_comp_R = nullptr;
    }
    Vcav.create(nspin, nrxx);
    Vel.create(nspin, nrxx);

    ModuleBase::GlobalFunc::ZEROS(delta_phi, nrxx);
    ModuleBase::GlobalFunc::ZEROS(TOTN_real, nrxx);
    ModuleBase::GlobalFunc::ZEROS(epspot, nrxx);
    ModuleBase::GlobalFunc::ZEROS(comp_real, nrxx);
    ModuleBase::GlobalFunc::ZEROS(phi_comp_R, nrxx);
    return;
}

surchem::~surchem()
{
    delete[] TOTN_real;
    delete[] delta_phi;
    delete[] epspot;
    delete[] comp_real;
    delete[] phi_comp_R;
}

void surchem::get_totn_reci(const UnitCell &cell, ModulePW::PW_Basis *rho_basis, complex<double> *totn_reci)
{
    double *tmp_totn_real = new double[rho_basis->nrxx];
    double *tmp_comp_real = new double[rho_basis->nrxx];
    complex<double> *comp_reci = new complex<double>[rho_basis->npw];
    ModuleBase::GlobalFunc::ZEROS(tmp_totn_real, rho_basis->nrxx);
    ModuleBase::GlobalFunc::ZEROS(tmp_comp_real, rho_basis->nrxx);
    ModuleBase::GlobalFunc::ZEROS(comp_reci, rho_basis->npw);
    add_comp_chg(cell, rho_basis, comp_q, comp_l, comp_center, comp_reci, comp_dim, false);
    rho_basis->recip2real(comp_reci, tmp_comp_real);

    for (int ir = 0; ir < rho_basis->nrxx;ir++)
    {
        tmp_totn_real[ir] = TOTN_real[ir] + tmp_comp_real[ir];
    }

    rho_basis->real2recip(tmp_totn_real, totn_reci);
    delete[] tmp_totn_real;
    delete[] tmp_comp_real;
    delete[] comp_reci;
}