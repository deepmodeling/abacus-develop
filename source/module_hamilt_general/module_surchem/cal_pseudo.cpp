#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "surchem.h"

// atom_in surchem::GetAtom;

void surchem::gauss_charge(const UnitCell& cell,
                           const ModulePW::PW_Basis* rho_basis,
                           complex<double>* N,
                           Structure_Factor* sf)
{
    sf->setup_structure_factor(&GlobalC::ucell, rho_basis); // call strucFac(ntype,ngmc)
    for (int it = 0; it < cell.ntype; it++)
    {
        double RCS = GetAtom.atom_RCS[cell.atoms[it].ncpp.psd];
        // double RCS = GetAtom.get_RCS(cell.atoms[it].psd);
        double sigma_rc_k = RCS / 2.5;
        double sigma_nc_k = RCS / 10.0;
        for (int ig = 0; ig < rho_basis->npw; ig++)
        {
            // G^2
            double gg = rho_basis->gg[ig];
            gg = gg * cell.tpiba2;

            N[ig].real(N[ig].real()
                       + (GetAtom.atom_Z[cell.atoms[it].ncpp.psd] - cell.atoms[it].ncpp.zv)
                             * sf->strucFac(it, ig).real() * exp(-0.5 * gg * (sigma_rc_k * sigma_rc_k)));
            N[ig].imag(N[ig].imag()
                       + (GetAtom.atom_Z[cell.atoms[it].ncpp.psd] - cell.atoms[it].ncpp.zv)
                             * sf->strucFac(it, ig).imag() * exp(-0.5 * gg * (sigma_rc_k * sigma_rc_k)));
        }
    }
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        N[ig] /= cell.omega;
    }
}

void surchem::cal_pseudo(const UnitCell& cell,
                         const ModulePW::PW_Basis* rho_basis,
                         const complex<double>* Porter_g,
                         complex<double>* PS_TOTN,
                         Structure_Factor* sf)
{
    complex<double>* N = new complex<double>[rho_basis->npw];
    ModuleBase::GlobalFunc::ZEROS(N, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(PS_TOTN, rho_basis->npw);

    gauss_charge(cell, rho_basis, N, sf);

    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        PS_TOTN[ig] = N[ig] + Porter_g[ig];
    }

    delete[] N;
}
