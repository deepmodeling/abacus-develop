#include "surchem.h"

double surchem::cal_Ael(const UnitCell &cell, ModulePW::PW_Basis *rho_basis)
{
    double Ael = 0.0;
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        Ael -= TOTN_real[ir] * delta_phi[ir];
    }
    Parallel_Reduce::reduce_double_pool(Ael);
    Ael = Ael * cell.omega / rho_basis->nxyz;
    // cout << "Ael: " << Ael << endl;
    return Ael;
}

double surchem::cal_Acav(const UnitCell &cell, ModulePW::PW_Basis *rho_basis)
{
    double Acav = 0.0;
    Acav = GlobalV::tau * qs;
    Acav = Acav * cell.omega / rho_basis->nxyz; // unit Ry
    Parallel_Reduce::reduce_double_pool(Acav);
    // cout << "Acav: " << Acav << endl;
    return Acav;
}