#include "surchem.h"

double surchem::cal_Ael(const UnitCell &cell, const int& nrxx, const int& nxyz)
{
    double Ael = 0.0;
    for (int ir = 0; ir < nrxx; ir++)
    {
        Ael -= TOTN_real[ir] * delta_phi[ir];
    }
    Parallel_Reduce::reduce_pool(Ael);
    Ael = Ael * cell.omega / nxyz;
    // cout << "Ael: " << Ael << endl;
    return Ael;
}

double surchem::cal_Acav(const UnitCell &cell, const int& nxyz)
{
    double Acav = 0.0;
    Acav = GlobalV::tau * qs;
    Acav = Acav * cell.omega / nxyz; // unit Ry
    Parallel_Reduce::reduce_pool(Acav);
    // cout << "Acav: " << Acav << endl;
    return Acav;
}