#include "surchem.h"

double surchem::cal_Ael(const UnitCell &cell, PW_Basis &pwb, const double *TOTN_real, const double *delta_phi_R)
{
    double Ael = 0.0;
    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        Ael -= TOTN_real[ir] * delta_phi_R[ir];
    }
    Parallel_Reduce::reduce_double_pool(Ael);
    Ael = Ael * cell.omega * 0.5 / pwb.nxyz;
    cout << "Ael: " << Ael << endl;
    return Ael;
}

double surchem::cal_Acav(const UnitCell &cell, PW_Basis &pwb, double qs)
{
    double Acav = 0.0;
    Acav = GlobalV::tau * qs;
    Acav = Acav * cell.omega * 0.5 / pwb.nxyz;
    Parallel_Reduce::reduce_double_pool(Acav);
    cout << "Acav: " << Acav << endl;
    return Acav;
}