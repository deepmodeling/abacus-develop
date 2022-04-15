#include "surchem.h"

void surchem::cal_totn(const UnitCell &cell, PW_Basis &pwb,
                       const complex<double> *Porter_g, complex<double> *N,
                       complex<double> *TOTN) {

    // vloc to N
    complex<double> *vloc_g = new complex<double>[pwb.ngmc];
    ModuleBase::GlobalFunc::ZEROS(vloc_g, pwb.ngmc);

    // method 1
    GlobalC::UFFT.ToReciSpace(GlobalC::pot.vltot,
                              vloc_g); // now n is vloc in Recispace
    for (int ig = pwb.gstart; ig < pwb.ngmc; ig++) {
        if (pwb.gg[ig] >= 1.0e-12) // LiuXh 20180410
        {
            const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI /
                               (cell.tpiba2 * pwb.gg[ig]);

            N[ig] = -vloc_g[ig] / fac;
        }
    }

    if (GlobalV::MY_RANK == 0) {
        N[0] = Porter_g[0];
    }

    for (int ig = 0; ig < pwb.ngmc; ig++) {
        TOTN[ig] = N[ig] - Porter_g[ig];
    }
    delete[] vloc_g;
    return;
}