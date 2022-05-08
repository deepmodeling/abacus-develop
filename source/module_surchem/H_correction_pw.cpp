// #include "../src_pw/diago_cg.h"
#include "../module_base/constants.h"
#include "../module_base/timer.h"
#include "../module_xc/xc_functional.h"
#include "../src_parallel/parallel_reduce.h"
#include "surchem.h"

#include <cmath>

ModuleBase::matrix surchem::v_correction(const UnitCell &cell,
                                         PW_Basis &pwb,
                                         const int &nspin,
                                         const double *const *const rho)
{
    ModuleBase::TITLE("surchem", "v_correction");
    ModuleBase::timer::tick("surchem", "v_correction");

    double *Porter = new double[pwb.nrxx];
    for (int i = 0; i < pwb.nrxx; i++)
        Porter[i] = 0.0;
    const int nspin0 = (nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin0; is++)
        for (int ir = 0; ir < pwb.nrxx; ir++)
            Porter[ir] += rho[is][ir];

    complex<double> *Porter_g = new complex<double>[pwb.ngmc];
    ModuleBase::GlobalFunc::ZEROS(Porter_g, pwb.ngmc);

    GlobalC::UFFT.ToReciSpace(Porter, Porter_g);

    complex<double> *N = new complex<double>[pwb.ngmc];
    complex<double> *TOTN = new complex<double>[pwb.ngmc];
    complex<double> *PS_TOTN = new complex<double>[pwb.ngmc];

    cal_totn(cell, pwb, Porter_g, N, TOTN);

    cal_pseudo(cell, pwb, Porter_g, PS_TOTN);

    ModuleBase::matrix v(nspin, pwb.nrxx);

    v += cal_vel(cell, pwb, TOTN, PS_TOTN, nspin);
    v += cal_vcav(cell, pwb, PS_TOTN, nspin);

    delete[] Porter;
    delete[] Porter_g;
    delete[] N;
    delete[] PS_TOTN;
    delete[] TOTN;

    ModuleBase::timer::tick("surchem", "v_correction");
    return v;
}