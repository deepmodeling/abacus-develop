// #include "../src_pw/diago_cg.h"
#include "../module_base/constants.h"
#include "../module_base/timer.h"
#include "../module_xc/xc_functional.h"
#include "../src_parallel/parallel_reduce.h"
#include "surchem.h"

#include <cmath>

ModuleBase::matrix surchem::v_correction(const UnitCell &cell,
                                         ModulePW::PW_Basis *rho_basis,
                                         const int &nspin,
                                         const double *const *const rho)
{
    ModuleBase::TITLE("surchem", "v_correction");
    ModuleBase::timer::tick("surchem", "v_correction");

    double *Porter = new double[rho_basis->nrxx];
    for (int i = 0; i < rho_basis->nrxx; i++)
        Porter[i] = 0.0;
    const int nspin0 = (nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin0; is++)
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
            Porter[ir] += rho[is][ir];

    complex<double> *Porter_g = new complex<double>[rho_basis->npw];
    ModuleBase::GlobalFunc::ZEROS(Porter_g, rho_basis->npw);

    rho_basis->real2recip(Porter, Porter_g);

    complex<double> *N = new complex<double>[rho_basis->npw];
    complex<double> *TOTN = new complex<double>[rho_basis->npw];
    complex<double> *PS_TOTN = new complex<double>[rho_basis->npw];

    cal_totn(cell, rho_basis, Porter_g, N, TOTN);

    cal_pseudo(cell, rho_basis, Porter_g, PS_TOTN);

    ModuleBase::matrix v(nspin, rho_basis->nrxx);

    v += cal_vel(cell, rho_basis, TOTN, PS_TOTN, nspin);
    v += cal_vcav(cell, rho_basis, PS_TOTN, nspin);

    delete[] Porter;
    delete[] Porter_g;
    delete[] N;
    delete[] PS_TOTN;
    delete[] TOTN;

    ModuleBase::timer::tick("surchem", "v_correction");
    return v;
}