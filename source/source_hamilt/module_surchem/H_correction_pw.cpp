#include <cmath>

#include "source_base/constants.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "surchem.h"

ModuleBase::matrix surchem::v_correction(const UnitCell& cell,
                                         const Parallel_Grid& pgrid,
                                         const ModulePW::PW_Basis* rho_basis,
                                         const int& nspin,
                                         const double* const* const rho,
                                         const double* vlocal,
                                         Structure_Factor* sf)
{
    ModuleBase::TITLE("surchem", "v_cor");
    ModuleBase::timer::tick("surchem", "v_cor");

    assert(rho_basis->nrxx>0);
   
    double* porter = new double[rho_basis->nrxx];
	for (int i = 0; i < rho_basis->nrxx; i++)
	{
		porter[i] = 0.0;
	}
    const int nspin0 = (nspin == 2) ? 2 : 1;
	for (int is = 0; is < nspin0; is++)
	{
		for (int ir = 0; ir < rho_basis->nrxx; ir++)
		{
			porter[ir] += rho[is][ir];
		}
	}

    complex<double>* porter_g = new complex<double>[rho_basis->npw];
    ModuleBase::GlobalFunc::ZEROS(porter_g, rho_basis->npw);

    rho_basis->real2recip(porter, porter_g);

    complex<double>* n = new complex<double>[rho_basis->npw];
    complex<double>* total_n = new complex<double>[rho_basis->npw];
    complex<double>* ps_totn = new complex<double>[rho_basis->npw];

    cal_totn(cell, rho_basis, porter_g, n, total_n, vlocal);

    cal_pseudo(cell, pgrid, rho_basis, porter_g, ps_totn, sf);

    ModuleBase::matrix v(nspin, rho_basis->nrxx);

    v += cal_vel(cell, rho_basis, total_n, ps_totn, nspin);
    v += cal_vcav(cell, rho_basis, ps_totn, nspin);

    delete[] porter;
    delete[] porter_g;
    delete[] n;
    delete[] ps_totn;
    delete[] total_n;

    ModuleBase::timer::tick("surchem", "v_cor");
    return v;
}
