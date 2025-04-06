#ifndef PDOS_H
#define PDOS_H

#include "module_base/matrix.h"

namespace ModuleIO
{

	void cal_pdos_gamma(
			const int nspin0,
			const double& emax,
			const double& emin,
			const double& dos_edelta_ev,
			const double& bcoeff,
			const double* sk,
			const psi::Psi<double>* psi,
			const Parallel_Orbitals& pv);

	void print_tdos_gamma(
			const ModuleBase::matrix& pdos,
			const int nlocal,
			const int npoints,
			const double& emin,
			const double& dos_edelta_ev);

	void print_pdos_gamma(
			const UnitCell& ucell,
			const ModuleBase::matrix& pdos,
			const int nlocal,
			const int npoints,
			const double& emin,
			const double& dos_edelta_ev);

}

#endif 
