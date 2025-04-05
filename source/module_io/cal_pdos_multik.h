#ifndef PDOS_H
#define PDOS_H

#include "module_base/matrix.h"

namespace ModuleIO
{

	void cal_pdos_multik(
			const int nspin0,
			const double& emax,
			const double& emin,
			const double& dos_edelta_ev,
			const double& bcoeff,
            hamilt::Hamilt<std::complex<double>>* p_ham,
			const psi::Psi<std::complex<double>>* psi,
			const Parallel_Orbitals& pv);

	void print_tdos_multik(
			const ModuleBase::matrix& pdos,
			const int nlocal,
			const int npoints,
			const double& emin,
			const double& dos_edelta_ev);

	void print_pdos_multik(
			const UnitCell& ucell,
			const ModuleBase::matrix& pdos,
			const int nlocal,
			const int npoints,
			const double& emin,
			const double& dos_edelta_ev);

}

#endif 
