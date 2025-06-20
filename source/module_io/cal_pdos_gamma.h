#ifndef CAL_PDOS_GAMMA_H
#define CAL_PDOS_GAMMA_H

#include "source_base/matrix.h"
#include "source_cell/klist.h"  // use K_Vectors
#include "module_psi/psi.h"     // use psi::Psi<T>
#include "module_hamilt_general/hamilt.h" // use hamilt::Hamilt<T>
#include "source_basis/module_ao/parallel_orbitals.h" // use Parallel_Orbitals

namespace ModuleIO
{

	void cal_pdos(
			const psi::Psi<double>* psi,
			hamilt::Hamilt<double>* p_ham,
			const Parallel_Orbitals& pv,
			const UnitCell& ucell,
			const K_Vectors& kv,
			const int nspin0,
			const int nbands,
			const ModuleBase::matrix& ekb,
			const double& emax,
			const double& emin,
			const double& dos_edelta_ev,
			const double& bcoeff);

	void print_tdos_gamma(
			const ModuleBase::matrix* pdos,
			const int nlocal,
			const int npoints,
			const double& emin,
			const double& dos_edelta_ev);

	void print_pdos_gamma(
			const UnitCell& ucell,
			const ModuleBase::matrix* pdos,
			const int nlocal,
			const int npoints,
			const double& emin,
			const double& dos_edelta_ev);

}

#endif 
