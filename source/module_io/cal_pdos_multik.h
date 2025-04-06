#ifndef CAL_PDOS_MULTIK_H
#define CAL_PDOS_MULTIK_H

#include "module_base/matrix.h"
#include "module_cell/klist.h"  // use K_Vectors
#include "module_psi/psi.h"     // use psi::Psi<T>
#include "module_hamilt_general/hamilt.h" // use hamilt::Hamilt<T>
#include "module_basis/module_ao/parallel_orbitals.h" // use Parallel_Orbitals

namespace ModuleIO
{

	void cal_pdos(
			const UnitCell& ucell,
			const int nspin0,
            const int nbands,
			const double& emax,
			const double& emin,
			const double& dos_edelta_ev,
			const double& bcoeff,
            hamilt::Hamilt<std::complex<double>>* p_ham,
			const psi::Psi<std::complex<double>>* psi,
			const Parallel_Orbitals& pv);

	void print_tdos_multik(
			const ModuleBase::matrix* pdos,
			const int nlocal,
			const int npoints,
			const double& emin,
			const double& dos_edelta_ev);

	void print_pdos_multik(
			const UnitCell& ucell,
			const ModuleBase::matrix* pdos,
			const int nlocal,
			const int npoints,
			const double& emin,
			const double& dos_edelta_ev);

}

#endif 
