#ifndef CAL_PDOS_MULTIK_H
#define CAL_PDOS_MULTIK_H

#include "source_base/matrix.h"
#include "source_cell/klist.h"  // use K_Vectors
#include "source_psi/psi.h"     // use psi::Psi<T>
#include "source_hamilt/hamilt.h" // use hamilt::Hamilt<T>
#include "source_basis/module_ao/parallel_orbitals.h" // use Parallel_Orbitals
#include <string>

namespace ModuleIO
{

    // pdos for multi-k point
	void cal_pdos(
			const psi::Psi<std::complex<double>>* psi,
			hamilt::Hamilt<std::complex<double>>* p_ham,
			const Parallel_Orbitals& pv,
			const UnitCell& ucell,
			const K_Vectors& kv,
			const int nspin0,
			const int nbands,
			const ModuleBase::matrix& ekb,
			const double& emax,
			const double& emin,
			const double& dos_edelta_ev,
		const double& bcoeff,
		const int istep,
		const int nlocal,
		const int nspin,
		const std::string& global_out_dir);

}

#endif 
