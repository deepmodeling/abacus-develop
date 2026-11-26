#ifndef CAL_PDOS_GAMMA_H
#define CAL_PDOS_GAMMA_H

#include "source_base/matrix.h"
#include "source_cell/klist.h"  // use K_Vectors
#include "source_psi/psi.h"     // use psi::Psi<T>
#include "source_hamilt/hamilt.h" // use hamilt::Hamilt<T>
#include "source_basis/module_ao/parallel_orbitals.h" // use Parallel_Orbitals
#include <vector>

namespace ModuleIO
{

    std::vector<ModuleBase::matrix> cal_pdos(
            const psi::Psi<double>* psi,
            hamilt::Hamilt<double>* p_ham,
            const Parallel_Orbitals& pv,
            const K_Vectors& kv,
            const int nspin0,
            const int nbands,
            const ModuleBase::matrix& ekb,
            const double& emax,
            const double& emin,
            const double& dos_edelta_ev,
            const double& bcoeff,
            const int nlocal,
            const int nspin);

}

#endif
