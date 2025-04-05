#ifndef WRITE_DOS_LCAO_H
#define WRITE_DOS_LCAO_H

#include "module_base/matrix.h" // use matrix
#include "module_cell/klist.h"  // use K_Vectors
#include "module_psi/psi.h"     // use psi::Psi<T>
#include "module_hamilt_general/hamilt.h" // use hamilt::Hamilt<T>
#include "module_basis/module_ao/parallel_orbitals.h" // use Parallel_Orbitals

namespace ModuleIO
{
	/// @brief calculate density of states(DOS), 
    /// partial density of states(PDOS),
    ///  and mulliken charge for LCAO base
    template <typename T>
    void write_dos_lcao(
        const UnitCell& ucell,
        const psi::Psi<T>* psi,
        const Parallel_Orbitals &pv, 
        const ModuleBase::matrix& ekb,
        const ModuleBase::matrix& wg,
        const double& dos_edelta_ev,
        const double& dos_scale,
        const double& bcoeff,
		const K_Vectors& kv,
		const int nbands,
		const elecstate::efermi &energy_fermi,
		hamilt::Hamilt<T>* p_ham,
        std::ofstream &ofs_running);
}
#endif
