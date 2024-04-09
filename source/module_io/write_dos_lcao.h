#ifndef WRITE_DOS_LCAO_H
#define WRITE_DOS_LCAO_H

#include "module_base/matrix.h"
#include "module_cell/klist.h"
#include "module_psi/psi.h"
#include "module_hamilt_general/hamilt.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"

namespace ModuleIO
{
	/// @brief calculate density of states(DOS) and partial density of states(PDOS) and mulliken charge for LCAO base
    template <typename T>
    void write_dos_lcao(
        const psi::Psi<T>* psi,
        LCAO_Matrix &lm,
        const Parallel_Orbitals &pv, 
        const ModuleBase::matrix& ekb,
        const ModuleBase::matrix& wg,
        const double& dos_edelta_ev,
        const double& dos_scale,
        const double& bcoeff,
        const K_Vectors& kv,
        hamilt::Hamilt<T>* p_ham);
}
#endif
