#ifndef WRITE_PROJ_BAND_LCAO_H
#define WRITE_PROJ_BAND_LCAO_H
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_estate/elecstate.h"
#include "source_psi/psi.h"
#include "source_hamilt/hamilt.h"
#include "source_basis/module_ao/parallel_orbitals.h"

namespace ModuleIO
{
    /// @param ucell Unit cell
    /// @param p_ham Hamiltonian used to build the projection overlaps
    /// @param istep Zero-based ionic step index; written one-based
    /// @param out_app_flag Output mode: 2 writes pbands{spin}g{step}.xml
    template <typename TK>
    void write_proj_band_lcao(
        const psi::Psi<TK>* psi,
        const Parallel_Orbitals &pv,
        const elecstate::ElecState* pelec,
        const K_Vectors& kv,
        const UnitCell &ucell, 
        hamilt::Hamilt<TK>* p_ham,
        const int istep,
        const int out_app_flag);
}

#endif
