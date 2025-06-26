#ifndef WRITE_PROJ_BAND_LCAO_H
#define WRITE_PROJ_BAND_LCAO_H
#include "source_basis/module_ao/ORB_read.h"
#include "source_cell/klist.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_estate/elecstate.h"
#include "source_psi/psi.h"
#include "source_hamilt/hamilt.h"
#include "source_basis/module_ao/parallel_orbitals.h"

namespace ModuleIO
{
    template <typename TK>
    void write_proj_band_lcao(
        const psi::Psi<TK>* psi,
		const Parallel_Orbitals &pv,
		const elecstate::ElecState* pelec,
		const K_Vectors& kv,
		const UnitCell &ucell, 
        hamilt::Hamilt<TK>* p_ham);
}

#endif
