#ifndef DFTU_NAO_FS_ACCUM_H
#define DFTU_NAO_FS_ACCUM_H

/// @file dftu_nao_fs_accum.h
/// @brief Diagonal-block accumulation helpers for DFT+U force/stress.
///
/// These templates are header-only so they can be unit-tested without
/// pulling in the full DftuFsEnv closure (ScalapackConnector, folding,
/// Plus_U_Base). accumulate_onsite_force stays in dftu_nao_fs_k.cpp
/// because it needs the Plus_U_Base occupation-matrix lookup.

#include "source_base/matrix.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/unitcell.h"

#include <cassert>
#include <complex>

namespace DFTU_LCAO
{

/// @brief Add the real part of diagonal local-block entries to one force component.
///
/// Sums dm(ir, ic) over local block pairs whose global orbital indices
/// coincide, attributing each entry to the atom owning the orbital along
/// Cartesian component dim.
template <typename T>
void accumulate_diag_force(const Parallel_Orbitals& pv,
                           const UnitCell& ucell,
                           const T* dm,
                           const int dim,
                           ModuleBase::matrix& force_dftu)
{
    assert(dm != nullptr);
    assert(dim >= 0 && dim < 3);
    for (int ir = 0; ir < pv.nrow; ir++)
    {
        const int iwt1 = pv.local2global_row(ir);
        const int iat1 = ucell.iwt2iat[iwt1];
        for (int ic = 0; ic < pv.ncol; ic++)
        {
            if (pv.local2global_col(ic) == iwt1)
            {
                force_dftu(iat1, dim) += std::real(dm[ic * pv.nrow + ir]);
            }
        }
    }
}

/// @brief Add the real part of diagonal local-block entries to one stress pair.
template <typename T>
void accumulate_diag_stress(const Parallel_Orbitals& pv,
                            const T* dm,
                            const int dim1,
                            const int dim2,
                            const double factor,
                            ModuleBase::matrix& stress_dftu)
{
    assert(dm != nullptr);
    assert(dim1 >= 0 && dim1 < 3);
    assert(dim2 >= 0 && dim2 < 3);
    for (int ir = 0; ir < pv.nrow; ir++)
    {
        const int iwt1 = pv.local2global_row(ir);
        for (int ic = 0; ic < pv.ncol; ic++)
        {
            if (pv.local2global_col(ic) == iwt1)
            {
                stress_dftu(dim1, dim2) += factor * std::real(dm[ic * pv.nrow + ir]);
            }
        }
    }
}

} // namespace DFTU_LCAO

#endif // DFTU_NAO_FS_ACCUM_H
