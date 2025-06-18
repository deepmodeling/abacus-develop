#ifndef EXX_ABFS_ABFS_INDEX_H
#define EXX_ABFS_ABFS_INDEX_H

#include "../module_basis/module_ao/ORB_atomic_lm.h"
#include "../source_base/element_basis_index.h"
#include "exx_abfs.h"

#include <vector>

class LCAO_Orbitals;

class Exx_Abfs::Abfs_Index
{
  public:
    static ModuleBase::Element_Basis_Index::Range construct_range(const LCAO_Orbitals& orb);
    static ModuleBase::Element_Basis_Index::Range construct_range(
        const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb);
};

#endif // EXX_ABFS_ABFS_INDEX_H