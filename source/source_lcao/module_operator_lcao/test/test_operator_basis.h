#ifndef TEST_OPERATOR_BASIS_H
#define TEST_OPERATOR_BASIS_H

#include "source_base/basis_flattened_indexing.h"
#include "source_cell/unitcell.h"

#include <cstddef>

namespace operator_test
{

// The operator fixtures use repeated pure-s shells with one radial function.
inline BasisFlattenedIndexing make_s_orbital_basis(const UnitCell& cell)
{
    const std::size_t norbital = static_cast<std::size_t>(cell.atoms[0].nw);
    ModuleBase::Element_Basis_Index::Range range(cell.ntype);
    for (int it = 0; it < cell.ntype; ++it)
    {
        range[it].resize(1);
        range[it][0].N = norbital;
        range[it][0].M = 1;
    }
    const ModuleBase::Element_Basis_Index::IndexLNM hierarchy =
        ModuleBase::Element_Basis_Index::construct_index(range);
    return BasisFlattenedIndexing(hierarchy, cell.iat2it, cell.get_npol());
}

} // namespace operator_test

#endif
