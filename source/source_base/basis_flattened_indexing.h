#ifndef ABACUS_SOURCE_SOURCE_BASE_BASIS_FLATTENED_INDEXING_H
#define ABACUS_SOURCE_SOURCE_BASE_BASIS_FLATTENED_INDEXING_H

#include "source_base/element_basis_index.h"

#include <cstddef>
#include <vector>

/**
 * @brief Half-open range of flattened basis indices belonging to one atom.
 */
struct BasisRange
{
    std::size_t ibas0;
    std::size_t ibas1;
};

/**
 * @brief Hierarchical coordinates of one flattened basis index.
 */
struct BasisLabel
{
    int iat;
    int itype;
    int l;
    int izeta;
    int im;
    int m;
    int spin;
};

/**
 * @brief Canonical atom-major flattening of the LCAO basis hierarchy.
 *
 * The flattened order is atom, angular momentum, radial function, magnetic
 * component, then spin component. This value defines global basis indexing
 * only; it contains no MPI, matrix-distribution, or descriptor state.
 */
class BasisFlattenedIndexing
{
  public:
    BasisFlattenedIndexing(const ModuleBase::Element_Basis_Index::IndexLNM& hierarchy,
                           const std::vector<int>& atom_types,
                           const int npol);

    std::size_t size() const
    {
        return labels_.size();
    }

    int nat() const
    {
        return nat_;
    }

    int npol() const
    {
        return npol_;
    }

    BasisRange range(const int iat) const;
    int atom_of(const std::size_t ibas) const;
    const BasisLabel& label(const std::size_t ibas) const;
    std::size_t ibas(const int iat,
                     const int l,
                     const int izeta,
                     const int im,
                     const int spin) const;

  private:
    std::vector<BasisRange> atom_ranges_;
    std::vector<BasisLabel> labels_;
    int nat_ = 0;
    int npol_ = 1;
};

#endif
