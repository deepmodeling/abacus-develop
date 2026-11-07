#ifndef ABACUS_SOURCE_BASE_ATOM_BASIS_VIEW_H
#define ABACUS_SOURCE_BASE_ATOM_BASIS_VIEW_H

#include <cstddef>

class BasisFlattenedIndexing;
class BlockCyclicDistribution;

/**
 * @brief Atom-local basis entries derived from basis indexing and distribution.
 *
 * Entries are atom-local basis numbers in local matrix order, materialized on
 * demand from the indexing and distribution inputs; the view allocates nothing.
 */
class AtomBasisView
{
  public:
    static AtomBasisView rows(const BasisFlattenedIndexing& basis,
                              const BlockCyclicDistribution& layout,
                              int iat);
    static AtomBasisView cols(const BasisFlattenedIndexing& basis,
                              const BlockCyclicDistribution& layout,
                              int iat);

    int nloc() const
    {
        return nloc_;
    }

    bool empty() const
    {
        return nloc_ == 0;
    }

    int local0() const
    {
        return local0_;
    }

    // index must be smaller than nloc(); out-of-range access is not checked.
    // This view borrows the distribution it was derived from; it must not
    // outlive it.
    int operator[](const std::size_t index) const;

  private:
    AtomBasisView(const BlockCyclicDistribution* layout,
                  const int ibas0,
                  const int local0,
                  const int nloc,
                  const bool col_axis);

    static AtomBasisView make(const BasisFlattenedIndexing& basis,
                              const BlockCyclicDistribution& layout,
                              const int iat,
                              const bool columns);

    const BlockCyclicDistribution* layout_ = nullptr;
    int ibas0_ = 0;
    int local0_ = 0;
    int nloc_ = 0;
    bool col_axis_ = false;
};

#endif
