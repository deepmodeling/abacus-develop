#include "source_base/atom_basis_view.h"

#include "source_base/basis_flattened_indexing.h"
#include "source_base/block_cyclic_distribution.h"

#include <stdexcept>

AtomBasisView::AtomBasisView(const BlockCyclicDistribution* layout,
                             const int ibas0,
                             const int local0,
                             const int nloc,
                             const bool col_axis)
    : layout_(layout), ibas0_(ibas0), local0_(local0), nloc_(nloc), col_axis_(col_axis)
{
}

AtomBasisView AtomBasisView::make(const BasisFlattenedIndexing& basis,
                                  const BlockCyclicDistribution& layout,
                                  const int iat,
                                  const bool columns)
{
    const int global_size = columns ? layout.ncglb() : layout.nrglb();
    if (basis.size() != static_cast<std::size_t>(global_size))
    {
        throw std::invalid_argument("basis indexing and matrix distribution sizes disagree");
    }

    const BasisRange basis_range = basis.range(iat);
    const int ibas0 = static_cast<int>(basis_range.ibas0);
    int local0 = -1;
    int last_local = -2;
    int nloc = 0;
    for (std::size_t ibas = basis_range.ibas0; ibas < basis_range.ibas1; ++ibas)
    {
        const int ibas_global = static_cast<int>(ibas);
        const int local = columns ? layout.icloc(ibas_global) : layout.irloc(ibas_global);
        if (local < 0)
        {
            continue;
        }
        if (nloc == 0)
        {
            local0 = local;
        }
        else if (local != last_local + 1)
        {
            throw std::invalid_argument("atom basis entries are not contiguous");
        }
        last_local = local;
        ++nloc;
    }
    return AtomBasisView(&layout, ibas0, local0, nloc, columns);
}

int AtomBasisView::operator[](const std::size_t index) const
{
    // Reproduces the eager value the old vector stored: the global basis index
    // of the owned local entry minus the atom's first basis index (ibas - ibas0).
    const int local = local0_ + static_cast<int>(index);
    const int ibas_global = col_axis_ ? layout_->icglb(local) : layout_->irglb(local);
    const int entry = ibas_global - ibas0_;
    return entry;
}

AtomBasisView AtomBasisView::rows(const BasisFlattenedIndexing& basis,
                                  const BlockCyclicDistribution& layout,
                                  const int iat)
{
    return make(basis, layout, iat, false);
}

AtomBasisView AtomBasisView::cols(const BasisFlattenedIndexing& basis,
                                  const BlockCyclicDistribution& layout,
                                  const int iat)
{
    return make(basis, layout, iat, true);
}
