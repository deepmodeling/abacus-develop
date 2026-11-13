#include "source_base/basis_flattened_indexing.h"

#include <stdexcept>

namespace
{
int magnetic_quantum_number(const int im)
{
    return im % 2 == 0 ? -im / 2 : (im + 1) / 2;
}
} // namespace

BasisFlattenedIndexing::BasisFlattenedIndexing(
    const ModuleBase::Element_Basis_Index::IndexLNM& hierarchy,
    const std::vector<int>& atom_types,
    const int npol)
    : nat_(static_cast<int>(atom_types.size())), npol_(npol)
{
    if (npol_ <= 0)
    {
        throw std::invalid_argument("BasisFlattenedIndexing requires positive npol");
    }

    atom_ranges_.reserve(atom_types.size());
    for (int iat = 0; iat < nat_; ++iat)
    {
        const std::size_t ibas0 = labels_.size();
        const int itype = atom_types[iat];
        if (itype < 0 || static_cast<std::size_t>(itype) >= hierarchy.size())
        {
            throw std::out_of_range("BasisFlattenedIndexing atom type is out of range");
        }

        for (std::size_t l = 0; l < hierarchy[itype].size(); ++l)
        {
            const ModuleBase::Element_Basis_Index::Index_TL& radial_index = hierarchy[itype][l];
            if (radial_index.M != 2 * l + 1)
            {
                throw std::invalid_argument(
                    "BasisFlattenedIndexing M axis does not contain 2*l+1 components");
            }
            if (radial_index.size() != radial_index.N)
            {
                throw std::invalid_argument(
                    "BasisFlattenedIndexing N axis does not match the radial index");
            }

            for (std::size_t izeta = 0; izeta < radial_index.N; ++izeta)
            {
                if (radial_index[izeta].size() != radial_index.M)
                {
                    throw std::invalid_argument(
                        "BasisFlattenedIndexing M axis size is inconsistent");
                }
                for (std::size_t im = 0; im < radial_index.M; ++im)
                {
                    const int m = magnetic_quantum_number(static_cast<int>(im));
                    for (int spin = 0; spin < npol_; ++spin)
                    {
                        const BasisLabel label = {iat,
                                                  itype,
                                                  static_cast<int>(l),
                                                  static_cast<int>(izeta),
                                                  static_cast<int>(im),
                                                  m,
                                                  spin};
                        labels_.push_back(label);
                    }
                }
            }
        }
        const BasisRange atom_range = {ibas0, labels_.size()};
        atom_ranges_.push_back(atom_range);
    }
}

BasisRange BasisFlattenedIndexing::range(const int iat) const
{
    if (iat < 0 || iat >= nat_)
    {
        throw std::out_of_range("BasisFlattenedIndexing atom index is out of range");
    }
    return atom_ranges_[iat];
}

int BasisFlattenedIndexing::atom_of(const std::size_t ibas) const
{
    return labels_.at(ibas).iat;
}

const BasisLabel& BasisFlattenedIndexing::label(const std::size_t ibas) const
{
    return labels_.at(ibas);
}

std::size_t BasisFlattenedIndexing::ibas(const int iat,
                                         const int l,
                                         const int izeta,
                                         const int im,
                                         const int spin) const
{
    const BasisRange atom_range = range(iat);
    for (std::size_t index = atom_range.ibas0; index < atom_range.ibas1; ++index)
    {
        const BasisLabel& basis = labels_[index];
        if (basis.l == l && basis.izeta == izeta && basis.im == im && basis.spin == spin)
        {
            return index;
        }
    }
    throw std::out_of_range("BasisFlattenedIndexing coordinate is out of range");
}
