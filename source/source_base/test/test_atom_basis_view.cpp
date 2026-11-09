#include "source_base/atom_basis_view.h"
#include "source_base/basis_flattened_indexing.h"
#include "source_base/parallel_blkcyc.h"

#include "gtest/gtest.h"

#include <stdexcept>
#include <vector>

namespace
{
ModuleBase::Element_Basis_Index::Range test_hierarchy_range()
{
    ModuleBase::Element_Basis_Index::Range range(2);
    range[0].resize(2);
    range[0][0].N = 2;
    range[0][0].M = 1;
    range[0][1].N = 1;
    range[0][1].M = 3;
    range[1].resize(1);
    range[1][0].N = 1;
    range[1][0].M = 1;
    return range;
}

ModuleBase::Element_Basis_Index::IndexLNM test_hierarchy()
{
    const ModuleBase::Element_Basis_Index::Range range = test_hierarchy_range();
    return ModuleBase::Element_Basis_Index::construct_index(range);
}

// Builds the entry list from the lazy view, one index at a time.
std::vector<int> collect_entries(const AtomBasisView& view)
{
    std::vector<int> entries;
    for (int k = 0; k < view.nloc(); ++k)
    {
        entries.push_back(view[k]);
    }
    return entries;
}
} // namespace

TEST(AtomBasisViewTest, DeriveSerialRowsAndCols)
{
    const std::vector<int> atom_types = {1, 0};
    const BasisFlattenedIndexing basis(test_hierarchy(), atom_types, 1);
    Parallel_Blkcyclic layout;
    layout.set_serial(6, 6);

    const AtomBasisView first_rows = AtomBasisView::rows(basis, layout, 0);
    const AtomBasisView second_rows = AtomBasisView::rows(basis, layout, 1);
    const AtomBasisView second_cols = AtomBasisView::cols(basis, layout, 1);

    EXPECT_EQ(first_rows.local0(), 0);
    EXPECT_EQ(collect_entries(first_rows), std::vector<int>({0}));
    EXPECT_EQ(second_rows.local0(), 1);
    EXPECT_EQ(collect_entries(second_rows), std::vector<int>({0, 1, 2, 3, 4}));
    EXPECT_EQ(second_cols.local0(), 1);
    EXPECT_EQ(collect_entries(second_cols), std::vector<int>({0, 1, 2, 3, 4}));
}

TEST(AtomBasisViewTest, RejectLayoutSizeMismatch)
{
    const std::vector<int> atom_types = {0};
    const BasisFlattenedIndexing basis(test_hierarchy(), atom_types, 1);
    Parallel_Blkcyclic layout;
    layout.set_serial(6, 6);

    EXPECT_THROW(AtomBasisView::rows(basis, layout, 0), std::invalid_argument);
    EXPECT_THROW(AtomBasisView::cols(basis, layout, 0), std::invalid_argument);
}

TEST(AtomBasisViewTest, RejectAtomOutOfRange)
{
    const std::vector<int> atom_types = {0};
    const BasisFlattenedIndexing basis(test_hierarchy(), atom_types, 1);
    Parallel_Blkcyclic layout;
    layout.set_serial(5, 5);

    EXPECT_THROW(AtomBasisView::rows(basis, layout, 1), std::out_of_range);
    EXPECT_THROW(AtomBasisView::cols(basis, layout, 1), std::out_of_range);
}
