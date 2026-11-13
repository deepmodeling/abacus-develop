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
} // namespace

TEST(BasisFlattenedIndexingTest, FlattenAtomOrderAndSpin)
{
    const std::vector<int> atom_types = {1, 0, 1};
    const BasisFlattenedIndexing indexing(test_hierarchy(), atom_types, 2);

    EXPECT_EQ(indexing.nat(), 3);
    EXPECT_EQ(indexing.npol(), 2);
    EXPECT_EQ(indexing.size(), 14);

    const BasisRange first = indexing.range(0);
    const BasisRange second = indexing.range(1);
    const BasisRange third = indexing.range(2);
    EXPECT_EQ(first.ibas0, 0);
    EXPECT_EQ(first.ibas1, 2);
    EXPECT_EQ(second.ibas0, 2);
    EXPECT_EQ(second.ibas1, 12);
    EXPECT_EQ(third.ibas0, 12);
    EXPECT_EQ(third.ibas1, 14);

    EXPECT_EQ(indexing.atom_of(0), 0);
    EXPECT_EQ(indexing.atom_of(1), 0);
    EXPECT_EQ(indexing.atom_of(11), 1);
    EXPECT_EQ(indexing.atom_of(13), 2);
}

TEST(BasisFlattenedIndexingTest, MapHierarchicalCoordinates)
{
    const std::vector<int> atom_types = {0, 1};
    const BasisFlattenedIndexing indexing(test_hierarchy(), atom_types, 2);

    EXPECT_EQ(indexing.ibas(0, 1, 0, 2, 1), 9);
    EXPECT_EQ(indexing.ibas(1, 0, 0, 0, 1), 11);

    const BasisLabel& label = indexing.label(9);
    EXPECT_EQ(label.iat, 0);
    EXPECT_EQ(label.itype, 0);
    EXPECT_EQ(label.l, 1);
    EXPECT_EQ(label.izeta, 0);
    EXPECT_EQ(label.im, 2);
    EXPECT_EQ(label.m, -1);
    EXPECT_EQ(label.spin, 1);
}

TEST(BasisFlattenedIndexingTest, RejectInvalidCoordinates)
{
    const std::vector<int> atom_types = {0};
    const BasisFlattenedIndexing indexing(test_hierarchy(), atom_types, 1);

    EXPECT_THROW(indexing.range(-1), std::out_of_range);
    EXPECT_THROW(indexing.range(1), std::out_of_range);
    EXPECT_THROW(indexing.atom_of(7), std::out_of_range);
    EXPECT_THROW(indexing.ibas(0, 0, 0, 0, 1), std::out_of_range);

    const std::vector<int> invalid_atom_types = {2};
    EXPECT_THROW(BasisFlattenedIndexing(test_hierarchy(), invalid_atom_types, 1),
                 std::out_of_range);
    EXPECT_THROW(BasisFlattenedIndexing(test_hierarchy(), atom_types, 0),
                 std::invalid_argument);

    ModuleBase::Element_Basis_Index::Range invalid = test_hierarchy_range();
    invalid[0][0].M = 2;
    const ModuleBase::Element_Basis_Index::IndexLNM invalid_hierarchy =
        ModuleBase::Element_Basis_Index::construct_index(invalid);
    EXPECT_THROW(BasisFlattenedIndexing(invalid_hierarchy, atom_types, 1),
                 std::invalid_argument);
}
