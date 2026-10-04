#include "source_cell/cell_tools.h"

#include <gtest/gtest.h>

TEST(CellTools, ExtractsPositionsMassesAndValenceChargesInAtomOrder)
{
    Atom atoms[2];
    atoms[0].na = 2;
    atoms[0].mass = 12.0;
    atoms[0].ncpp.zv = 4.0;
    atoms[0].tau = {ModuleBase::Vector3<double>(0.1, 0.2, 0.3),
                    ModuleBase::Vector3<double>(0.4, 0.5, 0.6)};
    atoms[1].na = 1;
    atoms[1].mass = 1.0;
    atoms[1].ncpp.zv = 1.0;
    atoms[1].tau = {ModuleBase::Vector3<double>(0.7, 0.8, 0.9)};

    const std::vector<unitcell::AtomData> data = unitcell::get_atom_data(atoms, 2, 10.0);
    ASSERT_EQ(data.size(), 3u);
    EXPECT_DOUBLE_EQ(data[0].position.x, 1.0);
    EXPECT_DOUBLE_EQ(data[0].position.y, 2.0);
    EXPECT_DOUBLE_EQ(data[0].position.z, 3.0);
    EXPECT_DOUBLE_EQ(data[1].position.x, 4.0);
    EXPECT_DOUBLE_EQ(data[1].position.y, 5.0);
    EXPECT_DOUBLE_EQ(data[1].position.z, 6.0);
    EXPECT_DOUBLE_EQ(data[2].position.x, 7.0);
    EXPECT_DOUBLE_EQ(data[2].position.y, 8.0);
    EXPECT_DOUBLE_EQ(data[2].position.z, 9.0);
    EXPECT_DOUBLE_EQ(data[0].mass, 12.0);
    EXPECT_DOUBLE_EQ(data[1].mass, 12.0);
    EXPECT_DOUBLE_EQ(data[2].mass, 1.0);
    EXPECT_DOUBLE_EQ(data[0].valence_charge, 4.0);
    EXPECT_DOUBLE_EQ(data[1].valence_charge, 4.0);
    EXPECT_DOUBLE_EQ(data[2].valence_charge, 1.0);
    EXPECT_DOUBLE_EQ(atoms[0].tau[0].x, 0.1);
}

TEST(CellTools, AcceptsNoAtomTypes)
{
    const std::vector<unitcell::AtomData> data = unitcell::get_atom_data(nullptr, 0, 10.0);
    EXPECT_TRUE(data.empty());
}
