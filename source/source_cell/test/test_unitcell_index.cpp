#include "../unitcell.h"

#include <gtest/gtest.h>

class TestUnitCellIndex : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        // Build a minimal UnitCell:
        // 2 species, 2 atoms of the first species, 1 atom of the second
        ucell.ntype = 2;
        ucell.nat = 3;
        ucell.atoms = new Atom[2];
        ucell.atoms[0].na = 2;
        ucell.atoms[0].nw = 4; // 4 orbitals
        ucell.atoms[1].na = 1;
        ucell.atoms[1].nw = 2; // 2 orbitals
        ucell.set_atom_flag = true;
    }

    void TearDown() override
    {
        // Do NOT delete[] ucell.atoms: ~UnitCell() frees it (set_atom_flag = true).
        // Do NOT delete[] iat2it/iat2ia either: they are owned by the internal
        // Statistics member, whose destructor releases them (see AGENTS.md).
    }

    UnitCell ucell;
};

TEST_F(TestUnitCellIndex, SetIat2itia)
{
    ucell.set_iat2itia();
    // iat2it: [0, 0, 1]
    EXPECT_EQ(ucell.iat2it[0], 0);
    EXPECT_EQ(ucell.iat2it[1], 0);
    EXPECT_EQ(ucell.iat2it[2], 1);
    // iat2ia: [0, 1, 0]
    EXPECT_EQ(ucell.iat2ia[0], 0);
    EXPECT_EQ(ucell.iat2ia[1], 1);
    EXPECT_EQ(ucell.iat2ia[2], 0);
}

TEST_F(TestUnitCellIndex, SetIat2itiaCalledTwice)
{
    // Cover the delete[] + new[] reallocation path: the second call must
    // produce the same mapping as the first.
    ucell.set_iat2itia();
    ucell.set_iat2itia();
    EXPECT_EQ(ucell.iat2it[0], 0);
    EXPECT_EQ(ucell.iat2it[1], 0);
    EXPECT_EQ(ucell.iat2it[2], 1);
    EXPECT_EQ(ucell.iat2ia[0], 0);
    EXPECT_EQ(ucell.iat2ia[1], 1);
    EXPECT_EQ(ucell.iat2ia[2], 0);
}

TEST_F(TestUnitCellIndex, SetIat2iwtNpol1)
{
    ucell.set_iat2iwt(1);
    // iat=0: iwt=0, iat=1: iwt=4, iat=2: iwt=8
    EXPECT_EQ(ucell.get_iat2iwt()[0], 0);
    EXPECT_EQ(ucell.get_iat2iwt()[1], 4);
    EXPECT_EQ(ucell.get_iat2iwt()[2], 8);
}

TEST_F(TestUnitCellIndex, SetIat2iwtNpol2)
{
    ucell.set_iat2iwt(2);
    // iat=0: iwt=0, iat=1: iwt=8, iat=2: iwt=16
    EXPECT_EQ(ucell.get_iat2iwt()[0], 0);
    EXPECT_EQ(ucell.get_iat2iwt()[1], 8);
    EXPECT_EQ(ucell.get_iat2iwt()[2], 16);
}
