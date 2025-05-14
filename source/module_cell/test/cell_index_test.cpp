#include "module_cell/cell_index.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <fstream>

/************************************************
 *  unit test of class CellIndex
 ***********************************************/

/**
 * - Tested Functions:
 *   - Indices used in UnitCell
 **/

class CellIndexTest : public testing::Test
{
  protected:
    std::vector<std::string> atom_labels = {"C", "H"};
    std::vector<int> atom_counts = {1, 2};
    std::vector<std::vector<int>> lnchi_counts = {{1, 1, 1}, {1, 1, 1}};
    CellIndex cell_index = CellIndex(atom_labels, atom_counts, lnchi_counts, 1);
};

TEST_F(CellIndexTest, EmptyTest)
{
    CellIndex cell_index1;
    EXPECT_EQ(0, cell_index1.get_ntype());
    EXPECT_EQ(0, cell_index1.get_nw());
    testing::internal::CaptureStdout();
    EXPECT_EXIT(cell_index1.get_atom_label(0), ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("iat out of range [0, nat)"));
}

TEST_F(CellIndexTest, Index)
{
    EXPECT_EQ(2, cell_index.get_ntype());
    EXPECT_EQ(3, cell_index.get_nat());
    EXPECT_EQ(1, cell_index.get_nat(0));
    EXPECT_EQ(2, cell_index.get_nat(1));
    EXPECT_EQ(27, cell_index.get_nw());
    EXPECT_EQ(9, cell_index.get_nw(0));
    EXPECT_EQ(9, cell_index.get_nw(1));
    EXPECT_EQ(0, cell_index.get_iwt(0, 0));
    EXPECT_EQ(9, cell_index.get_iwt(1, 0));
    EXPECT_EQ(18, cell_index.get_iwt(2, 0));
    EXPECT_EQ(2, cell_index.get_maxL(0));
    EXPECT_EQ(2, cell_index.get_maxL(1));
    EXPECT_EQ(2, cell_index.get_maxL(2));
    EXPECT_EQ(1, cell_index.get_nchi(0, 0));
    EXPECT_EQ("C", cell_index.get_atom_label(0));
    EXPECT_EQ("H", cell_index.get_atom_label(1));
}

TEST_F(CellIndexTest, WriteOrbInfo)
{
    cell_index.write_orb_info("./");
    std::ifstream ifs("./Orbital");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("#io    spec    l    m    z  sym"));
    EXPECT_THAT(str, testing::HasSubstr("0       C    2    4    1            dxy"));
    EXPECT_THAT(str, testing::HasSubstr("1       H    2    4    1            dxy"));
    EXPECT_THAT(str, testing::HasSubstr("2       H    2    4    1            dxy"));
    remove("./Orbital");
}