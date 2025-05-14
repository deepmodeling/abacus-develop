#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../sltk_atom.h"

/************************************************
 *  unit test of sltk_atom
 ***********************************************/

/**
 * - Tested Functions:
 *   - AllocateAdjacentSet
 *     - FAtom::allocate_AdjacentSet()
 *       - FAtom has a member "as", which is a shared_ptr of AdjacentSet
 *       - allocate_AdjacentSet() is a function of FAtom to
 *         allocate "as"
 *   - SetterGetters:
 *     - the setter and getter of FAtom
 *     - including d_x, d_y, d_z, type, natom
 */

class SltkAtomTest : public testing::Test
{
protected:
    FAtom test;
};


TEST_F(SltkAtomTest, SetterGetters)
{
    FAtom test_temp(1.0, 2.0, 3.0, 4, 5, 0, 1, 2);

    EXPECT_DOUBLE_EQ(test_temp.x, 1.0);
    EXPECT_DOUBLE_EQ(test_temp.y, 2.0);
    EXPECT_DOUBLE_EQ(test_temp.z, 3.0);
    EXPECT_EQ(test_temp.type, 4);
    EXPECT_EQ(test_temp.natom, 5);
    EXPECT_EQ(test_temp.cell_x, 0);
    EXPECT_EQ(test_temp.cell_y, 1);
    EXPECT_EQ(test_temp.cell_z, 2);
}
