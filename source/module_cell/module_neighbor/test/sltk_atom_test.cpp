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
 *     - FAtom::delete_vector()
 *       - this function will call AdjacentSet::delete_vector() to
 *       - delete the box and offset vectors in AdjacentSet
 *   - SetterGetters:
 *     - the setter and getter of FAtom
 *     - including d_x, d_y, d_z, type, natom
 */

class SltkAtomTest : public testing::Test
{
protected:
    FAtom test;
};

TEST_F(SltkAtomTest, AllocateAdjacentSet)
{
    test.clearAdjacent();
    FAtom test_temp(1.0, 2.0, 3.0, 4, 5, 0, 1, 2);
    test.addAdjacent(test_temp);
    EXPECT_EQ(test.getAdjacent().front()->getType(), 4);
}


TEST_F(SltkAtomTest, SetterGetters)
{
    FAtom test_temp(1.0, 2.0, 3.0, 4, 5, 0, 1, 2);

    EXPECT_DOUBLE_EQ(test_temp.x(), 1.0);
    EXPECT_DOUBLE_EQ(test_temp.y(), 2.0);
    EXPECT_DOUBLE_EQ(test_temp.z(), 3.0);
    EXPECT_EQ(test_temp.getType(), 4);
    EXPECT_EQ(test_temp.getNatom(), 5);
    EXPECT_EQ(test_temp.getCellX(), 0);
    EXPECT_EQ(test_temp.getCellY(), 1);
    EXPECT_EQ(test_temp.getCellZ(), 2);
}
