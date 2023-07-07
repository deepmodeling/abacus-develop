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

class sltkatom : public testing::Test
{
protected:
    FAtom test;
};

TEST_F(sltkatom, AllocateAdjacentSet)
{
    test.allocate_AdjacentSet();
    std::shared_ptr<AdjacentSet> a = test.getAdjacentSet();
    int x = 1;
    int y = 1;
    int z = 1;
    int offset = 1;
    int test_grid_in = 1;
    a->set(x, y, z, offset, test_grid_in);
    EXPECT_EQ(a->box[0], 26);

    test.delete_vector();
    std::shared_ptr<AdjacentSet> b = test.getAdjacentSet();
    EXPECT_EQ(b->box.size(), 0);
}

TEST_F(sltkatom, SetterGetters)
{
    test.setX(1.0);
    test.setY(2.0);
    test.setZ(3.0);
    test.setType(4);
    test.setNatom(5);
    EXPECT_DOUBLE_EQ(test.x(), 1.0);
    EXPECT_DOUBLE_EQ(test.y(), 2.0);
    EXPECT_DOUBLE_EQ(test.z(), 3.0);
    EXPECT_EQ(test.getType(), 4);
    EXPECT_EQ(test.getNatom(), 5);
}
