#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../sltk_adjacent_set.h"

/************************************************
 *  unit test of sltk_adjacent_set
 ***********************************************/

/**
 * - Tested Functions:
 *   - SetWithoutExpand()
 *     - set box index without expand
 *   - SetWithExpand()
 *     - set box index with expand
 */

class sltkAdjacentSet : public testing::Test
{
protected:
    AdjacentSet test;
};

TEST_F(sltkAdjacentSet, SetWithoutExpand)
{
    // set parameters
    test.setDx(3);
    test.setDy(3);
    test.setDz(3);
    test.setCenter(13);
    test.setExpandFlag(false);
    test.setTrueX(1);
    test.setTrueY(1);
    test.setTrueZ(1);
    // set box and offset
    int x=1;
    int y=1;
    int z=1;
    int offset=1;
    int test_grid_in=1;
    test.set(x,y,z,offset,test_grid_in);
    EXPECT_EQ(test.box[0],26);
    // set box and offset
    x=0;
    y=1;
    z=-1;
    test.set(x,y,z,offset,test_grid_in);
    EXPECT_EQ(test.box[1],15);
    // check length of box and offset
    EXPECT_EQ(test.box.size(), 2);
    EXPECT_EQ(test.offset.size(), 2);
    EXPECT_EQ(test.getLength(), 2);
    // getBox
    int box_x = 0;
    int box_y = 0;
    int box_z = 0;
    test.getBox(26, box_x, box_y, box_z);
    EXPECT_EQ(box_x, 1);
    EXPECT_EQ(box_y, 1);
    EXPECT_EQ(box_z, 1);
    test.getBox(15, box_x, box_y, box_z);
    EXPECT_EQ(box_x, 0);
    EXPECT_EQ(box_y, 1);
    EXPECT_EQ(box_z, -1);
    // delete_vector
    test.delete_vector();
    EXPECT_EQ(test.box.size(),0);
    EXPECT_EQ(test.offset.size(), 0);
}

TEST_F(sltkAdjacentSet, SetWithExpand)
{
    // set parameters
    test.setDx(6);
    test.setDy(6);
    test.setDz(6);
    test.setCenter(43);
    test.setExpandFlag(true);
    test.setTrueX(1);
    test.setTrueY(1);
    test.setTrueZ(1);
    // set box and offset
    int x = 0;
    int y = 0;
    int z = 0;
    int offset = 1;
    int test_grid_in = 1;
    test.set(x, y, z, offset, test_grid_in);
    EXPECT_EQ(test.box[0], 43);
    // set box and offset
    x = 1;
    y = 0;
    z = 0;
    test.set(x, y, z, offset, test_grid_in);
    EXPECT_EQ(test.box[1], 79);
    // check length of box and offset
    EXPECT_EQ(test.box.size(), 2);
    EXPECT_EQ(test.offset.size(), 2);
    EXPECT_EQ(test.getLength(), 2);
    // getBox
    int box_x;
    int box_y;
    int box_z;
    test.getBox(43, box_x, box_y, box_z);
    EXPECT_EQ(box_x, 0);
    EXPECT_EQ(box_y, 0);
    EXPECT_EQ(box_z, 0);
    test.getBox(79, box_x, box_y, box_z);
    EXPECT_EQ(box_x, 1);
    EXPECT_EQ(box_y, 0);
    EXPECT_EQ(box_z, 0);
    // delete_vector
    test.delete_vector();
    EXPECT_EQ(test.box.size(), 0);
    EXPECT_EQ(test.offset.size(), 0);
}
