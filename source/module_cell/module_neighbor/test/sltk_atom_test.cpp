#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../sltk_atom.h"

/************************************************
 *  unit test of sltk_atom
 ***********************************************/

/**
 * - Tested Functions:
 *   - FAtom::delete_vector(void)
 *     - delete the vector 
 */

class sltkatom : public testing::Test
{
};

TEST_F(sltkatom,deleteVector)
{
    FAtom test;
    test.allocate_AdjacentSet();
    std::shared_ptr<AdjacentSet> a=test.getAdjacentSet();
    int x=1;
    int y=1;
    int z=1;
    int offset=1;
    int test_grid_in=1;
    a->set(x,y,z,offset,test_grid_in);
    EXPECT_EQ(a->box[0],26);
    
    test.delete_vector();
    std::shared_ptr<AdjacentSet> b=test.getAdjacentSet();
    EXPECT_EQ(b->box.size(),0);
}

