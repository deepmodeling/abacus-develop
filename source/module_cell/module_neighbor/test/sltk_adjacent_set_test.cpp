#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../sltk_adjacent_set.h"

/************************************************
 *  unit test of sltk_adjacent_set
 ***********************************************/

/**
 * - Tested Functions:
 *   - delete_vector()
 *     - delete vector offset and box
 *   - assertCoordinateIsRight()
 *     - confirm that the coordinate is right
 *   - Index()
 *     - Use ternary to represent coordinates     
 */

class sltkAdjacentSet : public testing::Test
{
protected:
    AdjacentSet test;
};

TEST_F(sltkAdjacentSet,deleteVectorAndIndex)
{
    int x=1;
    int y=1;
    int z=1;
    int offset=1;
    int test_grid_in=1;
    test.set(x,y,z,offset,test_grid_in);
    EXPECT_EQ(test.box[0],26);
    x=0;
    y=1;
    z=-1;
    test.set(x,y,z,offset,test_grid_in);
    EXPECT_EQ(test.box[1],15);
    test.delete_vector();
    EXPECT_EQ(test.box.size(),0);
}
/*TEST_F(sltkAdjacentSet,assertCoordinateIsRight)
{
	std::string output;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ModuleBase::assertCoordinateIsRight(), ::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("
    
    "));
}
*/
