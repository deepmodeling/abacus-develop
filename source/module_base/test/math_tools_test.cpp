#include "../math_tools.h"
#include "gtest/gtest.h"

/************************************************
 *  unit test of math tools
 ***********************************************/

TEST(Tools,LinearFit)
{
    std::vector<double> x;
    std::vector<double> y;
    x = {1, 2, 3, 4, 5, 6, 7, 8};
    y = {2, 4, 6, 8, 10, 12, 14, 16}; // y = 2x
    std::pair<double, double> result = ModuleBase::linear_fit(x, y);
    EXPECT_DOUBLE_EQ(result.first, 2.0); // slope = 2.0
    EXPECT_DOUBLE_EQ(result.second, 0.0); // b = 0.0

    y = {1.8, 3.9, 6.1, 8.2, 9.8, 12.3, 13.9, 16.1};
    result = ModuleBase::linear_fit(x, y);
    EXPECT_NEAR(result.first, 2.0, 0.1); // slope ~ 2.0
    EXPECT_NEAR(result.second, 0.0, 0.12); // b ~ 0.0
}
