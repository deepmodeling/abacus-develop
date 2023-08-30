#include "../basic_funcs.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

/************************************************
 *  unit test of basic_funcs
 ***********************************************/

/**
 * - Tested Functions:
 *  - maxval_abs_2d()
 */

class BasicFuncsTest : public testing::Test
{
  protected:
    std::vector<std::vector<double>> array = {
        {1.0, 2.0, 3.0},
        {4.0, 5.0, 6.0},
        {7.0, 8.0, 9.0}
    };
    std::string output;
};

TEST_F(BasicFuncsTest, MaxvalAbs2d)
{
    EXPECT_DOUBLE_EQ(maxval_abs_2d(array), 9.0);
}

TEST_F(BasicFuncsTest, MaxlocAbs2d)
{
    std::vector<int> result(2);
    maxloc_abs_2d(array, result);
    EXPECT_EQ(result[0], 2);
    EXPECT_EQ(result[1], 2);
}

TEST_F(BasicFuncsTest, Sum2d)
{
    std::vector<ModuleBase::Vector3<double>> array_1;
    array_1.push_back(ModuleBase::Vector3<double>(1.0, 2.0, 3.0));
    array_1.push_back(ModuleBase::Vector3<double>(4.0, 5.0, 6.0));
    array_1.push_back(ModuleBase::Vector3<double>(7.0, 8.0, 9.0));
    EXPECT_DOUBLE_EQ(sum_2d(array_1), 45.0);
}

TEST_F(BasicFuncsTest, ScalarMul2d)
{
    std::vector<std::vector<double>> result;
    scalar_multiply_2d(array, 2.0, result);
    EXPECT_DOUBLE_EQ(result[0][0], 2.0);
    EXPECT_DOUBLE_EQ(result[0][1], 4.0);
    EXPECT_DOUBLE_EQ(result[0][2], 6.0);
    EXPECT_DOUBLE_EQ(result[1][0], 8.0);
    EXPECT_DOUBLE_EQ(result[1][1], 10.0);
    EXPECT_DOUBLE_EQ(result[1][2], 12.0);
    EXPECT_DOUBLE_EQ(result[2][0], 14.0);
    EXPECT_DOUBLE_EQ(result[2][1], 16.0);
    EXPECT_DOUBLE_EQ(result[2][2], 18.0);
}

TEST_F(BasicFuncsTest, AddScalarMul2d)
{
    std::vector<ModuleBase::Vector3<double>> array_1, array_2, result;
    array_1.push_back(ModuleBase::Vector3<double>(1.0, 2.0, 3.0));
    array_1.push_back(ModuleBase::Vector3<double>(4.0, 5.0, 6.0));
    array_1.push_back(ModuleBase::Vector3<double>(7.0, 8.0, 9.0));
    array_2.push_back(ModuleBase::Vector3<double>(1.0, 2.0, 3.0));
    array_2.push_back(ModuleBase::Vector3<double>(4.0, 5.0, 6.0));
    array_2.push_back(ModuleBase::Vector3<double>(7.0, 8.0, 9.0));
    add_scalar_multiply_2d(array_1, array_2, 2.0, result);
    EXPECT_DOUBLE_EQ(result[0][0], 3.0);
    EXPECT_DOUBLE_EQ(result[0][1], 6.0);
    EXPECT_DOUBLE_EQ(result[0][2], 9.0);
    EXPECT_DOUBLE_EQ(result[1][0], 12.0);
    EXPECT_DOUBLE_EQ(result[1][1], 15.0);
    EXPECT_DOUBLE_EQ(result[1][2], 18.0);
    EXPECT_DOUBLE_EQ(result[2][0], 21.0);
    EXPECT_DOUBLE_EQ(result[2][1], 24.0);
    EXPECT_DOUBLE_EQ(result[2][2], 27.0);
}

TEST_F(BasicFuncsTest, Subtract2d)
{
    std::vector<std::vector<double>> result;
    std::vector<std::vector<double>> array_2 = {
        {1.0, 2.0, 3.0},
        {4.0, 5.0, 6.0},
        {7.0, 8.0, 9.0}
    };
    subtract_2d(array, array_2, result);
    EXPECT_DOUBLE_EQ(result[0][0], 0.0);
    EXPECT_DOUBLE_EQ(result[0][1], 0.0);
    EXPECT_DOUBLE_EQ(result[0][2], 0.0);
    EXPECT_DOUBLE_EQ(result[1][0], 0.0);
    EXPECT_DOUBLE_EQ(result[1][1], 0.0);
    EXPECT_DOUBLE_EQ(result[1][2], 0.0);
    EXPECT_DOUBLE_EQ(result[2][0], 0.0);
    EXPECT_DOUBLE_EQ(result[2][1], 0.0);
    EXPECT_DOUBLE_EQ(result[2][2], 0.0);
}

TEST_F(BasicFuncsTest, FillScalar2d)
{
    std::vector<std::vector<double>> result;
    result.resize(3);
    for (int i = 0; i < 3; i++)
    {
        result[i].resize(3);
    }
    fill_scalar_2d(2.0, result);
    EXPECT_DOUBLE_EQ(result[0][0], 2.0);
    EXPECT_DOUBLE_EQ(result[0][1], 2.0);
    EXPECT_DOUBLE_EQ(result[0][2], 2.0);
    EXPECT_DOUBLE_EQ(result[1][0], 2.0);
    EXPECT_DOUBLE_EQ(result[1][1], 2.0);
    EXPECT_DOUBLE_EQ(result[1][2], 2.0);
    EXPECT_DOUBLE_EQ(result[2][0], 2.0);
    EXPECT_DOUBLE_EQ(result[2][1], 2.0);
    EXPECT_DOUBLE_EQ(result[2][2], 2.0);
}

TEST_F(BasicFuncsTest, WhereFillScalar2d)
{
    std::vector<std::vector<int>> array_mask = {
        {1, 0, 1},
        {0, 1, 0},
        {1, 0, 1}
    };
    std::vector<std::vector<double>> result;
    where_fill_scalar_2d(array_mask, 1, 2.0, result);
    EXPECT_DOUBLE_EQ(result[0][0], 2.0);
    EXPECT_DOUBLE_EQ(result[0][1], 0.0);
    EXPECT_DOUBLE_EQ(result[0][2], 2.0);
    EXPECT_DOUBLE_EQ(result[1][0], 0.0);
    EXPECT_DOUBLE_EQ(result[1][1], 2.0);
    EXPECT_DOUBLE_EQ(result[1][2], 0.0);
    EXPECT_DOUBLE_EQ(result[2][0], 2.0);
    EXPECT_DOUBLE_EQ(result[2][1], 0.0);
    EXPECT_DOUBLE_EQ(result[2][2], 2.0);
}

TEST_F(BasicFuncsTest, WhereFillScalarElse2d)
{
    std::vector<ModuleBase::Vector3<int>> array_mask;
    array_mask.push_back(ModuleBase::Vector3<int>(1,0,1));
    array_mask.push_back(ModuleBase::Vector3<int>(0,1,0));
    array_mask.push_back(ModuleBase::Vector3<int>(1,0,1));
    std::vector<ModuleBase::Vector3<double>> result;
    std::vector<ModuleBase::Vector3<double>> rest;
    rest.push_back(ModuleBase::Vector3<double>(1.0, 2.0, 3.0));
    rest.push_back(ModuleBase::Vector3<double>(4.0, 5.0, 6.0));
    rest.push_back(ModuleBase::Vector3<double>(7.0, 8.0, 9.0));
    where_fill_scalar_else_2d(array_mask, 1, 2.0, rest, result);
    EXPECT_DOUBLE_EQ(result[0][0], 2.0);
    EXPECT_DOUBLE_EQ(result[0][1], 2.0);
    EXPECT_DOUBLE_EQ(result[0][2], 2.0);
    EXPECT_DOUBLE_EQ(result[1][0], 4.0);
    EXPECT_DOUBLE_EQ(result[1][1], 2.0);
    EXPECT_DOUBLE_EQ(result[1][2], 6.0);
    EXPECT_DOUBLE_EQ(result[2][0], 2.0);
    EXPECT_DOUBLE_EQ(result[2][1], 8.0);
    EXPECT_DOUBLE_EQ(result[2][2], 2.0);
}

TEST_F(BasicFuncsTest, Print2D)
{
    std::vector<std::vector<double>> array = {
        {1.0, 2.0, 3.0},
        {4.0, 5.0, 6.0},
        {7.0, 8.0, 9.0}
    };
    testing::internal::CaptureStdout();
    print_2d(array);
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("1 2 3"));
    EXPECT_THAT(output, testing::HasSubstr("4 5 6"));
    EXPECT_THAT(output, testing::HasSubstr("7 8 9"));
}