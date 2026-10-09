#include "../matrix_methods.h"

#include <gtest/gtest.h>

#include <vector>

// Unit tests for the free matrix/vector helpers in matrix_methods.cpp.
// All functions are pure, so the tests assert exact values with hand-picked
// inputs. Two behaviors worth pinning down: DotInMAndV1 vs DotInMAndV2 differ
// in row/column orientation, and MPlus actually divides by its scalar.

TEST(MatrixMethods, ReshapeMToVFlattens)
{
    std::vector<ModuleBase::Vector3<double>> matrix;
    matrix.push_back(ModuleBase::Vector3<double>(1.0, 2.0, 3.0));
    matrix.push_back(ModuleBase::Vector3<double>(4.0, 5.0, 6.0));

    std::vector<double> result = ReshapeMToV(matrix);

    ASSERT_EQ(result.size(), 6);
    EXPECT_DOUBLE_EQ(result[0], 1.0);
    EXPECT_DOUBLE_EQ(result[1], 2.0);
    EXPECT_DOUBLE_EQ(result[2], 3.0);
    EXPECT_DOUBLE_EQ(result[3], 4.0);
    EXPECT_DOUBLE_EQ(result[4], 5.0);
    EXPECT_DOUBLE_EQ(result[5], 6.0);
}

TEST(MatrixMethods, ReshapeVToMGroupsByThree)
{
    std::vector<double> vec = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    std::vector<ModuleBase::Vector3<double>> result = ReshapeVToM(vec);

    ASSERT_EQ(result.size(), 2);
    EXPECT_DOUBLE_EQ(result[0].x, 1.0);
    EXPECT_DOUBLE_EQ(result[0].y, 2.0);
    EXPECT_DOUBLE_EQ(result[0].z, 3.0);
    EXPECT_DOUBLE_EQ(result[1].x, 4.0);
    EXPECT_DOUBLE_EQ(result[1].y, 5.0);
    EXPECT_DOUBLE_EQ(result[1].z, 6.0);
}

TEST(MatrixMethods, ReshapeRoundTrip)
{
    std::vector<ModuleBase::Vector3<double>> matrix;
    matrix.push_back(ModuleBase::Vector3<double>(1.5, -2.5, 3.5));
    matrix.push_back(ModuleBase::Vector3<double>(-4.5, 5.5, -6.5));

    std::vector<double> flat = ReshapeMToV(matrix);
    std::vector<ModuleBase::Vector3<double>> result = ReshapeVToM(flat);

    ASSERT_EQ(result.size(), matrix.size());
    for (size_t i = 0; i < matrix.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(result[i].x, matrix[i].x);
        EXPECT_DOUBLE_EQ(result[i].y, matrix[i].y);
        EXPECT_DOUBLE_EQ(result[i].z, matrix[i].z);
    }
}

TEST(MatrixMethods, MAddMElementwise)
{
    std::vector<std::vector<double>> a = {{1.0, 2.0}, {3.0, 4.0}};
    std::vector<std::vector<double>> b = {{5.0, 6.0}, {7.0, 8.0}};

    std::vector<std::vector<double>> result = MAddM(a, b);

    EXPECT_DOUBLE_EQ(result[0][0], 6.0);
    EXPECT_DOUBLE_EQ(result[0][1], 8.0);
    EXPECT_DOUBLE_EQ(result[1][0], 10.0);
    EXPECT_DOUBLE_EQ(result[1][1], 12.0);
}

TEST(MatrixMethods, MSubMElementwise)
{
    std::vector<std::vector<double>> a = {{5.0, 7.0}, {9.0, 11.0}};
    std::vector<std::vector<double>> b = {{1.0, 2.0}, {3.0, 4.0}};

    std::vector<std::vector<double>> result = MSubM(a, b);

    EXPECT_DOUBLE_EQ(result[0][0], 4.0);
    EXPECT_DOUBLE_EQ(result[0][1], 5.0);
    EXPECT_DOUBLE_EQ(result[1][0], 6.0);
    EXPECT_DOUBLE_EQ(result[1][1], 7.0);
}

TEST(MatrixMethods, VAddVElementwise)
{
    std::vector<double> a = {1.0, 2.0, 3.0};
    std::vector<double> b = {4.0, 5.0, 6.0};

    std::vector<double> result = VAddV(a, b);

    ASSERT_EQ(result.size(), 3);
    EXPECT_DOUBLE_EQ(result[0], 5.0);
    EXPECT_DOUBLE_EQ(result[1], 7.0);
    EXPECT_DOUBLE_EQ(result[2], 9.0);
}

TEST(MatrixMethods, VSubVElementwise)
{
    std::vector<double> a = {5.0, 7.0, 9.0};
    std::vector<double> b = {1.0, 2.0, 3.0};

    std::vector<double> result = VSubV(a, b);

    ASSERT_EQ(result.size(), 3);
    EXPECT_DOUBLE_EQ(result[0], 4.0);
    EXPECT_DOUBLE_EQ(result[1], 5.0);
    EXPECT_DOUBLE_EQ(result[2], 6.0);
}

TEST(MatrixMethods, DotInMAndV1RowMajor)
{
    // result[i] = sum_j matrix[i][j] * vec[j]
    std::vector<std::vector<double>> matrix = {{1.0, 2.0}, {3.0, 4.0}};
    std::vector<double> vec = {5.0, 6.0};

    std::vector<double> result = DotInMAndV1(matrix, vec);

    ASSERT_EQ(result.size(), 2);
    EXPECT_DOUBLE_EQ(result[0], 17.0); // 1*5 + 2*6
    EXPECT_DOUBLE_EQ(result[1], 39.0); // 3*5 + 4*6
}

TEST(MatrixMethods, DotInMAndV2ColumnMajor)
{
    // result[i] = sum_j matrix[j][i] * vec[j]
    std::vector<std::vector<double>> matrix = {{1.0, 2.0}, {3.0, 4.0}};
    std::vector<double> vec = {5.0, 6.0};

    std::vector<double> result = DotInMAndV2(matrix, vec);

    ASSERT_EQ(result.size(), 2);
    EXPECT_DOUBLE_EQ(result[0], 23.0); // 1*5 + 3*6
    EXPECT_DOUBLE_EQ(result[1], 34.0); // 2*5 + 4*6
}

TEST(MatrixMethods, DotInVAndVScalar)
{
    std::vector<double> vec1 = {1.0, 2.0, 3.0};
    std::vector<double> vec2 = {4.0, 5.0, 6.0};

    double result = DotInVAndV(vec1, vec2);

    EXPECT_DOUBLE_EQ(result, 32.0); // 4 + 10 + 18
}

TEST(MatrixMethods, OuterVAndVProducesMatrix)
{
    std::vector<double> a = {1.0, 2.0};
    std::vector<double> b = {3.0, 4.0};

    std::vector<std::vector<double>> result = OuterVAndV(a, b);

    EXPECT_DOUBLE_EQ(result[0][0], 3.0);
    EXPECT_DOUBLE_EQ(result[0][1], 4.0);
    EXPECT_DOUBLE_EQ(result[1][0], 6.0);
    EXPECT_DOUBLE_EQ(result[1][1], 8.0);
}

TEST(MatrixMethods, MPlusDividesByScalar)
{
    // Despite the name, MPlus divides every element by the scalar.
    std::vector<std::vector<double>> a = {{2.0, 4.0}, {6.0, 8.0}};

    std::vector<std::vector<double>> result = MPlus(a, 2.0);

    EXPECT_DOUBLE_EQ(result[0][0], 1.0);
    EXPECT_DOUBLE_EQ(result[0][1], 2.0);
    EXPECT_DOUBLE_EQ(result[1][0], 3.0);
    EXPECT_DOUBLE_EQ(result[1][1], 4.0);
}

TEST(MatrixMethods, DotInVAndFloatScales)
{
    std::vector<double> vec = {1.0, 2.0, 3.0};

    std::vector<double> result = DotInVAndFloat(vec, 3.0);

    ASSERT_EQ(result.size(), 3);
    EXPECT_DOUBLE_EQ(result[0], 3.0);
    EXPECT_DOUBLE_EQ(result[1], 6.0);
    EXPECT_DOUBLE_EQ(result[2], 9.0);
}
