#include <gtest/gtest.h>
#include "module_base/ndarray.h"
#include <utility>

// TEST(NDArray, Constructor)
// {
//     const NDArray<int> a;
//     EXPECT_EQ(a.size(), 0);
//     EXPECT_EQ(a.empty(), true);
// }

TEST(NDArray, InitializerListConstructor)
{
    const NDArray<int> a({1, 2, 3}); /* 1 * 2 * 3, 3d array */
    EXPECT_EQ(a.size(), 6);
    EXPECT_EQ(a.empty(), false);
}

TEST(NDArray, VariadicTemplateConstructor)
{
    const NDArray<int> a(1, 2, 3); /* 1 * 2 * 3, 3d array */
    EXPECT_EQ(a.size(), 6);
    EXPECT_EQ(a.empty(), false);
}

TEST(NDArray, CopyConstructor)
{
    const NDArray<int> a(1, 2, 3); /* 1 * 2 * 3, 3d array */
    const NDArray<int> b(a);
    EXPECT_EQ(b.size(), 6);
    EXPECT_EQ(b.empty(), false);
    // and a will be the same as b
    EXPECT_EQ(a.size(), 6);
    EXPECT_EQ(a.empty(), false);
}

TEST(NDArray, MoveConstructor)
{
    NDArray<int> a(1, 2, 3); /* 1 * 2 * 3, 3d array */
    const NDArray<int> b(std::move(a));
    EXPECT_EQ(b.size(), 6);
    EXPECT_EQ(b.empty(), false);
    // and a will be empty, but still valid (principle of std::move)
    EXPECT_EQ(a.size(), 0);
    EXPECT_EQ(a.empty(), true);
}

TEST(NDArray, CopyAssignment)
{
    const NDArray<int> a(1, 2, 3); /* 1 * 2 * 3, 3d array */
    NDArray<int> b(1);
    b = a;
    EXPECT_EQ(b.size(), 6);
    EXPECT_EQ(b.empty(), false);
    // and a will be the same as b
    EXPECT_EQ(a.size(), 6);
    EXPECT_EQ(a.empty(), false);
}

TEST(NDArray, MoveAssignment)
{
    NDArray<int> a(1, 2, 3); /* 1 * 2 * 3, 3d array */
    NDArray<int> b(1);
    b = std::move(a);
    EXPECT_EQ(b.size(), 6);
    EXPECT_EQ(b.empty(), false);
    // and a will be empty, but still valid (principle of std::move)
    EXPECT_EQ(a.size(), 0);
    EXPECT_EQ(a.empty(), true);
}

TEST(NDArray, EqualityOperator)
{
    const NDArray<int> a(1, 2, 3); /* 1 * 2 * 3, 3d array */
    const NDArray<int> b(1, 2, 3); /* 1 * 2 * 3, 3d array */
    const NDArray<int> c(1, 2, 4); /* 1 * 2 * 4, 3d array */
    EXPECT_EQ(a == b, true);
    EXPECT_EQ(a == c, false);
}

TEST(NDArray, InequalityOperator)
{
    const NDArray<int> a(1, 2, 3); /* 1 * 2 * 3, 3d array */
    const NDArray<int> b(1, 2, 3); /* 1 * 2 * 3, 3d array */
    const NDArray<int> c(1, 2, 4); /* 1 * 2 * 4, 3d array */
    EXPECT_EQ(a != b, false);
    EXPECT_EQ(a != c, true);
}

TEST(NDArray, Index)
{
    const NDArray<int> a(1, 2, 3); /* 1 * 2 * 3, 3d array */
    EXPECT_EQ(a.index(0, 0, 0), 0);
    EXPECT_EQ(a.index(0, 0, 1), 1);
    EXPECT_EQ(a.index(0, 0, 2), 2);
    EXPECT_EQ(a.index(0, 1, 0), 3);
    EXPECT_EQ(a.index(0, 1, 1), 4);
    EXPECT_EQ(a.index(0, 1, 2), 5);
}

TEST(NDArray, AtMethodMultiIndex)
{
    NDArray<int> a(1, 2, 3); /* 1 * 2 * 3, 3d array */
    a.at(0, 0, 0) = 1;
    a.at(0, 0, 1) = 2;
    a.at(0, 0, 2) = 3;
    a.at(0, 1, 0) = 4;
    a.at(0, 1, 1) = 5;
    a.at(0, 1, 2) = 6;
    EXPECT_EQ(a.at(0, 0, 0), 1);
    EXPECT_EQ(a.at(0, 0, 1), 2);
    EXPECT_EQ(a.at(0, 0, 2), 3);
    EXPECT_EQ(a.at(0, 1, 0), 4);
    EXPECT_EQ(a.at(0, 1, 1), 5);
    EXPECT_EQ(a.at(0, 1, 2), 6);
}

TEST(NDArray, IndexOperatorMultiIndex)
{
    NDArray<int> a(1, 2, 3);
    a(0, 0, 0) = 1;
    a(0, 0, 1) = 2;
    a(0, 0, 2) = 3;
    a(0, 1, 0) = 4;
    a(0, 1, 1) = 5;
    a(0, 1, 2) = 6;
    EXPECT_EQ(a(0, 0, 0), 1);
    EXPECT_EQ(a(0, 0, 1), 2);
    EXPECT_EQ(a(0, 0, 2), 3);
    EXPECT_EQ(a(0, 1, 0), 4);
    EXPECT_EQ(a(0, 1, 1), 5);
    EXPECT_EQ(a(0, 1, 2), 6);
}

TEST(NDArray, Reshape)
{
    NDArray<int> a(1, 2, 3); /* 1 * 2 * 3, 3d array */
    a.reshape(2, 3, 1); /* 2 * 3 * 1, 3d array */
    EXPECT_EQ(a.size(), 6);
    EXPECT_EQ(a.empty(), false);
    // expect assert error if the size is not the same
    EXPECT_DEATH(a.reshape(2, 3, 2), "");
}

TEST(NDArray, ReshapeValue)
{
    NDArray<int> a(1, 2, 3); /* 1 * 2 * 3, 3d array */
    /*
    [
        [
            [1, 2, 3],
            [4, 5, 6]
        ]
    ] // in sequence of 1, 2, 3, 4, 5, 6
    */
    a(0, 0, 0) = 1;
    a(0, 0, 1) = 2;
    a(0, 0, 2) = 3;
    a(0, 1, 0) = 4;
    a(0, 1, 1) = 5;
    a(0, 1, 2) = 6;
    a.reshape(2, 3, 1); /* 2 * 3 * 1, 3d array */
    /*
    [
        [
            [1],
            [2],
            [3]
        ],
        [
            [4],
            [5],
            [6]
        ]
    ]
    */
    EXPECT_EQ(a(0, 0, 0), 1);
    EXPECT_EQ(a(0, 1, 0), 2);
    EXPECT_EQ(a(0, 2, 0), 3);
    EXPECT_EQ(a(1, 0, 0), 4);
    EXPECT_EQ(a(1, 1, 0), 5);
    EXPECT_EQ(a(1, 2, 0), 6);
}

TEST(NDArray, ReshapeInfer)
{
    const NDArray<int> a(1, 2, 3); /* 1 * 2 * 3, 3d array */
    // infer the first dimension
    // infer the second dimension
    // infer the last dimension
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
