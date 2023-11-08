#include <gtest/gtest.h>
#include <ATen/core/tensor_accessor.h> // Include the header file you provided

namespace container {

// Test fixture to set up common data for tests
class TensorAccessorTest : public testing::Test {
protected:
    // Common setup code
    TensorAccessorTest() = default;

    // Common cleanup code
    virtual ~TensorAccessorTest() = default;
};

// Test the TensorAccessor class
TEST_F(TensorAccessorTest, TensorAccessorTest) {
    // Test data
    int data[6] = {1, 2, 3, 4, 5, 6};
    int sizes[3] = {2, 3, 1};
    int strides[3] = {3, 1, 1};

    TensorAccessor<int, 3, int> accessor(data, sizes, strides);

    // Test operator[] for 1D TensorAccessor
    EXPECT_EQ(accessor[0][0][0], 1);
    EXPECT_EQ(accessor[0][1][0], 2);
    EXPECT_EQ(accessor[0][2][0], 3);
    EXPECT_EQ(accessor[1][0][0], 4);
    EXPECT_EQ(accessor[1][1][0], 5);
    EXPECT_EQ(accessor[1][2][0], 6);

    // Test operator[] for 2D TensorAccessor
    auto sub_accessor_1 = accessor[1];
    EXPECT_EQ(sub_accessor_1[0][0], 4);
    EXPECT_EQ(sub_accessor_1[1][0], 5);
    EXPECT_EQ(sub_accessor_1[2][0], 6);

    auto sub_accessor_2 = accessor[1][0];
    EXPECT_EQ(sub_accessor_2[0], 4);
}

} // namespace container
