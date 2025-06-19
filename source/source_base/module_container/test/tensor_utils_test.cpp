#include <gtest/gtest.h>

#include <ATen/core/tensor_utils.h>


TEST(TensorUtils, _get_digit_places) {
    const int size = 6;
    float arr[size] = {-1.0, 2.5, 3.0, -0.5, 0.0, 1.234567};
    int int_count = 0, frac_count = 0;

    // Test for float type
    int total_digits = container::_get_digit_places(arr, size, int_count, frac_count);
    EXPECT_EQ(total_digits, 9);
    EXPECT_EQ(int_count, 2);
    EXPECT_EQ(frac_count, 7);
}

TEST(TensorUtils, _internal_output) {
    const int num_elements = 8;
    float* data = new float[num_elements];
    for (int ii = 0; ii < 8; ii++) {
        data[ii] = ii;
    }
    container::TensorShape shape1 {8}, shape2{2, 4}, shape3{2, 2, 2}, shape4{1, 2, 2, 2};
    // Test if the output operator produces the expected output
    std::ostringstream oss;
    container::_internal_output(oss, data, shape1, num_elements);
    container::_internal_output(oss, data, shape2, num_elements);
    container::_internal_output(oss, data, shape3, num_elements);
    container::_internal_output(oss, data, shape4, num_elements);
    const std::string expected_output = "[ 0, 1, 2, 3, 4, 5, 6, 7]";
    EXPECT_TRUE(oss.str().find(expected_output) == 0);
}

TEST(TensorUtils, removeTrailingZeros) {
    EXPECT_EQ(container::removeTrailingZeros(""), "0");
}