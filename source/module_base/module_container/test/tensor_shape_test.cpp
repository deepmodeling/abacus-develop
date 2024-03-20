#include <gtest/gtest.h>

#include <ATen/core/tensor_shape.h>


/**
 * @brief Test cases for constructors of container::TensorShape class.
 */
TEST(TensorShape, Constructor) {
    // Test default constructor
    container::TensorShape shape1;
    EXPECT_EQ(shape1.ndim(), 0);

    // Test initializer_list constructor
    container::TensorShape shape2({2, 3, 4});
    EXPECT_EQ(shape2.ndim(), 3);
    EXPECT_EQ(shape2.dim_size(0), 2);
    EXPECT_EQ(shape2.dim_size(1), 3);
    EXPECT_EQ(shape2.dim_size(2), 4);

    // Test vector constructor
    std::vector<int64_t> dims = {5, 6};
    container::TensorShape shape3(dims);
    EXPECT_EQ(shape3.ndim(), 2);
    EXPECT_EQ(shape3.dim_size(0), 5);
    EXPECT_EQ(shape3.dim_size(1), 6);
}

/**
 * @brief Test cases for size manipulation functions of container::TensorShape class.
 */
TEST(TensorShape, SizeManipulation) {
    // Test add_dim and dim_size
    container::TensorShape shape({2, 3});
    shape.add_dim(4);
    EXPECT_EQ(shape.ndim(), 3);
    EXPECT_EQ(shape.dim_size(0), 2);
    EXPECT_EQ(shape.dim_size(1), 3);
    EXPECT_EQ(shape.dim_size(2), 4);

    // Test remove_dim
    shape.remove_dim(1);
    EXPECT_EQ(shape.ndim(), 2);
    EXPECT_EQ(shape.dim_size(0), 2);
    EXPECT_EQ(shape.dim_size(1), 4);

    // Test set_dim_size
    shape.set_dim_size(1, 5);
    EXPECT_EQ(shape.dim_size(1), 5);

    // Test NumElements
    EXPECT_EQ(shape.NumElements(), 10);
}

/**
 * @brief Test cases for comparison operators of container::TensorShape class.
 */
TEST(TensorShape, Comparison) {
    container::TensorShape shape1({2, 3, 4});
    container::TensorShape shape2({2, 3, 4});
    container::TensorShape shape3({3, 3, 4});

    // Test == operator
    EXPECT_EQ(shape1, shape2);
    EXPECT_NE(shape1, shape3);

    // Test != operator
    EXPECT_NE(shape1, shape3);
    EXPECT_NE(shape2, shape3);
}

/**
 * @brief Test cases for output stream operator of container::TensorShape class.
 */
TEST(TensorShape, Output) {
    container::TensorShape shape({2, 3, 4});
    std::stringstream ss;
    ss << shape;
    EXPECT_EQ(ss.str(), "[2,3,4]");
}
