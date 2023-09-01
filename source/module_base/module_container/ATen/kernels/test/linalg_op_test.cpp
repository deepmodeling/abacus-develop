#include <gtest/gtest.h>

#include <ATen/core/tensor.h>
#include <ATen/kernels/linalg_op.h>
#include <ATen/kernels/test/op_test_utils.h>

namespace container {
namespace op {

template <typename T>
class LinalgOpTest : public testing::Test {
public:
    LinalgOpTest() {
        test_utils::init_blas_handle();
        test_utils::init_cusolver_handle();
    }
    ~LinalgOpTest() override {
        test_utils::delete_blas_handle();
        test_utils::delete_cusolver_handle();
    }
};

TYPED_TEST_SUITE(LinalgOpTest, test_utils::Types);

TYPED_TEST(LinalgOpTest, Transpose) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::transpose_op<Type, Device> transposeCalculator;

    const int dim = 3;
    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    A.reshape({-1, 3});
    Tensor B = A;
    B.zero();

    Tensor T = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(2.0), static_cast<Type>(4.0), static_cast<Type>(0.0),
                                 static_cast<Type>(3.0), static_cast<Type>(5.0), static_cast<Type>(6.0)}).to_device<Device>());
    T.reshape({-1, 3});

    transposeCalculator(A, {1, 0}, B);

    EXPECT_EQ(B, T);
}

TYPED_TEST(LinalgOpTest, Stride) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::stride_op<Type, Device> strideCalculator;

    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());

    Tensor T = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(4.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor B = T;
    B.zero();

    strideCalculator(A, {4}, B);

    EXPECT_EQ(B, T);
}

TYPED_TEST(LinalgOpTest, Inflate) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::inflate_op<Type, Device> inflateCalculator;

    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor B = A;
    B.zero();
    Tensor T = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(4.0), static_cast<Type>(6.0)}).to_device<Device>());
    

    inflateCalculator(T, {4}, B);

    EXPECT_EQ(A, B);
}

} // namespace op
} // namespace container
