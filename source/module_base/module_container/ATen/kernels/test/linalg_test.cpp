#include <gtest/gtest.h>

#include <ATen/core/tensor.h>
#include <ATen/kernels/linalg.h>
#include <test/test_utils.h>

namespace container {
namespace kernels {

template<typename T>
class LinalgTest : public testing::Test {
public:
    LinalgTest() = default;

    ~LinalgTest() override = default;
};

TYPED_TEST_SUITE(LinalgTest, test_utils::Types);

TYPED_TEST(LinalgTest, Add) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    Tensor A = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());
    Tensor B = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());

    Tensor expected = std::move(
        Tensor({static_cast<Type>(3.0), static_cast<Type>(6.0), static_cast<Type>(9.0),
                static_cast<Type>(12.0), static_cast<Type>(15.0), static_cast<Type>(18.0),
                static_cast<Type>(21.0), static_cast<Type>(24.0), static_cast<Type>(27.0)}).to_device<Device>());
    Tensor result = Tensor(expected.data_type(), expected.device_type(), expected.shape());
    kernels::add<Type, Device>()(
        A.NumElements(), static_cast<Type>(2.0), A.data<Type>(), static_cast<Type>(1.0), B.data<Type>(), result.data<Type>());
    EXPECT_EQ(result, expected);
}

TYPED_TEST(LinalgTest, Mul) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    Tensor A = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());
    Tensor B = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());

    Tensor expected = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(4.0), static_cast<Type>(9.0),
                static_cast<Type>(16.0), static_cast<Type>(25.0), static_cast<Type>(36.0),
                static_cast<Type>(49.0), static_cast<Type>(64.0), static_cast<Type>(81.0)}).to_device<Device>());
    Tensor result = Tensor(expected.data_type(), expected.device_type(), expected.shape());
    kernels::mul<Type, Device>()(
        A.NumElements(), static_cast<Type>(1.0), A.data<Type>(), B.data<Type>(), result.data<Type>());
    EXPECT_EQ(result, expected);
}

TYPED_TEST(LinalgTest, Div) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    Tensor A = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());
    Tensor B = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());

    Tensor expected = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(1.0), static_cast<Type>(1.0),
                static_cast<Type>(1.0), static_cast<Type>(1.0), static_cast<Type>(1.0),
                static_cast<Type>(1.0), static_cast<Type>(1.0), static_cast<Type>(1.0)}).to_device<Device>());

    Tensor result = Tensor(expected.data_type(), expected.device_type(), expected.shape());
    kernels::div<Type, Device>()(
        A.NumElements(), static_cast<Type>(1.0), A.data<Type>(), B.data<Type>(), result.data<Type>());
    EXPECT_EQ(result, expected);
}

TYPED_TEST(LinalgTest, Fma) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    Tensor A = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());
    Tensor B = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());
    Tensor C = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());

    Tensor expected = std::move(
        Tensor({static_cast<Type>(5.0),   static_cast<Type>(14.0),  static_cast<Type>(27.0),
                static_cast<Type>(44.0),  static_cast<Type>(65.0),  static_cast<Type>(90.0),
                static_cast<Type>(119.0), static_cast<Type>(152.0), static_cast<Type>(189.0)}).to_device<Device>());

    Tensor result = Tensor(expected.data_type(), expected.device_type(), expected.shape());
    kernels::fma<Type, Device>()(
        A.NumElements(), static_cast<Type>(2.0), A.data<Type>(), B.data<Type>(), static_cast<Type>(3.0), C.data<Type>(), result.data<Type>());
    EXPECT_EQ(result, expected);
}

TYPED_TEST(LinalgTest, Transpose) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    Tensor A = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());
    A.reshape({3, 3});
    Tensor A_transpose = A;
    Tensor expected = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(4.0), static_cast<Type>(7.0),
                static_cast<Type>(2.0), static_cast<Type>(5.0), static_cast<Type>(8.0),
                static_cast<Type>(3.0), static_cast<Type>(6.0), static_cast<Type>(9.0)}).to_device<Device>());
    expected.reshape({3, 3});
    std::vector<int> perm = {1, 0};

    kernels::transpose<Type, Device>()(
        perm, A.shape().dims(), A_transpose.shape().dims(), A.data<Type>(), A_transpose.data<Type>());
    EXPECT_EQ(A_transpose, expected);
}


TYPED_TEST(LinalgTest, Stride) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    Tensor A = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());
    A.reshape({-1});
    Tensor expected = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(5.0), static_cast<Type>(9.0)}).to_device<Device>());
    expected.reshape({-1});
    Tensor A_stride = expected;
    A_stride.zero();
    std::vector<int64_t> stride = {4};

    kernels::stride<Type, Device>()(
        stride, A.shape().dims(), A_stride.shape().dims(), A.data<Type>(), A_stride.data<Type>());
    EXPECT_EQ(A_stride, expected);
}


TYPED_TEST(LinalgTest, Inflate) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    Tensor expected = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                static_cast<Type>(0.0), static_cast<Type>(5.0), static_cast<Type>(0.0),
                static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(9.0)}).to_device<Device>());
    expected.reshape({-1});
    Tensor A = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(5.0), static_cast<Type>(9.0)}).to_device<Device>());
    A.reshape({-1});
    Tensor A_inflate = expected;
    A_inflate.zero();
    std::vector<int64_t> inflate = {4};

    kernels::inflate<Type, Device>()(
        inflate, A.shape().dims(), A_inflate.shape().dims(), A.data<Type>(), A_inflate.data<Type>());
    EXPECT_EQ(A_inflate, expected);
}


TYPED_TEST(LinalgTest, Reduce) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    Tensor A = std::move(
        Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());
    A.reshape({3, 3});
    Tensor expected = std::move(
        Tensor({static_cast<Type>(6.0), static_cast<Type>(15.0), static_cast<Type>(24.0)}).to_device<Device>());
    expected.reshape({-1});
    Tensor A_reduce = expected;
    A_reduce.zero();
    int64_t inner_most_dim = 3;

    kernels::reduce<Type, Device>()(
        A_reduce.NumElements(), inner_most_dim, A.data<Type>(), A_reduce.data<Type>());
    EXPECT_EQ(A_reduce, expected);
}

} // namespace kernels
} // namespace container