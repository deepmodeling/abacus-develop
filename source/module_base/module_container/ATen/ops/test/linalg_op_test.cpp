#include <gtest/gtest.h>

#include <ATen/core/tensor.h>
#include <ATen/ops/linalg_op.h>
#include <test/test_utils.h>

namespace container {
namespace op {

template <typename T>
class LinalgOpTest : public testing::Test {
public:
    LinalgOpTest() = default;
    ~LinalgOpTest() override = default;
};

TYPED_TEST_SUITE(LinalgOpTest, test_utils::Types);

TYPED_TEST(LinalgOpTest, Add) {
    using Type   = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::add_op addCalculator;

    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor B = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor C = A;
    C.zero();

    Tensor expected = std::move(Tensor(
                                {static_cast<Type>(2.0), static_cast<Type>(4.0), static_cast<Type>(6.0),
                                 static_cast<Type>(0.0), static_cast<Type>(8.0), static_cast<Type>(10.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(12.0)}).to_device<Device>());
    addCalculator(A, B, C);
    EXPECT_EQ(C, expected);

    C.zero();
    C = A + B;
    EXPECT_EQ(C, expected);

    A += B;
    EXPECT_EQ(A, expected);

    C.zero();
    C = expected - B;
    C -= B;
    expected.zero();
    EXPECT_EQ(C, expected);
}

TYPED_TEST(LinalgOpTest, Sub) {
    using Type   = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::add_op addCalculator;

    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor B = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor C = A;
    C.zero();

    Tensor expected = std::move(Tensor(
        {static_cast<Type>(2.0), static_cast<Type>(4.0), static_cast<Type>(6.0),
         static_cast<Type>(0.0), static_cast<Type>(8.0), static_cast<Type>(10.0),
         static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(12.0)}).to_device<Device>());

    C = expected - B;
    C -= B;
    expected.zero();
    EXPECT_EQ(C, expected);
}

TYPED_TEST(LinalgOpTest, AddScalar) {
    using Type   = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::add_op addCalculator;

    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor B = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor C = A;
    C.zero();

    Tensor expected = std::move(Tensor({static_cast<Type>(-1.0), static_cast<Type>(-2.0), static_cast<Type>(-3.0),
                                        static_cast<Type>(0.0),  static_cast<Type>(-4.0), static_cast<Type>(-5.0),
                                        static_cast<Type>(0.0),  static_cast<Type>(0.0),  static_cast<Type>(-6.0)}).to_device<Device>());

    auto alpha = static_cast<Type>(2.0);
    auto beta  = static_cast<Type>(-3.0);

    addCalculator(alpha, A, beta, B, C);
    EXPECT_EQ(C, expected);
}

TYPED_TEST(LinalgOpTest, Mul) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::mul_op mulCalculator;

    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());

    Tensor B = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor C = A;
    C.zero();
    Tensor expected = std::move(Tensor({
                                 static_cast<Type>(1.0), static_cast<Type>(4.0), static_cast<Type>(9.0),
                                 static_cast<Type>(0.0), static_cast<Type>(16.0),static_cast<Type>(25.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(36.0)}).to_device<Device>());

    mulCalculator(A, B, C);
    EXPECT_EQ(C, expected);

    C.zero();
    C = A * B;
    EXPECT_EQ(C, expected);

    A *= B;
    EXPECT_EQ(A, expected);
}

TYPED_TEST(LinalgOpTest, MulScalar) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::mul_op mulCalculator;

    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());

    Tensor C = A;
    C.zero();
    Tensor expected = std::move(Tensor({
                                 static_cast<Type>(2.0), static_cast<Type>(4.0), static_cast<Type>(6.0),
                                 static_cast<Type>(0.0), static_cast<Type>(8.0),static_cast<Type>(10.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(12.0)}).to_device<Device>());

    auto alpha = static_cast<Type>(2.0);
    mulCalculator(alpha, A, C);
    EXPECT_EQ(C, expected);
}

TYPED_TEST(LinalgOpTest, Div) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::div_op divCalculator;

    Tensor A = std::move(Tensor({static_cast<Type>(2.0), static_cast<Type>(4.0), static_cast<Type>(6.0),
                                 static_cast<Type>(0.0), static_cast<Type>(8.0), static_cast<Type>(10.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(12.0)}).to_device<Device>());

    Tensor B = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(1.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(1.0), static_cast<Type>(1.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor C = A;
    C.zero();
    Tensor expected = std::move(Tensor({
                                 static_cast<Type>(2.0), static_cast<Type>(2.0), static_cast<Type>(2.0),
                                 static_cast<Type>(0.0), static_cast<Type>(2.0), static_cast<Type>(2.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(2.0)}).to_device<Device>());

    divCalculator(A, B, C);
    EXPECT_EQ(C, expected);

    C.zero();
    C = A / B;
    EXPECT_EQ(C, expected);

    A /= B;
    EXPECT_EQ(A, expected);
}

TYPED_TEST(LinalgOpTest, Transpose) {
    using Type   = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::transpose_op<false> transposeCalculator;

    const int dim = 3;
    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    A.reshape({-1, 3});
    Tensor A_transpose = A;
    A_transpose.zero();

    Tensor expected = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(2.0), static_cast<Type>(4.0), static_cast<Type>(0.0),
                                 static_cast<Type>(3.0), static_cast<Type>(5.0), static_cast<Type>(6.0)}).to_device<Device>());
    expected.reshape({-1, 3});

    transposeCalculator(A, {1, 0}, A_transpose);

    EXPECT_EQ(A_transpose, expected);
}

TYPED_TEST(LinalgOpTest, Stride) {
    using Type   = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::stride_op strideCalculator;

    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());

    Tensor expected = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(4.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor A_stride = expected;
    A_stride.zero();

    strideCalculator(A, {4}, A_stride);

    EXPECT_EQ(A_stride, expected);
}

TYPED_TEST(LinalgOpTest, Inflate) {
    using Type   = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::inflate_op inflateCalculator;

    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor A_inflate = A;
    A_inflate.zero();
    Tensor expected = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(4.0), static_cast<Type>(6.0)}).to_device<Device>());


    inflateCalculator(expected, {4}, A_inflate);

    EXPECT_EQ(A, A_inflate);
}

TYPED_TEST(LinalgOpTest, Reduce) {
    using Type   = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::reduce_op reduceCalculator;

    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    A.reshape({-1, 3});

    Tensor expected = std::move(
        Tensor({static_cast<Type>(6.0), static_cast<Type>(9.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor A_reduce = expected;
    A_reduce.zero();

    reduceCalculator(A, 3, A_reduce);

    EXPECT_EQ(A_reduce, expected);
}

} // namespace op
} // namespace container
