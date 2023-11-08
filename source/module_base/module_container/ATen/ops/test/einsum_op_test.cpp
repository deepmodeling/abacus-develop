#include <gtest/gtest.h>

#include <ATen/core/tensor.h>
#include <ATen/ops/einsum_op.h>
#include <test/test_utils.h>

namespace container {
namespace op {

template <typename T>
class EinsumOpTest : public testing::Test {
public:
    EinsumOpTest() {
        test_utils::init_blas_handle();
        test_utils::init_cusolver_handle();
    }
    ~EinsumOpTest() override {
        test_utils::delete_blas_handle();
        test_utils::delete_cusolver_handle();
    }
};

TYPED_TEST_SUITE(EinsumOpTest, test_utils::Types);

TYPED_TEST(EinsumOpTest, Transform) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    const int dim = 3;
    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    A.reshape({-1, dim});
    Tensor expected = std::move(Tensor(
                                {static_cast<Type>(1.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(2.0), static_cast<Type>(4.0), static_cast<Type>(0.0),
                                 static_cast<Type>(3.0), static_cast<Type>(5.0), static_cast<Type>(6.0)}).to_device<Device>());
    expected.reshape({-1, dim});
    // const Tensor expected = std::move(Tensor({static_cast<Type>(21.0)}).to_device<Device>());

    Tensor A_transformed = op::einsum("ij->ji", A);
    EXPECT_EQ(A_transformed, expected);
}

TYPED_TEST(EinsumOpTest, Reduce) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    const int dim = 3;
    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    A.reshape({-1, dim});
    Tensor expected_1 = std::move(Tensor(
                                {static_cast<Type>(6.0), static_cast<Type>(9.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor expected_2 = std::move(Tensor(
                                {static_cast<Type>(1.0), static_cast<Type>(6.0), static_cast<Type>(14.0)}).to_device<Device>());
    // const Tensor expected = std::move(Tensor({static_cast<Type>(21.0)}).to_device<Device>());

    // Case 1: Normal reduction
    Tensor A_reduced = op::einsum("ij->i", A);
    EXPECT_EQ(A_reduced, expected_1);

    // Case 2: Transpose reduction
    A_reduced = op::einsum("ij->j", A);
    EXPECT_EQ(A_reduced, expected_2);

    // Case 3: All reduction
    A_reduced = op::einsum("ij->", A);
    EXPECT_EQ(A_reduced, Tensor({static_cast<Type>(21.0)}).to_device<Device>());
}

TYPED_TEST(EinsumOpTest, Stride) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    const int dim = 3;
    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    A.reshape({-1, dim});
    Tensor expected = std::move(Tensor(
                                {static_cast<Type>(1.0), static_cast<Type>(4.0), static_cast<Type>(6.0)}).to_device<Device>());

    Tensor A_strided = op::einsum("ii->i", A);
    EXPECT_EQ(A_strided, expected);
}

TYPED_TEST(EinsumOpTest, Inflate) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    const int dim = 3;
    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(4.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor expected = std::move(Tensor(
                                {static_cast<Type>(1.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    expected.reshape({-1, dim});

    Tensor A_inflated = op::einsum("i->ii", A);
    EXPECT_EQ(A_inflated, expected);
}

TYPED_TEST(EinsumOpTest, ContractDot) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    const int n = 4;
    const Tensor x = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0), static_cast<Type>(4.0)}).to_device<Device>());
    const Tensor y = std::move(Tensor({static_cast<Type>(4.0), static_cast<Type>(3.0), static_cast<Type>(2.0), static_cast<Type>(1.0)}).to_device<Device>());
    
    const Tensor expected = std::move(Tensor({static_cast<Type>(20.0)}).to_device<Device>());

    Tensor z = op::einsum("i,i->", x, y);
    EXPECT_EQ(z, expected);
}

TYPED_TEST(EinsumOpTest, ContractGemv) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    const int m = 2, n = 4;
    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0), static_cast<Type>(4.0),
                                 static_cast<Type>(5.0), static_cast<Type>(6.0), static_cast<Type>(7.0), static_cast<Type>(8.0)}).to_device<Device>());
    A.reshape({m, n});
    const Tensor x1 = std::move(Tensor(
                                {static_cast<Type>(4.0), static_cast<Type>(3.0), static_cast<Type>(2.0), static_cast<Type>(1.0)}).to_device<Device>());
    const Tensor x2 = std::move(Tensor(
                                {static_cast<Type>(1.0), static_cast<Type>(2.0)}).to_device<Device>());
    
    const Tensor expected_1 = std::move(Tensor(
                                {static_cast<Type>(20.0),static_cast<Type>(60.0)}).to_device<Device>());
    const Tensor expected_2 = std::move(Tensor(
                                {static_cast<Type>(11.0),static_cast<Type>(14.0),static_cast<Type>(17.0), static_cast<Type>(20.0)}).to_device<Device>());

    Tensor y = op::einsum("ij,j->i", A, x1);
    EXPECT_EQ(y, expected_1);
    y = op::einsum("ij,i->j", A, x2);
    EXPECT_EQ(y, expected_2);
}

TYPED_TEST(EinsumOpTest, ContractGemm) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    const int m = 2, k = 4, n = 2;
    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0), static_cast<Type>(4.0),
                                 static_cast<Type>(5.0), static_cast<Type>(6.0), static_cast<Type>(7.0), static_cast<Type>(8.0)}).to_device<Device>());
    A.reshape({m, k});
    Tensor B = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), 
                                 static_cast<Type>(3.0), static_cast<Type>(4.0),
                                 static_cast<Type>(5.0), static_cast<Type>(6.0), 
                                 static_cast<Type>(7.0), static_cast<Type>(8.0)}).to_device<Device>());
    B.reshape({k, n});
    Tensor expected_1 = std::move(Tensor(
                                {static_cast<Type>(50.0), static_cast<Type>(60.0), 
                                 static_cast<Type>(114.0),static_cast<Type>(140.0)}).to_device<Device>());
    expected_1.reshape({m, n});
    Tensor expected_2 = std::move(Tensor(
                                {static_cast<Type>(11.0), static_cast<Type>(23.0), static_cast<Type>(35.0), static_cast<Type>(47.0),
                                 static_cast<Type>(14.0), static_cast<Type>(30.0), static_cast<Type>(46.0), static_cast<Type>(62.0),
                                 static_cast<Type>(17.0), static_cast<Type>(37.0), static_cast<Type>(57.0), static_cast<Type>(77.0),
                                 static_cast<Type>(20.0), static_cast<Type>(44.0), static_cast<Type>(68.0), static_cast<Type>(92.0)}).to_device<Device>());
    expected_2.reshape({k, k});

    Tensor C = op::einsum("ij,jk->ik", A, B);
    EXPECT_EQ(C, expected_1);
    C = op::einsum("ij,ki->jk", A, B);
    EXPECT_EQ(C, expected_2);
}

TYPED_TEST(EinsumOpTest, TransformEllipsis) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    const int dim = 3;
    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0),
                                 static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    A.reshape({-1, dim, dim});
    Tensor expected = std::move(Tensor(
                                {static_cast<Type>(1.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(2.0), static_cast<Type>(4.0), static_cast<Type>(0.0),
                                 static_cast<Type>(3.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                                 static_cast<Type>(1.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(2.0), static_cast<Type>(4.0), static_cast<Type>(0.0),
                                 static_cast<Type>(3.0), static_cast<Type>(5.0), static_cast<Type>(6.0)}).to_device<Device>());
    expected.reshape({-1, dim, dim});
    Tensor expected_ellipsis = std::move(Tensor(
                                {static_cast<Type>(1.0), static_cast<Type>(1.0), 
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), 
                                 static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(2.0), static_cast<Type>(2.0), 
                                 static_cast<Type>(4.0), static_cast<Type>(4.0), 
                                 static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(3.0), static_cast<Type>(3.0), 
                                 static_cast<Type>(5.0), static_cast<Type>(5.0), 
                                 static_cast<Type>(6.0), static_cast<Type>(6.0)}).to_device<Device>());
    expected_ellipsis.reshape({dim, dim, -1});

    Tensor A_transformed = op::einsum("ijk->ikj", A);
    EXPECT_EQ(A_transformed, expected);
    Tensor A_transformed_ellipsis = op::einsum("...ij->...ji", A);
    EXPECT_EQ(A_transformed_ellipsis, expected);
    A_transformed_ellipsis = op::einsum("i...j->j...i", A);
    EXPECT_EQ(A_transformed_ellipsis, expected_ellipsis);
}

TYPED_TEST(EinsumOpTest, ReduceEllipsis) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    const int dim = 3;
    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0),
                                 static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0),
                                 static_cast<Type>(0.0), static_cast<Type>(10.0),static_cast<Type>(11.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(12.0)}).to_device<Device>());
    A.reshape({-1, dim, dim});
    Tensor expected_1 = std::move(Tensor(
                                {static_cast<Type>(6.0), static_cast<Type>(9.0), static_cast<Type>(6.0),
                                 static_cast<Type>(24.0),static_cast<Type>(21.0),static_cast<Type>(12.0)}).to_device<Device>());
    expected_1.reshape({-1, dim});
    Tensor expected_2 = std::move(Tensor(
                                {static_cast<Type>(1.0), static_cast<Type>(6.0), static_cast<Type>(14.0),
                                 static_cast<Type>(7.0), static_cast<Type>(18.0),static_cast<Type>(32.0)}).to_device<Device>());
    expected_2.reshape({-1, dim});

    // Case 1: Normal reduction
    Tensor A_reduced = op::einsum("ijk->ij", A);
    EXPECT_EQ(A_reduced, expected_1);
    Tensor A_reduced_ellipsis = op::einsum("...i->...", A);
    EXPECT_EQ(A_reduced_ellipsis, expected_1);

    // Case 2: Transpose reduction
    A_reduced = op::einsum("ijk->ik", A);
    EXPECT_EQ(A_reduced, expected_2);
    A_reduced_ellipsis =  op::einsum("...jk->...k", A);
    EXPECT_EQ(A_reduced_ellipsis, expected_2);

    // Case 3: All reduction
    A_reduced = op::einsum("ijk->", A);
    EXPECT_EQ(A_reduced, Tensor({static_cast<Type>(78.0)}).to_device<Device>());
    // Not available
    // A_reduced = op::einsum("...->", A);
}

TYPED_TEST(EinsumOpTest, StrideEllipsis) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    const int dim = 3;
    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0),
                                 static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0),
                                 static_cast<Type>(0.0), static_cast<Type>(10.0),static_cast<Type>(11.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(12.0),
                                 static_cast<Type>(13.0),static_cast<Type>(14.0),static_cast<Type>(15.0),
                                 static_cast<Type>(0.0), static_cast<Type>(16.0),static_cast<Type>(17.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(18.0)}).to_device<Device>());
    A.reshape({-1, dim, dim});
    Tensor expected_1 = std::move(Tensor(
                                {static_cast<Type>(1.0), static_cast<Type>(4.0), static_cast<Type>(6.0),
                                 static_cast<Type>(7.0), static_cast<Type>(10.0),static_cast<Type>(12.0),
                                 static_cast<Type>(13.0),static_cast<Type>(16.0),static_cast<Type>(18.0)}).to_device<Device>());
    expected_1.reshape({-1, dim});
    Tensor expected_2 = std::move(Tensor(
                                {static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(10.0),static_cast<Type>(11.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(18.0)}).to_device<Device>());
    expected_2.reshape({-1, dim});

    // Case 1:
    Tensor A_strided = op::einsum("ijj->ij", A);
    EXPECT_EQ(A_strided, expected_1);
    Tensor A_strided_ellipsis = op::einsum("...jj->...j", A);
    EXPECT_EQ(A_strided_ellipsis, expected_1);

    // Case 2:
    A_strided = op::einsum("iij->ij", A);
    EXPECT_EQ(A_strided, expected_2);
    A_strided_ellipsis = op::einsum("ii...->i...", A);
    EXPECT_EQ(A_strided_ellipsis, expected_2);
    A_strided_ellipsis = op::einsum("iij...->ij...", A);
    EXPECT_EQ(A_strided_ellipsis, expected_2);
    A_strided_ellipsis = op::einsum("...iij->ij...", A);
    EXPECT_EQ(A_strided_ellipsis, expected_2);
}

TYPED_TEST(EinsumOpTest, InflateEllipsis) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    const int dim = 3;
    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(4.0), static_cast<Type>(6.0)}).to_device<Device>());
    Tensor expected = std::move(Tensor(
                                {static_cast<Type>(1.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    expected.reshape({-1, dim, dim});

    Tensor A_inflated = op::einsum("i->iii", A);
    EXPECT_EQ(A_inflated, expected);
    Tensor A_inflated_ellipsis = op::einsum("...i->...iii", A);
    EXPECT_EQ(A_inflated_ellipsis, expected);
}

TYPED_TEST(EinsumOpTest, ContractGemmEllipsis) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    const int m = 2, k = 4, n = 2, batch_size = 2;
    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0), static_cast<Type>(4.0),
                                 static_cast<Type>(5.0), static_cast<Type>(6.0), static_cast<Type>(7.0), static_cast<Type>(8.0),
                                 static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0), static_cast<Type>(4.0),
                                 static_cast<Type>(5.0), static_cast<Type>(6.0), static_cast<Type>(7.0), static_cast<Type>(8.0)}).to_device<Device>());
    A.reshape({batch_size, m, k});
    Tensor B = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), 
                                 static_cast<Type>(3.0), static_cast<Type>(4.0),
                                 static_cast<Type>(5.0), static_cast<Type>(6.0), 
                                 static_cast<Type>(7.0), static_cast<Type>(8.0)}).to_device<Device>());
    B.reshape({k, n});
    Tensor expected = std::move(Tensor(
                                {static_cast<Type>(50.0), static_cast<Type>(60.0), 
                                 static_cast<Type>(114.0),static_cast<Type>(140.0),
                                 static_cast<Type>(50.0), static_cast<Type>(60.0), 
                                 static_cast<Type>(114.0),static_cast<Type>(140.0)}).to_device<Device>());
    expected.reshape({batch_size, m, n});

    Tensor C = op::einsum("ijk,...kl->i...jl", A, B);
    EXPECT_EQ(C, expected);

    B = std::move(Tensor({       static_cast<Type>(1.0), static_cast<Type>(2.0), 
                                 static_cast<Type>(3.0), static_cast<Type>(4.0),
                                 static_cast<Type>(5.0), static_cast<Type>(6.0), 
                                 static_cast<Type>(7.0), static_cast<Type>(8.0),
                                 static_cast<Type>(1.0), static_cast<Type>(2.0), 
                                 static_cast<Type>(3.0), static_cast<Type>(4.0),
                                 static_cast<Type>(5.0), static_cast<Type>(6.0), 
                                 static_cast<Type>(7.0), static_cast<Type>(8.0)}).to_device<Device>());
    
    B.reshape({batch_size, k, n});
    C = op::einsum("ijk,ikl->ijl", A, B);
    EXPECT_EQ(C, expected);
    C = op::einsum("...jk,...kl->...jl", A, B);
    EXPECT_EQ(C, expected);
}

} // namespace op
} // namespace container
