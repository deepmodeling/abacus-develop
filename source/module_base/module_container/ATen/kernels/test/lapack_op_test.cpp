#include <gtest/gtest.h>

#include <ATen/core/tensor.h>
#include <ATen/kernels/lapack_op.h>
#include <ATen/kernels/test/op_test_utils.h>

namespace container {
namespace op {

template <typename T>
class LapackOpTest : public testing::Test {
public:
    LapackOpTest() {
        test_utils::init_blas_handle();
        test_utils::init_cusolver_handle();
    }
    ~LapackOpTest() override {
        test_utils::delete_blas_handle();
        test_utils::delete_cusolver_handle();
    }
};

TYPED_TEST_SUITE(LapackOpTest, test_utils::Types);

TYPED_TEST(LapackOpTest, Trtri) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = DEVICE_CPU;

    blas_gemm<Type, Device> gemmCalculator;
    lapack_trtri<Type, Device> trtriCalculator;

    const int dim = 3;
    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0),
                                 static_cast<Type>(0.0), static_cast<Type>(4.0), static_cast<Type>(5.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(6.0)}).to_device<Device>());
    
    Tensor I = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(1.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(1.0)}).to_device<Device>());
    Tensor B = A;
    Tensor C = B;
    C.zero();
    
    const char trans = 'N';
    const int m = 3;
    const int n = 3;
    const int k = 3;
    const Type alpha = static_cast<Type>(1.0);
    const Type beta  = static_cast<Type>(0.0);
    // Note all blas and lapack operators within container are column major!
    // For this reason, we should employ 'L' instead of 'U' in the subsequent line.
    trtriCalculator('L', 'N', dim, B.data<Type>(), dim);
    gemmCalculator(trans, trans, m, n, k, &alpha, B.data<Type>(), k, A.data<Type>(), n, &beta, C.data<Type>(), n);
    
    EXPECT_EQ(C, I);
}

TYPED_TEST(LapackOpTest, Potrf) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    blas_gemm<Type, Device> gemmCalculator;
    lapack_potrf<Type, Device> potrfCalculator;
    set_matrix<Type, Device> setMatrixCalculator;

    const int dim = 3;
    Tensor A = std::move(Tensor({static_cast<Type>(4.0), static_cast<Type>(1.0), static_cast<Type>(2.0),
                                 static_cast<Type>(1.0), static_cast<Type>(5.0), static_cast<Type>(3.0),
                                 static_cast<Type>(2.0), static_cast<Type>(3.0), static_cast<Type>(6.0)}).to_device<Device>());
    
    Tensor B = A;
    Tensor C = B;
    C.zero();
    
    const char transa = 'N';
    const char transb = 'C';
    const int m = 3;
    const int n = 3;
    const int k = 3;
    const Type alpha = static_cast<Type>(1.0);
    const Type beta  = static_cast<Type>(0.0);
    // Note all blas and lapack operators within container are column major!
    // For this reason, we should employ 'L' instead of 'U' in the subsequent line.
    potrfCalculator('L', dim, B.data<Type>(), dim);
    // Keep the upper triangle of B
    setMatrixCalculator('U', B.data<Type>(), dim);
    // A = U**T * U
    gemmCalculator(transa, transb, m, n, k, &alpha, B.data<Type>(), k, B.data<Type>(), n, &beta, C.data<Type>(), n);

    EXPECT_EQ(A, C);
}

TYPED_TEST(LapackOpTest, dnevd) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Real = typename PossibleComplexToReal<Type>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;
    
    blas_gemm<Type, Device> gemmCalculator;
    blas_axpy<Type, Device> axpyCalculator;
    lapack_dnevd<Type, Device> dnevdCalculator;

    const int dim = 3;
    Tensor A = std::move(Tensor({static_cast<Type>(4.0), static_cast<Type>(1.0), static_cast<Type>(1.0),
                                 static_cast<Type>(1.0), static_cast<Type>(5.0), static_cast<Type>(3.0),
                                 static_cast<Type>(1.0), static_cast<Type>(3.0), static_cast<Type>(6.0)}).to_device<Device>());
    
    Tensor E = std::move(Tensor({static_cast<Real>(0.0), static_cast<Real>(0.0), static_cast<Real>(0.0)}).to_device<Device>());
    Tensor B = A;
    Tensor expected_C1 = A;
    Tensor expected_C2 = A;
    expected_C1.zero();
    expected_C2.zero();
    
    const char trans = 'N';
    const int m = 3;
    const int n = 3;
    const int k = 3;
    const Type alpha = static_cast<Type>(1.0);
    const Type beta  = static_cast<Type>(0.0);
    // Note all blas and lapack operators within container are column major!
    // For this reason, we should employ 'L' instead of 'U' in the subsequent line.
    dnevdCalculator('V', 'U', B.data<Type>(), dim, E.data<Real>());
    
    E = E.to_device<DEVICE_CPU>();
    const Tensor Alpha = std::move(Tensor({
            static_cast<Type>(E.data<Real>()[0]), 
            static_cast<Type>(E.data<Real>()[1]),
            static_cast<Type>(E.data<Real>()[2])}));

    // Check the eigenvalues and eigenvectors
    // A * x = lambda * x
    gemmCalculator(trans, trans, m, n, k, &alpha, A.data<Type>(), m, B.data<Type>(), k, &beta, expected_C1.data<Type>(), m);
    for (int ii = 0; ii < dim; ii++) {
        axpyCalculator(dim, Alpha.data<Type>() + ii, B.data<Type>() + ii * dim, 1, expected_C2.data<Type>() + ii * dim, 1);
    }
    EXPECT_EQ(expected_C1, expected_C2);
}


TYPED_TEST(LapackOpTest, dngvd) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Real = typename PossibleComplexToReal<Type>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;
    
    blas_gemm<Type, Device> gemmCalculator;
    blas_axpy<Type, Device> axpyCalculator;
    lapack_dngvd<Type, Device> dngvdCalculator;

    const int dim = 3;
    Tensor A = std::move(Tensor({static_cast<Type>(4.0), static_cast<Type>(1.0), static_cast<Type>(1.0),
                                 static_cast<Type>(1.0), static_cast<Type>(5.0), static_cast<Type>(3.0),
                                 static_cast<Type>(1.0), static_cast<Type>(3.0), static_cast<Type>(6.0)}).to_device<Device>());
    
    Tensor I = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(1.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(1.0)}).to_device<Device>());
    
    Tensor E = std::move(Tensor({static_cast<Real>(0.0), static_cast<Real>(0.0), static_cast<Real>(0.0)}).to_device<Device>());
    Tensor B = A;
    Tensor expected_C1 = A;
    Tensor expected_C2 = A;
    expected_C1.zero();
    expected_C2.zero();
    
    const char trans = 'N';
    const int m = 3;
    const int n = 3;
    const int k = 3;
    const Type alpha = static_cast<Type>(1.0);
    const Type beta  = static_cast<Type>(0.0);
    // Note al<l blas and lapack operators within container are column major!
    // For this reason, we should employ 'L' instead of 'U' in the subsequent line.
    dngvdCalculator(1, 'V', 'U', B.data<Type>(), I.data<Type>(), dim, E.data<Real>());

    E = E.to_device<DEVICE_CPU>();
    const Tensor Alpha = std::move(Tensor({
            static_cast<Type>(E.data<Real>()[0]), 
            static_cast<Type>(E.data<Real>()[1]),
            static_cast<Type>(E.data<Real>()[2])}));

    // Check the eigenvalues and eigenvectors
    // A * x = lambda * x
    gemmCalculator(trans, trans, m, n, k, &alpha, A.data<Type>(), m, B.data<Type>(), k, &beta, expected_C1.data<Type>(), m);
    for (int ii = 0; ii < dim; ii++) {
        axpyCalculator(dim, Alpha.data<Type>() + ii, B.data<Type>() + ii * dim, 1, expected_C2.data<Type>() + ii * dim, 1);
    }
    EXPECT_EQ(expected_C1, expected_C2);
}

} // namespace op
} // namespace container
