#include <gtest/gtest.h>

#include <ATen/core/tensor.h>
#include <ATen/kernels/lapack.h>
#include <base/utils/gtest.h>

namespace container {
namespace kernels {

template <typename T>
class LapackTest : public testing::Test {
public:
    LapackTest() {
        base::utils::init_blas_handle();
        base::utils::init_cusolver_handle();
    }
    ~LapackTest() override {
        base::utils::delete_blas_handle();
        base::utils::delete_cusolver_handle();
    }
};

TYPED_TEST_SUITE(LapackTest, base::utils::Types);

TYPED_TEST(LapackTest, Trtri) {
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

TYPED_TEST(LapackTest, Potrf) {

    return;
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
    gemmCalculator(transa, transb, m, n, k, &alpha, B.to_device<DEVICE_CPU>().data<Type>(), k, B.to_device<DEVICE_CPU>().data<Type>(), n, &beta, C.to_device<DEVICE_CPU>().data<Type>(), n);

    EXPECT_EQ(A, C);
}

// lapack_geqrf_inplace,
// check that QtQ = I
TYPED_TEST(LapackTest, GeqrfInPlace) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    lapack_geqrf_inplace<Type, Device> geqrfCalculator;

    const int m = 4;
    const int n = 3;  // m >= n，Q is m x n column-orthogonal matrix
    const int lda = m;

    Tensor A_input = std::move(Tensor({
        static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0), static_cast<Type>(4.0),
        static_cast<Type>(5.0), static_cast<Type>(6.0), static_cast<Type>(7.0), static_cast<Type>(8.0),
        static_cast<Type>(9.0), static_cast<Type>(10.0), static_cast<Type>(11.0), static_cast<Type>(12.0)
    }).to_device<Device>());

    Tensor A = A_input; // will be overwritten as Q

    // do geqrf -> get orthogonal Q
    geqrfCalculator(m, n, A.data<Type>(), lda);

    // check on CPU
    Tensor Q = A.to_device<DEVICE_CPU>();
    const Type* Q_data = Q.data<Type>();

    // compute QtQ = Q^T * Q (n x n)
    Tensor QtQ = Q; // std::move(Tensor(std::vector<Type>(n * n, static_cast<Type>(0.0))).to_device<DEVICE_CPU>());
    const Type alpha = static_cast<Type>(1.0);
    const Type beta  = static_cast<Type>(0.0);

    blas_gemm<Type, DEVICE_CPU> gemm;
    gemm('C', 'N',           // Q^T * Q
         n, n, m,            //  n x n
         &alpha,
         Q_data, lda,        // Q^T
         Q_data, lda,        // Q
         &beta,
         QtQ.data<Type>(), n);

    // To print value: first to_device CPU, then print
    // // Test code: print A
    // std::cout << "A = " << std::endl;
    // for (int i = 0; i < m; ++i) {
    //     for (int j = 0; j < n; ++j) {
    //         std::cout << A_input.to_device<DEVICE_CPU>().data<Type>()[i + j * m] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // // Test code: print Q
    // std::cout << "Q = " << std::endl;
    // for (int i = 0; i < m; ++i) {
    //     for (int j = 0; j < n; ++j) {
    //         std::cout << Q.data<Type>()[i + j * m] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // // Test code: print QtQ
    // std::cout << "QtQ = " << std::endl;
    // for (int i = 0; i < n; ++i) {
    //     for (int j = 0; j < n; ++j) {
    //         std::cout << QtQ.data<Type>()[i + j * n] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // check QtQ
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Type expected = (i == j) ? static_cast<Type>(1.0) : static_cast<Type>(0.0);
            EXPECT_NEAR(std::abs(QtQ.data<Type>()[i + j * n]), std::abs(expected), 1e-5)
                << "Q^T * Q not identity at (" << i << "," << j << ")";
        }
    }
}

// Test for lapack_heevd and lapack_heevx:
// Solve a standard eigenvalue problem
// and check that A*V = V*E
TYPED_TEST(LapackTest, heevd) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Real = typename GetTypeReal<Type>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    blas_gemm<Type, Device> gemmCalculator;
    blas_axpy<Type, Device> axpyCalculator;
    lapack_heevd<Type, Device> heevdCalculator;

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
    // heevdCalculator('V', 'U', B.data<Type>(), dim, E.data<Real>());
    heevdCalculator(dim, B.data<Type>(), dim, E.data<Real>());

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

TYPED_TEST(LapackTest, heevx) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Real = typename GetTypeReal<Type>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    blas_gemm<Type, Device> gemmCalculator;
    blas_axpy<Type, Device> axpyCalculator;
    lapack_heevx<Type, Device> heevxCalculator;

    const int dim = 3;
    const int neig = 2;  // Compute first 2 eigenvalues

    Tensor A = std::move(Tensor({static_cast<Type>(4.0), static_cast<Type>(1.0), static_cast<Type>(1.0),
                                 static_cast<Type>(1.0), static_cast<Type>(5.0), static_cast<Type>(3.0),
                                 static_cast<Type>(1.0), static_cast<Type>(3.0), static_cast<Type>(6.0)}).to_device<Device>());

    Tensor E = std::move(Tensor({static_cast<Real>(0.0), static_cast<Real>(0.0)}).to_device<Device>());
    Tensor V = A;
    Tensor expected_C1 = std::move(Tensor({static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                           static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(0.0)}).to_device<Device>());
    Tensor expected_C2 = expected_C1;
    expected_C1.zero();
    expected_C2.zero();

    const char trans = 'N';
    const int m = 3;
    const int n = neig;
    const int k = 3;
    const Type alpha = static_cast<Type>(1.0);
    const Type beta  = static_cast<Type>(0.0);

    // Compute first neig eigenvalues and eigenvectors using heevx
    heevxCalculator(dim, dim, A.data<Type>(), neig, E.data<Real>(), V.data<Type>());

    E = E.to_device<ct::DEVICE_CPU>();
    const Tensor Alpha = std::move(Tensor({
            static_cast<Type>(E.data<Real>()[0]),
            static_cast<Type>(E.data<Real>()[1])}));

    // Check the eigenvalues and eigenvectors
    // A * x = lambda * x for the first neig eigenvectors
    // check that A * V = V * E
    // get A * V
    gemmCalculator(trans, trans, m, n, k, &alpha, A.data<Type>(), m, V.data<Type>(), k, &beta, expected_C1.data<Type>(), m);
    // get V * E
    for (int ii = 0; ii < neig; ii++) {
        axpyCalculator(dim, Alpha.data<Type>() + ii, V.data<Type>() + ii * dim, 1, expected_C2.data<Type>() + ii * dim, 1);
    }

    EXPECT_EQ(expected_C1, expected_C2);
}

// Test for lapack_hegvd and lapack_hegvx
// Solve a generalized eigenvalue problem
// and check that A * v = e * B * v
TYPED_TEST(LapackTest, hegvd) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Real = typename GetTypeReal<Type>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    blas_gemm<Type, Device> gemmCalculator;
    blas_axpy<Type, Device> axpyCalculator;
    lapack_hegvd<Type, Device> hegvdCalculator;

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
    // Note all blas and lapack operators within container are column major!
    // For this reason, we should employ 'L' instead of 'U' in the subsequent line.
    hegvdCalculator(dim, dim, A.data<Type>(), I.data<Type>(), E.data<Real>(), B.data<Type>());

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

TYPED_TEST(LapackTest, hegvx) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Real = typename GetTypeReal<Type>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    blas_gemm<Type, Device> gemmCalculator;
    blas_axpy<Type, Device> axpyCalculator;
    lapack_hegvx<Type, Device> hegvxCalculator;

    const int dim = 3;
    const int neig = 2;  // Compute first 2 eigenvalues

    Tensor A = std::move(Tensor({static_cast<Type>(4.0), static_cast<Type>(1.0), static_cast<Type>(1.0),
                                 static_cast<Type>(1.0), static_cast<Type>(5.0), static_cast<Type>(3.0),
                                 static_cast<Type>(1.0), static_cast<Type>(3.0), static_cast<Type>(6.0)}).to_device<Device>());

    Tensor B = std::move(Tensor({static_cast<Type>(2.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(2.0), static_cast<Type>(0.0),
                                 static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(2.0)}).to_device<Device>());

    Tensor E = std::move(Tensor({static_cast<Real>(0.0), static_cast<Real>(0.0)}).to_device<Device>());
    Tensor V = A;
    Tensor expected_C1 = std::move(Tensor({static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(0.0),
                                           static_cast<Type>(0.0), static_cast<Type>(0.0), static_cast<Type>(0.0)}).to_device<Device>());
    Tensor expected_C2 = expected_C1;
    Tensor C_temp = expected_C1;
    expected_C1.zero();
    expected_C2.zero();

    const char trans = 'N';
    const int m = 3;
    const int n = neig;
    const int k = 3;
    const Type alpha = static_cast<Type>(1.0);
    const Type beta  = static_cast<Type>(0.0);

    // Compute first neig eigenvalues and eigenvectors using hegvx
    hegvxCalculator(dim, dim, A.data<Type>(), B.data<Type>(), neig, E.data<Real>(), V.data<Type>());

    E = E.to_device<ct::DEVICE_CPU>();
    const Tensor Alpha = std::move(Tensor({
            static_cast<Type>(E.data<Real>()[0]),
            static_cast<Type>(E.data<Real>()[1])}));

    // Check the eigenvalues and eigenvectors
    // A * x = lambda * B * x for the first neig eigenvectors
    // check that A * V = E * B * V
    // get A * V
    gemmCalculator(trans, trans, m, n, k, &alpha, A.data<Type>(), m, V.data<Type>(), k, &beta, expected_C1.data<Type>(), m);
    // get E * B * V
    // where B is 2 * eye(3,3)
    // get C_temp = B * V first
    gemmCalculator(trans, trans, m, n, k, &alpha, B.data<Type>(), m, V.data<Type>(), k, &beta, C_temp.data<Type>(), m);
    // then compute C2 = E * B * V
    for (int ii = 0; ii < neig; ii++) {
        axpyCalculator(dim, Alpha.data<Type>() + ii, C_temp.data<Type>() + ii * dim, 1, expected_C2.data<Type>() + ii * dim, 1);
    }

    EXPECT_EQ(expected_C1, expected_C2);
}

} // namespace kernels
} // namespace container
