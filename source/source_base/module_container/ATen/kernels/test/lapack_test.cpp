#include <gtest/gtest.h>

#include <ATen/core/tensor.h>
#include <ATen/kernels/lapack.h>
#include <base/utils/gtest.h>

#include <complex>
#include <limits>
#include <stdexcept>
#include <vector>

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

#if defined(__CUDA)
namespace {
template <typename T>
T trtri_value(double real, double imag)
{
    return static_cast<T>(real);
}

template <>
std::complex<float> trtri_value(double real, double imag)
{
    return std::complex<float>(real, imag);
}

template <>
std::complex<double> trtri_value(double real, double imag)
{
    return std::complex<double>(real, imag);
}
} // namespace

template <typename T>
class CudaTrtriTest : public LapackTest<T> {};

using CudaTrtriTypes = ::testing::Types<float, double, std::complex<float>, std::complex<double>>;
TYPED_TEST_SUITE(CudaTrtriTest, CudaTrtriTypes);

TYPED_TEST(CudaTrtriTest, InverseAndUntouchedStorage)
{
    using T = TypeParam;
    using Real = typename GetTypeReal<T>::type;
    for (const int n : {1, 3, 7})
    {
        const int padded = n + 3;
        for (const int lda : {n, padded})
        {
            for (const char uplo : {'U', 'L'})
            {
                for (const char diag : {'N', 'U'})
                {
                    SCOPED_TRACE(::testing::Message() << "n=" << n << " lda=" << lda << " uplo=" << uplo << " diag=" << diag);
                    Tensor input(DataTypeToEnum<T>::value, {lda, n});
                    T* data = input.data<T>();
                    std::vector<std::complex<double>> matrix(n * n, 0.0);
                    for (int col = 0; col < n; ++col)
                    {
                        for (int row = 0; row < lda; ++row)
                        {
                            data[col * lda + row] = T(-9);
                            const bool selected = uplo == 'U' ? row <= col : row >= col;
                            if (row < n && selected)
                            {
                                const double real = row == col ? 2.0 + row : 0.05 * (row + col + 1);
                                const double imag = row == col ? 0.0 : 0.03 * (row - col);
                                const T value = trtri_value<T>(real, imag);
                                const bool implicit = row == col && diag == 'U';
                                // Zero stored diagonals prove that unit-diagonal TRSM ignores them.
                                data[col * lda + row] = implicit ? T(0) : value;
                                matrix[col * n + row] = implicit ? T(1) : value;
                            }
                        }
                    }
                    Tensor device = input.to_device<DEVICE_GPU>();
                    lapack_trtri<T, DEVICE_GPU>()(uplo, diag, n, device.data<T>(), lda);
                    const Tensor result = device.to_device<DEVICE_CPU>();
                    const T* inverse = result.data<T>();
                    std::vector<std::complex<double>> unpacked(n * n, 0.0);
                    for (int col = 0; col < n; ++col)
                    {
                        for (int row = 0; row < lda; ++row)
                        {
                            const bool selected = uplo == 'U' ? row <= col : row >= col;
                            if (row < n && selected)
                            {
                                const bool implicit = row == col && diag == 'U';
                                unpacked[col * n + row] = implicit ? T(1) : inverse[col * lda + row];
                            }
                            else
                            {
                                EXPECT_EQ(inverse[col * lda + row], data[col * lda + row]);
                            }
                        }
                    }
                    const double tolerance = 50 * n * std::numeric_limits<Real>::epsilon();
                    for (int col = 0; col < n; ++col)
                    {
                        for (int row = 0; row < n; ++row)
                        {
                            std::complex<double> left = 0.0;
                            std::complex<double> right = 0.0;
                            // Independent host products, without calling the BLAS implementation under test.
                            for (int k = 0; k < n; ++k)
                            {
                                left += matrix[k * n + row] * unpacked[col * n + k];
                                right += unpacked[k * n + row] * matrix[col * n + k];
                            }
                            const double expected = row == col ? 1.0 : 0.0;
                            EXPECT_LE(std::abs(left - expected), tolerance);
                            EXPECT_LE(std::abs(right - expected), tolerance);
                        }
                    }
                }
            }
        }
    }
}

TYPED_TEST(CudaTrtriTest, RejectsSingularDiagonal)
{
    using T = TypeParam;
    const int n = 3;
    const int lda = 6;
    const char diag = 'N';
    for (const char uplo : {'U', 'L'})
    {
        for (const int zero : {0, n - 1})
        {
            Tensor input(DataTypeToEnum<T>::value, {lda, n});
            input.zero();
            for (int i = 0; i < n; ++i)
            {
                input.data<T>()[i * lda + i] = i == zero ? T(0) : T(2);
            }
            Tensor device = input.to_device<DEVICE_GPU>();
            EXPECT_THROW((lapack_trtri<T, DEVICE_GPU>()(uplo, diag, n, device.data<T>(), lda)), std::runtime_error);
        }
    }
}
#endif

TYPED_TEST(LapackTest, Potrf) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    lapack_potrf<Type, Device> potrfCalculator;
    set_matrix<Type, Device> setMatrixCalculator;

    const int dim = 3;
    Tensor A = std::move(Tensor({static_cast<Type>(4.0), static_cast<Type>(1.0), static_cast<Type>(2.0),
                                 static_cast<Type>(1.0), static_cast<Type>(5.0), static_cast<Type>(3.0),
                                 static_cast<Type>(2.0), static_cast<Type>(3.0), static_cast<Type>(6.0)}).to_device<Device>());

    Tensor B = A;
    // Note all blas and lapack operators within container are column major!
    // For this reason, we should employ 'L' instead of 'U' in the subsequent line.
    potrfCalculator('L', dim, B.data<Type>(), dim);
    // B may live on an accelerator, so pull it back before inspecting elements on the host.
    const Tensor factorized = B.to_device<DEVICE_CPU>();
    EXPECT_GT(std::abs(factorized.data<Type>()[0]), 0.0);
    EXPECT_GT(std::abs(factorized.data<Type>()[4]), 0.0);
    EXPECT_GT(std::abs(factorized.data<Type>()[8]), 0.0);

    setMatrixCalculator('L', B.data<Type>(), dim);
    const Tensor masked = B.to_device<DEVICE_CPU>();
    EXPECT_EQ(masked.data<Type>()[1], static_cast<Type>(0.0));
    EXPECT_EQ(masked.data<Type>()[2], static_cast<Type>(0.0));
    EXPECT_EQ(masked.data<Type>()[5], static_cast<Type>(0.0));
}

TYPED_TEST(LapackTest, GetrfGetriGetrs) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    // This test drives the wrappers with host stack buffers, so it is pinned to the CPU
    // backend; handing these pointers to the cuSolver path would be an invalid device pointer.
    using Device = DEVICE_CPU;

    const int dim = 2;
    const int rhs_count = 1;
    const int workspace_size = 8;
    Type matrix[4] = {static_cast<Type>(4.0),
                      static_cast<Type>(2.0),
                      static_cast<Type>(1.0),
                      static_cast<Type>(3.0)};
    int pivots[dim] = {0};

    lapack_getrf<Type, Device>()(dim, dim, matrix, dim, pivots);

    Type factorized[4] = {matrix[0], matrix[1], matrix[2], matrix[3]};
    Type rhs[2] = {static_cast<Type>(1.0), static_cast<Type>(1.0)};
    lapack_getrs<Type, Device>()('N', dim, rhs_count, factorized, dim, pivots, rhs, dim);
    EXPECT_NEAR(std::abs(rhs[0] - static_cast<Type>(0.2)), 0.0, 1.0e-6);
    EXPECT_NEAR(std::abs(rhs[1] - static_cast<Type>(0.2)), 0.0, 1.0e-6);

    Type workspace[workspace_size];
    lapack_getri<Type, Device>()(dim, matrix, dim, pivots, workspace, workspace_size);
    EXPECT_NEAR(std::abs(matrix[0] - static_cast<Type>(0.3)), 0.0, 1.0e-6);
    EXPECT_NEAR(std::abs(matrix[1] - static_cast<Type>(-0.2)), 0.0, 1.0e-6);
    EXPECT_NEAR(std::abs(matrix[2] - static_cast<Type>(-0.1)), 0.0, 1.0e-6);
    EXPECT_NEAR(std::abs(matrix[3] - static_cast<Type>(0.4)), 0.0, 1.0e-6);
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
