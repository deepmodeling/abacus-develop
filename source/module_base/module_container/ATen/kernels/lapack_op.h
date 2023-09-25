#ifndef ATEN_KERNELS_LAPACK_OP_H_
#define ATEN_KERNELS_LAPACK_OP_H_

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_types.h>

#include <base/third_party/lapack.h>

namespace container {
namespace op {

/**
 * @brief A functor to set a square matrix to upper or lower triangular format.
 *
 * @tparam T The floating-point type used for matrix elements.
 * @tparam Device The type of computing device (not used in the function).
 */
template <typename T, typename Device>
struct set_matrix {
    /**
     * @brief Set a square matrix to upper or lower format.
     *
     * @param uplo Format, 'U' for upper, 'L' for lower.
     * @param A Input matrix A [dim: row x col, column major, lda = lda].
     * @param dim Size of the matrix A.
     */
    void operator() (
        const char& uplo,
        T* A,
        const int& dim);
};


/**
 * @brief Inverse of the triangular factor of a complex matrix.
 *
 * This struct defines an operator that computes the inverse of the triangular factor
 * of a complex matrix. The result is stored in the input matrix 'A', overwriting its content.
 *
 * @tparam T The data type used for the computation (e.g., float, double, etc.).
 * @tparam Device The type of device used for the computation (e.g., CPU or CUDA).
 */
template <typename T, typename Device>
struct lapack_trtri {
    /**
     * @brief Computes the inverse of the triangular factor of a complex matrix.
     *
     * This function computes the inverse of the triangular factor of the given complex matrix A.
     * The result is stored in the input matrix 'A', overwriting its content.
     *
     * @param uplo Specifies whether the matrix is upper or lower triangular ('U' or 'L').
     * @param diag Specifies whether the matrix has unit diagonal ('U' for unit diagonal, 'N' otherwise).
     * @param dim The dimension of the matrix A.
     * @param Mat The complex matrix A (column major) whose triangular factor's inverse will be computed.
     * @param lda Leading dimension of matrix A (stride between elements in the same column).
     * @note This operator is column-major.
     */
    void operator()(
        const char& uplo,
        const char& diag,
        const int& dim,
        T* Mat,
        const int& lda
    );
};


/**
 * @brief Cholesky decomposition for a complex Hermitian-definite matrix.
 *
 * This struct defines an operator that performs Cholesky decomposition on a complex Hermitian-definite matrix.
 * The result is stored in the lower triangular part of the matrix 'A'.
 *
 * @tparam T The data type used for the computation (e.g., float, double, etc.).
 * @tparam Device The type of device used for the computation (e.g., CPU or CUDA).
 */
template <typename T, typename Device>
struct lapack_potrf {
    /**
     * @brief Performs Cholesky decomposition on a complex Hermitian matrix.
     *
     * This function computes the Cholesky decomposition of the given complex Hermitian matrix A.
     * The lower triangular part of the matrix 'A' is overwritten with the result.
     *
     * @param uplo Specifies whether the matrix is stored in the upper or lower triangle ('U' or 'L').
     * @param dim The dimension of the Hermitian matrix A.
     * @param Mat The complex Hermitian matrix A (column major) to be decomposed. The result is stored in A.
     * @param lda Leading dimension of matrix A (stride between elements in the same column).
     */
    void operator()(
        const char& uplo,
        const int& dim,
        T* Mat,
        const int& lda
    );
};


/**
 * @brief Eigenvalue and Eigenvector computation for a complex Hermitian-definite eigenproblem.
 *
 * This struct defines an operator that computes all the eigenvalues and eigenvectors of a complex
 * Hermitian-definite eigenproblem. It uses a divide and conquer algorithm for the eigenvectors.
 * The CPU version is implemented using the `evd` interface, and the CUDA version is implemented
 * through the `cusolverDnZheevd` interface.
 *
 * @tparam T The data type used for the computation (e.g., float, double, etc.).
 * @tparam Device The type of device used for the computation (e.g., CPU or CUDA).
 */
template <typename T, typename Device>
struct lapack_dnevd {
    /**
     * @brief Computes eigenvalues and eigenvectors of a complex Hermitian matrix.
     *
     * This function computes eigenvalues and eigenvectors of the given complex Hermitian matrix A.
     * The eigenvalues are stored in the 'W' array.
     * Check the API documentation for more information about the underlying implementations:
     * https://netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga77dfa610458b6c9bd7db52533bfd53a1.html
     * https://netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_ga9b3e110476166e66f2f62fa1fba6344a.html 
     *
     * @param jobz 'N' to compute eigenvalues only, 'V' to compute eigenvalues and eigenvectors.
     * @param uplo 'U' if the upper triangle of A should be used, 'L' if the lower triangle of A should be used.
     * @param Mat The Hermitian matrix A (col major), where eigenvectors will also be stored.
     * @param dim The dimension of the Hermitian matrix A.
     * @param eigen_val The array to store the calculated eigenvalues.
     * 
     * @note For more information about the underlying implementations, refer to the provided API documentation.
     * @note The type of eigen_val is always real, even when the type of Mat is complex. So we use a struct
     *      GetTypeReal to convert the type of T to real: 
     *      GetTypeReal<std::complex<float>>::type = float
     *      GetTypeReal<std::complex<double>>::type = double
     *      Otherwise, GetTypeReal<T>::type = T
     */
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const char& jobz,
        const char& uplo,
        T* Mat,
        const int& dim,
        Real* eigen_val);
};


/**
 * @brief A struct representing the DNGVD operation for computing eigenvalues and eigenvectors of a complex generalized Hermitian-definite eigenproblem.
 *
 * The DNGVD operation computes all the eigenvalues and eigenvectors of a complex generalized Hermitian-definite eigenproblem.
 * If eigenvectors are desired, it uses a divide and conquer algorithm. The operation supports both CPU and CUDA versions.
 *
 * API Documentation:
 * 1. LAPACK zhegvd: https://netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_ga74fdf9b5a16c90d8b7a589dec5ca058a.html
 * 2. cuSOLVER cusolverDnZhegvd: https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdn-t-sygvd
 *
 * @tparam T The type of the elements in the matrices (e.g., float, double, etc.).
 * @tparam Device The device where the matrices reside (e.g., CPU, GPU, etc.).
 */
template <typename T, typename Device>
struct lapack_dngvd {
    /**
     * @brief Compute all the eigenvalues and eigenvectors of a complex generalized Hermitian-definite eigenproblem.
     *
     * This function computes all the eigenvalues and eigenvectors of a complex generalized Hermitian-definite eigenproblem.
     * The problem is defined as A * x = lambda * B * x, where A is the hermitian matrix, B is the overlap matrix, x is the eigenvector,
     * and lambda is the eigenvalue. The operation supports both CPU and CUDA versions.
     *
     * @param itype Type of generalized problem to solve (1: A*x = lambda*B*x, 2: A*B*x = lambda*x).
     * @param jobz 'N' to compute eigenvalues only, 'V' to compute eigenvalues and eigenvectors.
     * @param uplo 'U' if the upper triangle of matrices A and B should be used, 'L' if the lower triangle should be used.
     * @param Mat_A A pointer to the hermitian matrix A (col major).
     * @param Mat_B A pointer to the overlap matrix B (col major).
     * @param dim The dimension of the matrices A and B.
     * @param eigen_val A pointer to store the calculated eigenvalues.
     *
     * @note The leading dimension ldh should be at least max(1, nstart) for row-major storage and column-major storage.
     * @note The length of the A and B arrays should be at least nstart * ldh.
     * @note The type of eigen_val is always real, even when the type of Mat is complex. So we use a struct
     *      GetTypeReal to convert the type of T to real: 
     *      GetTypeReal<std::complex<float>>::type = float
     *      GetTypeReal<std::complex<double>>::type = double
     *      Otherwise, GetTypeReal<T>::type = T
     */
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const int& itype,
        const char& jobz,
        const char& uplo,
        T* Mat_A,
        T* Mat_B,
        const int& dim,
        Real* eigen_val);
};

#if defined(__CUDA) || defined(__ROCM)
void createGpuSolverHandle();  // create cusolver handle
void destroyGpuSolverHandle(); // destroy cusolver handle
#endif

} // namespace container
} // namespace op

#endif // ATEN_KERNELS_LAPACK_OP_H_