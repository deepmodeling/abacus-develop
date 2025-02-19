#ifndef ATEN_KERNELS_LAPACK_H_
#define ATEN_KERNELS_LAPACK_H_

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_types.h>

#include <base/third_party/lapack.h>

namespace container {
namespace kernels {


template <typename T, typename Device>
struct set_matrix {
    void operator() (
        const char& uplo,
        T* A,
        const int& dim);
};


/**
 * @brief A functor for LAPACK TRTRI operation, which computes the inverse of a triangular matrix.
 * 
 * trtri: triangular inverse
 */
template <typename T, typename Device>
struct lapack_trtri {
    /**
     * @brief Computes the inverse of a triangular matrix.
     *
     * @param uplo Specifies whether the matrix is upper or lower triangular.
     *             'U' or 'u' indicates upper triangular, 'L' or 'l' indicates lower triangular.
     * @param diag Specifies whether the matrix is unit triangular.
     *             'N' or 'n' indicates non-unit triangular, 'U' or 'u' indicates unit triangular.
     * @param dim The order of the matrix (number of rows and columns).
     * @param Mat Pointer to the matrix data.
     * @param lda The leading dimension of the matrix.
     */
    void operator()(
        const char& uplo,
        const char& diag,
        const int& dim,
        T* Mat,
        const int& lda);
};


/**
 * @brief Performs the Cholesky decomposition of a given matrix using the LAPACK POTRF routine.
 * The decomposition is performed in-place on the input matrix.
 * 
 */
template <typename T, typename Device>
struct lapack_potrf {
    /**
     * @brief Performs the Cholesky decomposition of a symmetric positive-definite matrix.
     *
     * @param uplo Specifies whether the upper or lower triangular part of the matrix is stored.
     *             'U' for upper triangular, 'L' for lower triangular.
     *          The factorization has the form
     *              A = U^T * U,  if UPLO = 'U', or
     *              A = L * L^T,  if UPLO = 'L',
     *          where U is an upper triangular matrix and L is lower triangular.
     * @param dim The order of the matrix (number of rows and columns).
     * @param Mat Pointer to the matrix data. The matrix is stored in column-major order.
     * @param lda The leading dimension of the matrix. It must be at least max(1, dim).
     */
    void operator()(
        const char& uplo,
        const int& dim,
        T* Mat, 
        const int& lda);
};


/**
 * @brief A wrapper for LAPACK eigenvalue decomposition.
 * 
 * LAPACK
 * - syevd: computes all eigenvalues and, optionally, eigenvectors of a real symmetric matrix.
 * - heevd: computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix.
 *
 */
template <typename T, typename Device>
struct lapack_dnevd {
    using Real = typename GetTypeReal<T>::type;
    /**
     * @brief Perform eigenvalue decomposition.
     *
     * This function performs eigenvalue decomposition on a given matrix.
     *
     * @param jobz Specifies whether to compute eigenvectors ('V') or not ('N').
     * @param uplo Specifies whether the upper ('U') or lower ('L') triangular part of the matrix is used.
     * @param Mat Pointer to the matrix data.
     * @param dim The dimension of the matrix.
     * @param eigen_val Pointer to the array where the computed eigenvalues will be stored.
     */
    void operator()(
        const char& jobz,
        const char& uplo,
        T* Mat,
        const int& dim,
        Real* eigen_val);
};


/**
 * @brief A wrapper for computing the \b generalized eigenvalues and eigenvectors of a pair of matrices using LAPACK.
 * 
 * LAPACK
 * - sygvd: computes all generalized eigenvalues and, optionally, eigenvectors of
 *          a real symmetric-definite generalized eigenproblem.
 * - hegvd: computes all generalized eigenvalues and, optionally, eigenvectors of
 *          a complex Hermitian-definite generalized eigenproblem.
 */
template <typename T, typename Device>
struct lapack_dngvd {
    using Real = typename GetTypeReal<T>::type;
    /**
     * @brief Computes the generalized eigenvalues and eigenvectors of a pair of matrices.
     * 
     * @param itype Specifies the problem type to be solved:
     *              1: A*x = lambda*B*x
     *              2: A*B*x = lambda*x
     *              3: B*A*x = lambda*x
     * @param jobz Specifies whether to compute eigenvectors:
     *             'N': Compute eigenvalues only.
     *             'V': Compute eigenvalues and eigenvectors.
     * @param uplo Specifies whether the upper or lower triangular part of the matrices is used:
     *             'U': Upper triangular part.
     *             'L': Lower triangular part.
     * @param Mat_A Pointer to the first matrix (A).
     * @param Mat_B Pointer to the second matrix (B).
     * @param dim The dimension of the matrices (number of rows/columns).
     * @param eigen_val Pointer to the array where the computed eigenvalues will be stored.
     */
    void operator()(
        const int& itype,
        const char& jobz,
        const char& uplo,
        T* Mat_A,
        T* Mat_B,
        const int& dim,
        Real* eigen_val);
};


/**
 * @brief Functor for performing LU factorization using LAPACK's getrf routine.
 *
 * Performs LU factorization of a general m-by-n matrix using partial pivoting with row interchanges
 * using the LAPACK getrf routine.
 * The factorization is performed in-place on the input matrix.
 */
template <typename T, typename Device>
struct lapack_getrf {
    /**
     * @brief Perform a LU factorization on a matrix.
     * 
     * The factorization has the form
     *     A = P * L * U,
     * where
     * P is a permutation matrix,
     * L is lower triangular with unit diagonal elements, and
     * U is upper triangular.
     *
     * @param m The number of rows of the matrix.
     * @param n The number of columns of the matrix.
     * @param Mat Pointer to the matrix data.
     * @param lda Leading dimension of the matrix.
     * @param ipiv Pointer to the array of pivot indices.
     */
    void operator()(
        const int& m,
        const int& n,
        T* Mat,
        const int& lda,
        int* ipiv);
};


/**
 * @brief Computes the inverse of a matrix using LAPACK's GETRI routine.
 * 
 * GETRI computes the inverse of a matrix using the LU factorization computed by GETRF.
 * 
 * @warning cuSOLVER does \b not provide LU-based matrix inversion interface (getri).
 *          To compute the inverse on \b GPU, use getrs instead.
 */
template <typename T, typename Device>
struct lapack_getri {
    /**
     * @brief Computes the inverse of a matrix.
     * 
     * Inverts U and then computes inv(A) by solving the system inv(A)*L = inv(U) for inv(A).
     * 
     * Perform LU factorization by getrf first!
     * 
     * @param n The order of the matrix (number of rows and columns).
     * @param Mat Pointer to the matrix data stored in column-major order.
     * @param lda The leading dimension of the matrix.
     * @param ipiv Pointer to the pivot indices from the GETRF LU factorization.
     * @param work Pointer to the workspace array.
     * @param lwork The size of the workspace array.
     */
    void operator()(
        const int& n,
        T* Mat,
        const int& lda,
        const int* ipiv,
        T* work,
        const int& lwork);
};

/**
 * @brief A functor for solving a system of linear equations using LAPACK's getrs function.
 *
 * GETRS: triangular solve using factor
 * Solves a system of linear equations
 *      A * X = B or A^T * X = B
 * with a general N-by-N matrix A using the LU factorization computed by getrf.
 *
 */
template <typename T, typename Device>
struct lapack_getrs {
    /**
     * @brief Solves a system of linear equations with a general N-by-N matrix A.
     *
     * @param trans Transpose operation.
     *              'N' or 'n' indicates no transpose.
     *              'T' or 't' indicates transpose.
     * @param n Number of rows of the matrix A.
     * @param nrhs Number of right-hand sides, i.e., the number of columns of the matrix B.
     * @param A Pointer to the matrix A.
     * @param lda leading dimension of the matrix A.
     * @param ipiv Pointer to the pivot indices from the getrf LU factorization.
     * @param B Pointer to the matrix B.
     * @param ldb Leading dimension of the matrix B.
     */
    void operator()(
        const char& trans,
        const int& n,
        const int& nrhs,
        T* A,
        const int& lda,
        const int* ipiv,
        T* B,
        const int& ldb);
};

/**
 * @brief Functor for performing QR factorization using LAPACK's GEQRF routine.
 * 
 * GEQRF: QR factorization
 * Computes the QR factorization of a general m-by-n matrix A using Householder reflectors.
 * The factorization has the form
 *    A = Q * R,
 * where Q is orthogonal and R is upper triangular.
 */
template <typename T, typename Device>
struct lapack_geqrf {
    /**
     * @brief Perform QR factorization on a matrix.
     * 
     * The factorization has the form
     *    A = Q * R,
     * where Q is orthogonal and R is upper triangular.
     * 
     * @param m The number of rows of the matrix.
     * @param n The number of columns of the matrix.
     * @param A Pointer to the matrix data.
     * @param lda Leading dimension of the matrix.
     * @param tau Pointer to the array of scalar factors of the elementary reflectors.
     * @param work Pointer to the workspace array.
     * @param lwork The size of the workspace array.
     */
    void operator()(
        int m,
        int n,
        T* A,
        int lda,
        T* tau,
        T* work,
        int lwork);
};

/**
 * @brief Functor for solving a system of linear equations with a triangular matrix using LAPACK's TRSM routine.
 * 
 * TRSM: triangular solve matrix
 * Solves one of the following matrix equations:
 * - op(A) * X = alpha * B
 * - X * op(A) = alpha * B
 * where op(A) is either A, A^T, A^H, or A^T.
 * 
 */
template <typename T, typename Device>
struct lapack_trsm {
    /**
     * @brief Solve a system of linear equations with a triangular matrix.
     * 
     * Solves one of the following matrix equations:
     * - op(A) * X = alpha * B
     * - X * op(A) = alpha * B
     * where op(A) is either A, A^T, A^H, or A^T.
     *
     * @param side Specifies whether op(A) multiplies B from the left or right.
     *             'L' or 'l' for left, 'R' or 'r' for right.
     * @param uplo Specifies whether the matrix A is an upper or lower triangular matrix.
     *             'U' or 'u' for upper, 'L' or 'l' for lower.
     * @param transA Specifies the form of op(A) to be used in the matrix multiplication.
     *               'N' or 'n' for no transpose, 'T' or 't' for transpose, 'C' or 'c' for conjugate transpose.
     * @param diag Specifies whether or not A is unit triangular.
     *             'U' or 'u' for unit triangular, 'N' or 'n' for non-unit triangular.
     * @param m The number of rows of the matrix B. m >= 0.
     * @param n The number of columns of the matrix B. n >= 0.
     * @param alpha Scalar multiplier applied to op(A) * B.
     * @param A Pointer to the matrix A.
     * @param lda Leading dimension of A. lda >= max(1, m) if side == 'L' or lda >= max(1, n) if side == 'R'.
     * @param B Pointer to the matrix B.
     * @param ldb Leading dimension of B. ldb >= max(1, m).
     */
    void operator()(
        char side,
        char uplo,
        char transA,
        char diag,
        int m,
        int n,
        T alpha,
        T* A,
        int lda,
        T* B,
        int ldb);
};


#if defined(__CUDA) || defined(__ROCM)
// TODO: Use C++ singleton to manage the GPU handles
void createGpuSolverHandle();  // create cusolver handle
void destroyGpuSolverHandle(); // destroy cusolver handle
#endif

} // namespace container
} // namespace kernels

#endif // ATEN_KERNELS_LAPACK_H_