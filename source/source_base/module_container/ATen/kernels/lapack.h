#ifndef ATEN_KERNELS_LAPACK_H_
#define ATEN_KERNELS_LAPACK_H_

#include "source_base/macros.h"
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


// --- 1. Matrix Decomposition ---
template <typename T, typename Device>
struct lapack_trtri {
    void operator()(
        const char& uplo,
        const char& diag,
        const int& dim,
        T* Mat,
        const int& lda);
};


template <typename T, typename Device>
struct lapack_potrf {
    void operator()(
        const char& uplo,
        const int& dim,
        T* Mat,
        const int& lda);
};

template <typename T, typename Device>
struct lapack_getrf {
    void operator()(
        const int& m,
        const int& n,
        T* Mat,
        const int& lda,
        int* ipiv);
};


template <typename T, typename Device>
struct lapack_getri {
    void operator()(
        const int& n,
        T* Mat,
        const int& lda,
        const int* ipiv,
        T* work,
        const int& lwork);
};

// This is QR factorization in-place
// that will change input Mat A to orthogonal/unitary matrix Q
template <typename T, typename Device>
struct lapack_geqrf_inplace {
    /**
     * @brief Perform in-place QR factorization of a matrix using LAPACK's geqrf function.
     *
     * This function computes the QR factorization of an m-by-n matrix A as A = Q * R,
     * where Q is an orthogonal/unitary matrix and R is an upper triangular matrix.
     * The factorization is performed in-place, meaning the input matrix A will be modified.
     *
     * On exit: A is overwritten with the QR factorization Q orthogonal/unitary matrix
     *
     * @param m The number of rows in the matrix A. m >= 0
     * @param n The number of columns in the matrix A. n >= 0
     * @param A Pointer to the matrix A to be factorized. On exit, contains the QR factorization
     * @param lda The leading dimension of the matrix A. lda >= max(1, m)
     */
    void operator()(
        const int m,
        const int n,
        T *A,
        const int lda);
};

// This is QR factorization
// where [in]Mat will be kept and the results are stored in separate matrix Q
// template <typename T, typename Device>
// struct lapack_geqrf{
//     /**
//      * Perform QR factorization of a matrix using LAPACK's geqrf function.
//      *
//      * @param m The number of rows in the matrix.
//      * @param n The number of columns in the matrix.
//      * @param Mat The matrix to be factorized.
//      *        On exit, the upper triangle contains the upper triangular matrix R,
//      *        and the elements below the diagonal, with the array TAU, represent
//      *        the unitary matrix Q as a product of min(m,n) elementary reflectors.
//      * @param lda The leading dimension of the matrix.
//      * @param tau Array of size min(m,n) containing the Householder reflectors.
//      */
//     void operator()(
//         const int m,
//         const int n,
//         T *Mat,
//         const int lda,
//         T *tau);
// };


// --- 2. Linear System Solvers ---
template <typename T, typename Device>
struct lapack_getrs {
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



// --- 3. Standard & Generalized Eigenvalue ---

// ============================================================================
// Standard Hermitian Eigenvalue Problem Solvers
// ============================================================================
// The following structures (lapack_heevd and lapack_heevx) implement solvers
// for standard Hermitian eigenvalue problems of the form:
//      A * x = lambda * x
// where:
//   - A is a Hermitian matrix
//   - lambda are the eigenvalues to be computed
//   - x are the corresponding eigenvectors
//
// ============================================================================
template <typename T, typename Device>
struct lapack_heevd {
    // !> ZHEEVD computes all eigenvalues and, optionally, eigenvectors of a
    // !> complex Hermitian matrix A.  If eigenvectors are desired, it uses a
    // !> divide and conquer algorithm.
    // !>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
    // !>          orthonormal eigenvectors of the matrix A.
    /**
     * @brief Computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix.
     *
     * This function solves the standard Hermitian eigenvalue problem A*x = lambda*x,
     * where A is a Hermitian matrix. It computes all eigenvalues and optionally
     * the corresponding eigenvectors using a divide and conquer algorithm.
     *
     * @param[in] dim   The order of the matrix A. dim >= 0.
     * @param[in,out] Mat   On entry, the Hermitian matrix A.
     *              On exit, if eigenvectors are computed, A contains the
     *              orthonormal eigenvectors of the matrix A.
     * @param[in] lda   The leading dimension of the array Mat. lda >= max(1, dim).
     * @param[out] eigen_val Array of size at least dim. On normal exit, contains the
     *                  eigenvalues in ascending order.
     *
     * @note
     * See LAPACK ZHEEVD or CHEEVD documentation for more details.
     * The matrix is assumed to be stored in upper or lower triangular form
     * according to the uplo parameter (not shown here but typically passed
     * to the actual implementation).
     */
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const int dim,
        T* Mat,
        const int lda,
        Real* eigen_val);
};

template <typename T, typename Device>
struct lapack_heevx {
    using Real = typename GetTypeReal<T>::type;
    /**
     * @brief Computes selected eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix.
     *
     * This function solves the problem A*x = lambda*x, where A is a Hermitian matrix.
     * It computes a subset of eigenvalues and, optionally, the corresponding eigenvectors.
     *
     * @param dim   The order of the matrix A. dim >= 0.
     * @param lda   The leading dimension of the array Mat. lda >= max(1, dim).
     * @param[in] Mat   On entry, the Hermitian matrix A. On exit, A is kept.
     *                  Only used to provide values of matrix.
     * @param neig  The number of eigenvalues to be found. 0 <= neig <= dim.
     * @param eigen_val On normal exit, the first \p neig elements contain the selected
     *                  eigenvalues in ascending order.
     * @param eigen_vec If eigen_vec is not nullptr, then on exit it contains the
     *                  orthonormal eigenvectors of the matrix A. The eigenvectors are stored in
     *                  the columns of eigen_vec, in the same order as the eigenvalues.
     *
     * @note
     * See LAPACK ZHEEVX or CHEEVX documentation for more details.
     * This routine allocates auxiliary memory inside to prevent input matrix from being destroyed.
     */
    void operator()(
        const int dim,
        const int lda,
        const T *Mat,
        const int neig,
        Real *eigen_val,
        T *eigen_vec);
};


// ============================================================================
// Generalized Hermitian-definite Eigenvalue Problem Solvers
// ============================================================================
// The following structures (lapack_hegvd and lapack_hegvx) implement solvers
// for generalized Hermitian-definite eigenvalue problems of the form:
//      A * x = lambda * B * x
// where:
//   - A is a Hermitian matrix
//   - B is a Hermitian positive definite matrix
//   - lambda are the eigenvalues to be computed
//   - x are the corresponding eigenvectors
//
// ============================================================================

template <typename T, typename Device>
struct lapack_hegvd {
    using Real = typename GetTypeReal<T>::type;
    /**
     * @brief Computes all the eigenvalues and, optionally, the eigenvectors of a complex generalized Hermitian-definite eigenproblem.
     *
     * This function solves the problem A*x = lambda*B*x, where A and B are Hermitian matrices, and B is also positive definite.
     *
     * @param n The order of the matrices Mat_A and Mat_B. n >= 0.
     * @param lda The leading dimension of the arrays Mat_A and Mat_B. lda >= max(1, n).
     * @param Mat_A On entry, the Hermitian matrix A. On exit, it may be overwritten.
     * @param Mat_B On entry, the Hermitian positive definite matrix B. On exit, it may be overwritten.
     * @param eigen_val Array to store the computed eigenvalues in ascending order.
     * @param eigen_vec If not nullptr, array to store the computed eigenvectors.
     *
     * @note
     * See LAPACK ZHEGVD or CHEGVD documentation for more details.
     * This function assumes that A and B have the same leading dimensions, lda.
     * This function copies B to auxiliary memory to avoid being overwritten.
     */
    void operator()(
        const int n,
        const int lda,
        T *Mat_A,
        T *Mat_B,
        Real *eigen_val,
        T *eigen_vec);
};

template <typename T, typename Device>
struct lapack_hegvx {
    using Real = typename GetTypeReal<T>::type;
    /**
     * @ brief hegvx computes the first m eigenvalues and their corresponding eigenvectors of
     * a complex generalized Hermitian-definite eigenproblem.
     *
     * In this op, the CPU version is implemented through the `hegvx` interface, and the CUDA version
     * is implemented through the `evd` interface and acquires the first m eigenpairs
     *
     * hegvx 'V' 'I' 'U'  is used to compute the first m eigenpairs of the problem
     *
     * @param n The order of the matrices A and B. n >= 0.
     * @param lda The leading dimension of the array A and B. lda >= max(1, n).
     * @param A On entry, the Hermitian matrix A. On exit, if info = 0, A contains the matrix Z of eigenvectors.
     * @param B On entry, the Hermitian positive definite matrix B. On exit, the triangular factor from the Cholesky factorization of B.
     * @param m The number of eigenvalues and eigenvectors to be found. 0 < m <= n.
     * @param eigen_val The first m eigenvalues in ascending order.
     * @param eigen_vec The first m columns contain the orthonormal eigenvectors of the matrix A corresponding to the selected eigenvalues.
     *
     * @note
     * See LAPACK ZHEGVX doc for more details.
     * This routine allocates auxiliary memory inside to prevent input matrix from being destroyed.
     */
    void operator()(
        const int n,
        const int lda,
        T *Mat_A,
        T *Mat_B,
        const int m,
        Real *eigen_val,
        T *eigen_vec);
};


#if defined(__CUDA) || defined(__ROCM)
// TODO: Use C++ singleton to manage the GPU handles
void createGpuSolverHandle();  // create cusolver handle
void destroyGpuSolverHandle(); // destroy cusolver handle
#endif

} // namespace container
} // namespace kernels

#endif // ATEN_KERNELS_LAPACK_H_
