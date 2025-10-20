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
struct lapack_heevd {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const char& jobz,
        const char& uplo,
        T* Mat,
        const int& dim,
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
     * @param jobz  'N': Compute eigenvalues only; 'V': Compute eigenvalues and eigenvectors.
     * @param range 'A': All eigenvalues; 'V': Eigenvalues in the half-open interval (vl, vu]; 'I': Eigenvalues with indices il through iu.
     * @param uplo  'U': Upper triangle of A is stored; 'L': Lower triangle is stored.
     * @param dim   The order of the matrix A. dim >= 0.
     * @param Mat   On entry, the Hermitian matrix A. On exit, it may be overwritten.
     * @param vl    Lower bound of the interval to search for eigenvalues if range == 'V'.
     * @param vu    Upper bound of the interval to search for eigenvalues if range == 'V'.
     * @param il    Index of the smallest eigenvalue to be returned if range == 'I'.
     * @param iu    Index of the largest eigenvalue to be returned if range == 'I'.
     * @param m     Output: The total number of found eigenvalues.
     * @param eigen_val Array to store the computed eigenvalues in ascending order.
     * @param eigen_vec If not nullptr and jobz == 'V', array to store the computed eigenvectors.
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

#if defined(__CUDA) || defined(__ROCM)
// TODO: Use C++ singleton to manage the GPU handles
void createGpuSolverHandle();  // create cusolver handle
void destroyGpuSolverHandle(); // destroy cusolver handle
#endif

} // namespace container
} // namespace kernels

#endif // ATEN_KERNELS_LAPACK_H_
