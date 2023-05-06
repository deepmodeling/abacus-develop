#ifndef CONTAINER_KERNELS_LAPACK_OP_H_
#define CONTAINER_KERNELS_LAPACK_OP_H_

#include "../tensor.h"
#include "../tensor_types.h"
#include "third_party/lapack_connector.h"

#if defined(__CUDA) || defined(__UT_USE_CUDA)
#include <cusolverDn.h>
// cuSOLVER API errors
static const char* cusolverGetErrorEnum(cusolverStatus_t error) {
    switch (error) {
        case CUSOLVER_STATUS_SUCCESS:
            return "CUSOLVER_STATUS_SUCCESS";
        case CUSOLVER_STATUS_NOT_INITIALIZED:
            return "CUSOLVER_STATUS_NOT_INITIALIZED";
        case CUSOLVER_STATUS_ALLOC_FAILED:
            return "CUSOLVER_STATUS_ALLOC_FAILED";
        case CUSOLVER_STATUS_INVALID_VALUE:
            return "CUSOLVER_STATUS_INVALID_VALUE";
        case CUSOLVER_STATUS_ARCH_MISMATCH:
            return "CUSOLVER_STATUS_ARCH_MISMATCH";
        case CUSOLVER_STATUS_MAPPING_ERROR:
            return "CUSOLVER_STATUS_MAPPING_ERROR";
        case CUSOLVER_STATUS_EXECUTION_FAILED:
            return "CUSOLVER_STATUS_EXECUTION_FAILED";
        case CUSOLVER_STATUS_INTERNAL_ERROR:
            return "CUSOLVER_STATUS_INTERNAL_ERROR";
        case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
            return "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
        case CUSOLVER_STATUS_NOT_SUPPORTED:
            return "CUSOLVER_STATUS_NOT_SUPPORTED ";
        case CUSOLVER_STATUS_ZERO_PIVOT:
            return "CUSOLVER_STATUS_ZERO_PIVOT";
        case CUSOLVER_STATUS_INVALID_LICENSE:
            return "CUSOLVER_STATUS_INVALID_LICENSE";
        default:
            return "Unknown cusolverStatus_t message";
    }
}

inline void cusolverAssert(cusolverStatus_t code, const char* file, int line, bool abort = true)
{
    if (code != CUSOLVER_STATUS_SUCCESS)
    {
        fprintf(stderr, "cuSOLVER Assert: %s %s %d\n", cusolverGetErrorEnum(code), file, line);
        if (abort)
            exit(code);
    }
}

#define cusolverErrcheck(res) {                    \
        cusolverAssert((res), __FILE__, __LINE__); \
    }
#endif // __CUDA || __UT_USE_CUDA

namespace container {
namespace op {

template <typename T, typename Device>
struct dngvd_op {
    /// @brief DNGVD computes all the eigenvalues and eigenvectors of a complex generalized
    /// Hermitian-definite eigenproblem. If eigenvectors are desired, it uses a divide and conquer algorithm.
    ///
    /// In this op, the CPU version is implemented through the `gvd` interface, and the CUDA version
    /// is implemented through the `gvd` interface.
    /// API doc:
    /// 1. zhegvd: https://netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_ga74fdf9b5a16c90d8b7a589dec5ca058a.html
    /// 2. cusolverDnZhegvd: https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdn-t-sygvd
    ///
    /// Input Parameters
    ///     @param d : the type of device
    ///     @param nstart : the number of cols of the matrix
    ///     @param ldh : the number of rows of the matrix
    ///     @param A : the hermitian matrix A in A x=lambda B x (col major)
    ///     @param B : the overlap matrix B in A x=lambda B x (col major)
    /// Output Parameter
    ///     @param W : calculated eigenvalues
    ///     @param V : calculated eigenvectors (col major)
    void operator()(
            const int nstart,
            const int ldh,
            const std::complex<T>* A,
            const std::complex<T>* B,
            T* W,
            std::complex<T>* V);
};


template <typename T, typename Device>
struct dnevx_op {
    /// @brief DNEVX computes the first m eigenvalues ​​and their corresponding eigenvectors of
    /// a complex generalized Hermitian-definite eigenproblem
    ///
    /// In this op, the CPU version is implemented through the `evx` interface, and the CUDA version
    /// is implemented through the `evd` interface and acquires the first m eigenpairs.
    /// API doc:
    /// 1. zheevx: https://netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_gaabef68a9c7b10df7aef8f4fec89fddbe.html
    /// 2. cusolverDnZheevd: https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdn-t-syevd
    ///
    /// Input Parameters
    ///     @param d : the type of device
    ///     @param nstart : the number of cols of the matrix
    ///     @param ldh : the number of rows of the matrix
    ///     @param A : the hermitian matrix A in A x=lambda B x (row major)
    /// Output Parameter
    ///     @param W : calculated eigenvalues
    ///     @param V : calculated eigenvectors (row major)
    void operator()(
            const int nstart,
            const int ldh,
            const std::complex<T>* A,
            const int m,
            T* W,
            std::complex<T>* V);
};


#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
void createCusolverHandle();
void destroyCusolverHandle();
#endif

} // namespace container
} // namespace op

#endif // CONTAINER_KERNELS_LAPACK_OP_H_