// TODO: This is a temperary location for these functions.
// And will be moved to a global module(module base) later.
#ifndef MODULE_HSOLVER_DNGVD_H
#define MODULE_HSOLVER_DNGVD_H

#include "math_kernel_op.h"
#include "module_base/lapack_connector.h"

namespace hsolver
{

template <typename T, typename Device>
struct dngvd_op
{
    using Real = typename GetTypeReal<T>::type;
    /// @brief DNGVD computes all the eigenvalues and eigenvectors of a complex generalized
    /// Hermitian-definite eigenproblem. If eigenvectors are desired, it uses a divide and conquer algorithm.
    ///
    /// In this op, the CPU version is implemented through the `gvd` interface, and the CUDA version
    /// is implemented through the `gvd` interface.
    /// API doc:
    /// 1. zhegvd:
    /// https://netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_ga74fdf9b5a16c90d8b7a589dec5ca058a.html
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
    void operator()(const Device* d, const int nstart, const int ldh, const T* A, const T* B, Real* W, T* V);
};

template <typename T, typename Device>
struct dngv_op
{
    using Real = typename GetTypeReal<T>::type;
    /// @brief DNGVX computes first m eigenvalues and eigenvectors of a complex generalized
    /// Input Parameters
    ///     @param d : the type of device
    ///     @param nbase : the number of dim of the matrix
    ///     @param ldh : the number of dmx of the matrix
    ///     @param A : the hermitian matrix A in A x=lambda B x (col major)
    ///     @param B : the overlap matrix B in A x=lambda B x (col major)
    /// Output Parameter
    ///     @param W : calculated eigenvalues
    ///     @param V : calculated eigenvectors (col major)
    void operator()(const Device* d, const int nstart, const int ldh, const T* A, T* B, Real* W, T* V);
};

template <typename T, typename Device>
struct dngvx_op
{
    using Real = typename GetTypeReal<T>::type;
    /// @brief DNGVX computes first m eigenvalues and eigenvectors of a complex generalized
    /// Input Parameters
    ///     @param d : the type of device
    ///     @param nbase : the number of dim of the matrix
    ///     @param ldh : the number of dmx of the matrix
    ///     @param A : the hermitian matrix A in A x=lambda B x (col major)
    ///     @param B : the overlap matrix B in A x=lambda B x (col major)
    ///     @param m : the number of eigenpair
    /// Output Parameter
    ///     @param W : calculated eigenvalues
    ///     @param V : calculated eigenvectors (col major)
    void operator()(const Device* d, const int nstart, const int ldh, T* A, T* B, const int m, Real* W, T* V);
};

template <typename T, typename Device>
struct dnevx_op
{
    using Real = typename GetTypeReal<T>::type;
    /// @brief DNEVX computes the first m eigenvalues and their corresponding eigenvectors of
    /// a complex generalized Hermitian-definite eigenproblem
    ///
    /// In this op, the CPU version is implemented through the `evx` interface, and the CUDA version
    /// is implemented through the `evd` interface and acquires the first m eigenpairs.
    /// API doc:
    /// 1. zheevx:
    /// https://netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_gaabef68a9c7b10df7aef8f4fec89fddbe.html
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
    void operator()(const Device* d, const int nstart, const int ldh, const T* A, const int m, Real* W, T* V);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM

void createGpuSolverHandle();
void destroyGpuSolverHandle();

#endif

} // namespace hsolver

#endif // !MODULE_HSOLVER_DNGVD_H