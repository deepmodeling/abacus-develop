// TODO: This is a temperary location for these functions.
// And will be moved to a global module(module base) later.
#ifndef MODULE_HSOLVER_DNGVD_H
#define MODULE_HSOLVER_DNGVD_H

#include "module_base/complexmatrix.h"
#include "module_base/lapack_connector.h"
#include "module_hsolver/include/math_kernel.h"
#include "module_psi/include/memory.h"

namespace hsolver
{

template <typename FPTYPE, typename Device> struct dngvx_op
{
    /// @brief DNGVX computes the first m eigenvalues ​​and their corresponding eigenvectors of
    /// a complex generalized Hermitian-definite eigenproblem
    ///
    /// In this op, the CPU version is implemented through the `gvx` interface, and the CUDA version
    /// is implemented through the `gvd` interface and acquires the first m eigenpairs.
    /// API doc:
    /// 1. zhegvx: https://netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_ga8ea76cbbb14edb5a22069e203fc8e8b2.html
    /// 2. cusolverDnZhegvd: https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdn-t-sygvd
    ///
    /// Input Parameters
    ///     @param d : the type of device
    ///     @param nstart : the number of cols of the matrix
    ///     @param ldh : the number of rows of the matrix
    ///     @param A : the hermitian matrix A in A x=lambda B x (row major)
    ///     @param B : the overlap matrix B in A x=lambda B x (row major)
    ///     @param m : the number of the first m eigenvalues to calculate
    /// Output Parameter
    ///     @param W : calculated eigenvalues
    ///     @param V : calculated eigenvectors (row major)
    void operator()(const Device* d,
                    const int nstart,
                    const int ldh,
                    const std::complex<FPTYPE>* A,
                    const std::complex<FPTYPE>* B,
                    const int m,
                    double* W,
                    std::complex<FPTYPE>* V);
};

template <typename FPTYPE, typename Device> struct dngv_op
{
    /// @brief DNGV computes all the eigenvalues and eigenvectors of a complex generalized
    /// Hermitian-definite eigenproblem
    ///
    /// In this op, the CPU version is implemented through the `gv` interface, and the CUDA version
    /// is implemented through the `gvd` interface.
    /// API doc:
    /// 1. zhegv: https://netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_gaf7b790b3b89de432a423c9006c1cc1ac.html
    /// 2. cusolverDnZhegvd: https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdn-t-sygvd
    ///
    /// Input Parameters
    ///     @param d : the type of device
    ///     @param nstart : the number of cols of the matrix
    ///     @param ldh : the number of rows of the matrix
    ///     @param A : the hermitian matrix A in A x=lambda B x (row major)
    ///     @param B : the overlap matrix B in A x=lambda B x (row major)
    /// Output Parameter
    ///     @param W : calculated eigenvalues
    ///     @param V : calculated eigenvectors (row major)
    void operator()(const Device* d,
                    const int nstart,
                    const int ldh,
                    const std::complex<FPTYPE>* A,
                    const std::complex<FPTYPE>* B,
                    double* W,
                    std::complex<FPTYPE>* V);
};

template <typename FPTYPE, typename Device> struct dngvd_op
{
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
    ///     @param A : the hermitian matrix A in A x=lambda B x (row major)
    ///     @param B : the overlap matrix B in A x=lambda B x (row major)
    /// Output Parameter
    ///     @param W : calculated eigenvalues
    ///     @param V : calculated eigenvectors (row major)
    void operator()(const Device* d,
                    const int nstart,
                    const int ldh,
                    const std::complex<FPTYPE>* A,
                    const std::complex<FPTYPE>* B,
                    double* W,
                    std::complex<FPTYPE>* V);
};


template <typename FPTYPE, typename Device> struct dnevx_op
{
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
    void operator()(const Device* d,
                    const int nstart,
                    const int ldh,
                    const std::complex<FPTYPE>* A,
                    const int m,
                    double* W,
                    std::complex<FPTYPE>* V);
};


#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM

void createCUSOLVERhandle();
void destoryCUSOLVERhandle();

#endif

} // namespace hsolver

#endif // !MODULE_HSOLVER_DNGVD_H