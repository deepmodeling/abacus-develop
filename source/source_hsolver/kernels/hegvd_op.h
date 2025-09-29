/// This is the module for wrapper of 
/// DeNse Generalized eigenValue (eXtended)
/// HErmitian / SYmmetric

// named HEGVD, actually includes HE/SY GV/GVD/GVX

#ifndef MODULE_HSOLVER_HEGVD_H
#define MODULE_HSOLVER_HEGVD_H

// Note:
// names follow the same style as standard LAPACK APIs:
// -----------------------------------
// he stands for Hermitian
// sy stands for Symmetric
// gv stands for Generalized eigenValue problem
// ev stands for EigenValues
// dn stands for dense, maybe, who knows?
// x stands for compute a subset of the eigenvalues and, optionally,
// their corresponding eigenvectors
// d for all, x for selected
// gv: all, gvd: all/devide-and-conquer, x: selected eigenvalues
// -----------------------------------
// search for docs using the op function name as keywords.

// TODO: This is a temperary location for these functions.
// And will be moved to a global module(module base) later.

#include "source_base/macros.h"
#include "source_base/module_external/lapack_wrapper.h"
#include "source_base/parallel_reduce.h"
#include "source_base/module_device/types.h"

namespace hsolver
{

inline double get_real(const std::complex<double> &x) { return x.real(); }

inline float get_real(const std::complex<float> &x) { return x.real(); }

inline double get_real(const double &x) { return x; }

inline float get_real(const float &x) { return x; }


template <typename T, typename Device>
struct hegvd_op
{
    using Real = typename GetTypeReal<T>::type;
    /// @brief HEGVD computes all the eigenvalues and eigenvectors of a complex generalized
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
struct hegv_op
{
    using Real = typename GetTypeReal<T>::type;
    /// @brief HEGV computes first m eigenvalues and eigenvectors of a complex generalized
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
struct hegvx_op
{
    using Real = typename GetTypeReal<T>::type;
    /// @brief HEGVX computes first m eigenvalues and eigenvectors of a complex generalized
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
struct heevx_op
{
    using Real = typename GetTypeReal<T>::type;
    /// @brief heevx computes the first m eigenvalues and their corresponding eigenvectors of
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
    ///     @param ndim : the size of square matrix
    ///     @param lda : leading dimension of the matrix
    ///     @param A : the hermitian matrix A in A x=lambda x
    ///     @param neig : the number of eigenpairs to be calculated
    /// Output Parameter
    ///     @param w: calculated eigenvalues
    ///     @param z: calculated eigenvectors
    void operator()(const Device *d, const int ndim, const int lda, const T *A, const int neig, Real *w, T *z);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM

void createGpuSolverHandle();
void destroyGpuSolverHandle();

#endif

} // namespace hsolver

#endif // !MODULE_HSOLVER_HEGVD_H