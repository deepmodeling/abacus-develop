// TODO: This is a temperary location for these functions.
// And will be moved to a global module(module base) later.
#ifndef MODULE_HSOLVER_MATH_KERNEL_H
#define MODULE_HSOLVER_MATH_KERNEL_H

#include "module_base/blas_connector.h"
#include "module_psi/psi.h"
#include "module_base/parallel_reduce.h"
#include "module_psi/kernels/memory_op.h"
#include <module_base/macros.h>


#if defined(__CUDA) || defined(__UT_USE_CUDA)
#include <cuda_runtime.h>
#include "cublas_v2.h"

#define cublasErrcheck(res) { cublasAssert((res), __FILE__, __LINE__); }

static const char *_cublasGetErrorEnum(cublasStatus_t error) {
    switch (error) {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";

        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";

        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";

        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";

        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";

        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";

        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";

        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
    }
    return "<unknown>";
}

inline void cublasAssert(cublasStatus_t code, const char *file, int line, bool abort=true) {
    if (code != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr,"cuBLAS Assert: %s %s %d\n", _cublasGetErrorEnum(code), file, line);
        if (abort) exit(code);
    }
}
#endif //__CUDA || __UT_USE_CUDA

namespace hsolver
{

template <typename T, typename Device>
struct line_minimize_with_block_op {
    /// @brief dot_real_op computes the dot product of the given complex arrays(treated as float arrays).
    /// And there's may have MPI communications while enabling planewave parallization strategy.
    ///
    /// Input Parameters
    /// \param dev : the type of computing device
    /// \param A : input array arr
    /// \param dim : size of A
    /// \param lda : leading dimention of A
    /// \param batch : batch size, the size of the result array res
    ///
    /// \return res : the result vector
    /// T : dot product result
    void operator() (
        T* grad_out,
        T* hgrad_out,
        T* psi_out,
        T* hpsi_out,
        const int &n_basis,
        const int &n_basis_max,
        const int &n_band);
};


template <typename T, typename Device>
struct calc_grad_with_block_op {
    /// @brief dot_real_op computes the dot product of the given complex arrays(treated as float arrays).
    /// And there's may have MPI communications while enabling planewave parallization strategy.
    ///
    /// Input Parameters
    /// \param dev : the type of computing device
    /// \param A : input array arr
    /// \param dim : size of A
    /// \param lda : leading dimention of A
    /// \param batch : batch size, the size of the result array res
    ///
    /// \return res : the result vector
    /// T : dot product result
    using Real = typename GetTypeReal<T>::type;
    void operator() (
        const Real* prec_in,
        Real* err_out,
        Real* beta_out,
        T* psi_out,
        T* hpsi_out,
        T* grad_out,
        T* grad_old_out,
        const int &n_basis,
        const int &n_basis_max,
        const int &n_band);
};

template <typename T, typename Device>
struct dot_real_op {
    using Real = typename GetTypeReal<T>::type;
    /// @brief dot_real_op computes the dot product of the given complex arrays(treated as float arrays).
  /// And there's may have MPI communications while enabling planewave parallization strategy.
  ///
  /// Input Parameters
  /// \param d : the type of computing device
  /// \param dim : array size
  /// \param psi_L : input array A
  /// \param psi_R : input array B
  /// \param reduce : flag to control whether to perform the MPI communications
  ///
  /// \return
  /// FPTYPE : dot product result
    Real operator() (
      const Device* d,
      const int& dim,
        const T* psi_L,
        const T* psi_R,
      const bool reduce = true);
};

// vector operator: result[i] = vector[i] / constant
template <typename T, typename Device> struct vector_div_constant_op
{
    using Real = typename GetTypeReal<T>::type;
    /// @brief result[i] = vector[i] / constant
    ///
    /// Input Parameters
    /// \param d : the type of computing device
    /// \param dim : array size
    /// \param vector : input array
    /// \param constant : input constant
    ///
    /// Output Parameters
    /// \param result : output array
    void operator()(const Device* d,
                    const int dim,
        T* result,
        const T* vector,
        const Real constant);
};

// replace vector_div_constant_op : x = alpha * x
template <typename FPTYPE, typename Device> struct scal_op
{
    /// @brief x = alpha * x
    ///
    /// Input Parameters
    /// \param d : the type of computing device
    /// \param N : array size
    /// \param alpha : input constant
    /// \param X : input array
    /// \param incx : computing strip of array X
    ///
    /// Output Parameters
    /// \param X : output array
    void operator()(const Device* d,
                    const int& N,
                    const std::complex<FPTYPE>* alpha,
                    std::complex<FPTYPE>* X,
                    const int& incx);
};

// vector operator: result[i] = vector1[i](complex) * vector2[i](not complex)
template <typename T, typename Device> struct vector_mul_vector_op
{
    using Real = typename GetTypeReal<T>::type;
    /// @brief result[i] = vector1[i](complex) * vector2[i](not complex)
    ///
    /// Input Parameters
    /// \param d : the type of computing device
    /// \param dim : array size
    /// \param vector1 : input array A
    /// \param vector2 : input array B
    ///
    /// Output Parameters
    /// \param result : output array
    void operator()(const Device* d,
        const int& dim,
        T* result,
        const T* vector1,
        const Real* vector2);
};

// vector operator: result[i] = vector1[i](complex) / vector2[i](not complex)
template <typename T, typename Device> struct vector_div_vector_op
{
    using Real = typename GetTypeReal<T>::type;
    /// @brief result[i] = vector1[i](complex) / vector2[i](not complex)
    ///
    /// Input Parameters
    /// \param d : the type of computing device
    /// \param dim : array size
    /// \param vector1 : input array A
    /// \param vector2 : input array B
    ///
    /// Output Parameters
    /// \param result : output array
    void operator()(const Device* d,
        const int& dim,
        T* result,
        const T* vector1,
        const Real* vector2);
};

// vector operator: result[i] = vector1[i] * constant1 + vector2[i] * constant2
template <typename FPTYPE, typename Device> struct constantvector_addORsub_constantVector_op
{
    /// @brief result[i] = vector1[i] * constant1 + vector2[i] * constant2
    ///
    /// Input Parameters
    /// \param d : the type of computing device
    /// \param dim : array size
    /// \param vector1 : input array A
    /// \param constant1 : input constant a
    /// \param vector2 : input array B
    /// \param constant2 : input constant b
    ///
    /// Output Parameters
    /// \param result : output array
    void operator()(const Device* d,
                    const int& dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector1,
                    const FPTYPE constant1,
                    const std::complex<FPTYPE>* vector2,
                    const FPTYPE constant2);
};

//  compute Y = alpha * X + Y
template <typename FPTYPE, typename Device> struct axpy_op
{
    /// @brief Y = alpha * X + Y
    ///
    /// Input Parameters
    /// \param d : the type of computing device
    /// \param N : array size
    /// \param alpha : input constant alpha
    /// \param X : input array X
    /// \param incX : computing strip of X
    /// \param Y : computing strip of Y
    /// \param incY : computing strip of Y
    ///
    /// Output Parameters
    /// \param Y : output array Y
    void operator()(const Device* d,
                    const int& N,
                    const std::complex<FPTYPE>* alpha,
                    const std::complex<FPTYPE>* X,
                    const int& incX,
                    std::complex<FPTYPE>* Y,
                    const int& incY);
};

// compute y = alpha * op(A) * x + beta * y
template <typename T, typename Device> struct gemv_op
{
    /// @brief y = alpha * op(A) * x + beta * y
    ///
    /// Input Parameters
    /// \param d : the type of computing device
    /// \param trans : whether to transpose A
    /// \param m : first dimension of matrix
    /// \param n : second dimension of matrix
    /// \param alpha : input constant alpha
    /// \param A : input matrix A
    /// \param lda : leading dimention of A
    /// \param X : input array X
    /// \param incx : computing strip of X
    /// \param beta : input constant beta
    /// \param Y : input array Y
    /// \param incy : computing strip of Y
    ///
    /// Output Parameters
    /// \param Y : output array Y
    void operator()(const Device* d,
                    const char& trans,
                    const int& m,
                    const int& n,
        const T* alpha,
        const T* A,
                    const int& lda,
        const T* X,
                    const int& incx,
        const T* beta,
        T* Y,
                    const int& incy);
};


// compute C = alpha * op(A) * op(B) + beta * C
template <typename T, typename Device> struct gemm_op
{
    /// @brief C = alpha * op(A) * op(B) + beta * C
    ///
    /// Input Parameters
    /// \param d : the type of computing device
    /// \param transa : whether to transpose matrix A
    /// \param transb : whether to transpose matrix B
    /// \param m : first dimension of matrix mulplication
    /// \param n : second dimension of matrix mulplication
    /// \param k : third dimension of matrix mulplication
    /// \param alpha : input constant alpha
    /// \param a : input matrix A
    /// \param lda : leading dimention of A
    /// \param b : input matrix B
    /// \param ldb : leading dimention of A
    /// \param beta : input constant beta
    /// \param c : input matrix C
    /// \param ldc : leading dimention of C
    ///
    /// Output Parameters
    /// \param c : output matrix C
    void operator()(const Device* d,
        const char& transa,
        const char& transb,
        const int& m,
        const int& n,
        const int& k,
        const T* alpha,
        const T* a,
        const int& lda,
        const T* b,
        const int& ldb,
        const T* beta,
        T* c,
        const int& ldc);
};

template <typename T, typename Device> struct matrixTranspose_op
{
    /// @brief transpose the input matrix
    ///
    /// Input Parameters
    /// \param d : the type of computing device
    /// \param row : first dimension of matrix
    /// \param col : second dimension of matrix
    /// \param input_matrix : input matrix
    ///
    /// Output Parameters
    /// \param output_matrix : output matrix
    void operator()(const Device* d,
        const int& row,
        const int& col,
        const T* input_matrix,
        T* output_matrix);
};

template <typename T, typename Device> struct matrixSetToAnother
{
    /// @brief initialize matrix B with A
    ///
    /// Input Parameters
    /// \param d : the type of computing device
    /// \param n : first dimension of matrix
    /// \param A : input matrix A
    /// \param LDA : leading dimension of A
    /// \param LDB : leading dimension of B
    ///
    /// Output Parameters
    /// \param B : output matrix B
    void operator()(const Device* d,
        const int& n,
        const T* A,
        const int& LDA,
        T* B,
        const int& LDB);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM

template <typename T>
struct line_minimize_with_block_op<T, psi::DEVICE_GPU> {
  using Real = typename GetTypeReal<T>::type;
  void operator()(
    T* grad_out,
    T* hgrad_out,
    T* psi_out,
    T* hpsi_out,
    const int &n_basis,
    const int &n_basis_max,
    const int &n_band);
};

template <typename T>
struct calc_grad_with_block_op<T, psi::DEVICE_GPU> {
  using Real = typename GetTypeReal<T>::type;
  void operator()(
    const Real* prec_in,
    Real* err_out,
    Real* beta_out,
    T* psi_out,
    T* hpsi_out,
    T* grad_out,
    T* grad_old_out,
    const int &n_basis,
    const int &n_basis_max,
    const int &n_band);
};

// Partially specialize functor for psi::GpuDevice.
template <typename T>
struct dot_real_op<T, psi::DEVICE_GPU> {
    using Real = typename GetTypeReal<T>::type;
    Real operator()(
        const psi::DEVICE_GPU* d,
        const int& dim,
        const T* psi_L,
        const T* psi_R,
        const bool reduce = true);
};

// vector operator: result[i] = vector[i] / constant
template <typename T> struct vector_div_constant_op<T, psi::DEVICE_GPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const psi::DEVICE_GPU* d,
        const int dim,
        T* result,
        const T* vector,
        const Real constant);
};

// vector operator: result[i] = vector1[i](complex) * vector2[i](not complex)
template <typename T> struct vector_mul_vector_op<T, psi::DEVICE_GPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const psi::DEVICE_GPU* d,
        const int& dim,
        T* result,
        const T* vector1,
        const  Real* vector2);
};

// vector operator: result[i] = vector1[i](complex) / vector2[i](not complex)
template <typename T> struct vector_div_vector_op<T, psi::DEVICE_GPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const psi::DEVICE_GPU* d,
        const int& dim,
        T* result,
        const T* vector1,
        const Real* vector2);
};

// vector operator: result[i] = vector1[i] * constant1 + vector2[i] * constant2
template <typename FPTYPE> struct constantvector_addORsub_constantVector_op<FPTYPE, psi::DEVICE_GPU>
{
    void operator()(const psi::DEVICE_GPU* d,
                    const int& dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector1,
                    const FPTYPE constant1,
                    const std::complex<FPTYPE>* vector2,
                    const FPTYPE constant2);
};

template <typename T> struct matrixSetToAnother<T, psi::DEVICE_GPU>
{
    void operator()(const psi::DEVICE_GPU* d,
        const int& n,
        const T* A, // input
        const int& LDA,
        T* B, // output
        const int& LDB);
};


void createGpuBlasHandle();
void destoryBLAShandle();

#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
} // namespace hsolver

#endif // MODULE_HSOLVER_MATH_KERNEL_H