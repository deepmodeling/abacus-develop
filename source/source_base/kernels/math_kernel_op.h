// TODO: This is a temperary location for these functions.
// And will be moved to a global module(module base) later.
#ifndef MODULE_HSOLVER_MATH_KERNEL_H
#define MODULE_HSOLVER_MATH_KERNEL_H

#include "source_base/macros.h"

#include "source_base/module_external/blas_connector.h"
#include "source_base/parallel_reduce.h"

#include "source_base/module_device/memory_op.h"
#include "source_base/module_device/types.h"

#if defined(__CUDA) || defined(__UT_USE_CUDA)
#include <cuda_runtime.h>

#include "cublas_v2.h"
#endif //__CUDA || __UT_USE_CUDA

namespace ModuleBase {

//---------------------------------------------------------------------------------
//-----------------------------0. Tool Functions-----------------------------------
//---------------------------------------------------------------------------------
inline std::complex<double> set_real_tocomplex(const std::complex<double> &x) {
  return {x.real(), 0.0};
}

inline std::complex<float> set_real_tocomplex(const std::complex<float> &x) {
  return {x.real(), 0.0};
}

inline double set_real_tocomplex(const double &x) { return x; }

inline float set_real_tocomplex(const float &x) { return x; }

inline std::complex<double> get_conj(const std::complex<double> &x) {
  return {x.real(), -x.imag()};
}

inline std::complex<float> get_conj(const std::complex<float> &x) {
  return {x.real(), -x.imag()};
}

inline double get_conj(const double &x) { return x; }

inline float get_conj(const float &x) { return x; }


//---------------------------------------------------------------------------------
//-----------------------------1. Vector Operations--------------------------------
//---------------------------------------------------------------------------------
template <typename FPTYPE, typename Device> struct scal_op {
  /// @brief x = alpha * x, where alpha and x are complex numbers
  ///
  /// Input Parameters
  /// \param N : array size
  /// \param alpha : input constant
  /// \param X : input array
  /// \param incx : computing strip of array X
  ///
  /// Output Parameters
  /// \param X : output array
  void operator()(const int &N,
                  const std::complex<FPTYPE> *alpha, std::complex<FPTYPE> *X,
                  const int &incx);
};

template <typename T, typename Device> struct vector_mul_real_op {
  using Real = typename GetTypeReal<T>::type;
  /// @brief result[i] = vector[i] * constant, where vector is complex number and constant is real number。
  ///        It is different from the scal_op, which is used to multiply a complex number by a complex number.
  ///
  /// Input Parameters
  /// \param dim : array size
  /// \param vector : input array
  /// \param constant : input constant
  ///
  /// Output Parameters
  /// \param result : output array
  /// \note Use mulitple instead of divide. It is faster.
  void operator()(const int dim, T* result, const T* vector, const Real constant);
};

// vector operator: result[i] = vector1[i](complex) * vector2[i](not complex)
template <typename T, typename Device> struct vector_mul_vector_op {
  using Real = typename GetTypeReal<T>::type;
  /// @brief result[i] = vector1[i](complex) * vector2[i](not complex)
  ///
  /// Input Parameters
  /// \param dim : array size
  /// \param vector1 : input array A
  /// \param vector2 : input array B
  /// \param add : flag to control whether to add the result to the output array
  ///
  /// Output Parameters
  /// \param result : output array
  void operator()(const int& dim, T* result, const T* vector1, const Real* vector2, const bool& add = false);
};

// vector operator: result[i] = vector[i] / constant
template <typename T, typename Device> struct vector_div_constant_op {
  using Real = typename GetTypeReal<T>::type;
  /// @brief result[i] = vector[i] / constant
  ///
  /// Input Parameters
  /// \param dim : array size
  /// \param vector : input array 
  /// \param constant : input constant
  ///
  /// Output Parameters
  /// \param result : output array
  void operator()(const int& dim, T* result, const T* vector, const Real constant);
};

// vector operator: result[i] = vector1[i](complex) / vector2[i](not complex)
template <typename T, typename Device> struct vector_div_vector_op {
  using Real = typename GetTypeReal<T>::type;
  /// @brief result[i] = vector1[i](complex) / vector2[i](not complex)
  ///
  /// Input Parameters
  /// \param dim : array size
  /// \param vector1 : input array A
  /// \param vector2 : input array B
  ///
  /// Output Parameters
  /// \param result : output array
  void operator()(const int &dim, T *result, const T *vector1,
                  const Real *vector2);
};

//  compute Y = alpha * X + Y
template <typename T, typename Device> struct axpy_op {
  /// @brief Y = alpha * X + Y
  ///
  /// Input Parameters
  /// \param N : array size
  /// \param alpha : input constant alpha
  /// \param X : input array X
  /// \param incX : computing strip of X
  /// \param Y : computing strip of Y
  /// \param incY : computing strip of Y
  ///
  /// Output Parameters
  /// \param Y : output array Y
  void operator()(const int &N, const T *alpha, const T *X,
                  const int &incX, T *Y, const int &incY);
};

// vector operator: result[i] = vector1[i] * constant1 + vector2[i] * constant2
template <typename T, typename Device>
struct vector_add_vector_op {
  using Real = typename GetTypeReal<T>::type;
  /// @brief result[i] = vector1[i] * constant1 + vector2[i] * constant2
  ///
  /// Input Parameters
  /// \param dim : array size
  /// \param vector1 : input array A
  /// \param constant1 : input constant a
  /// \param vector2 : input array B
  /// \param constant2 : input constant b
  ///
  /// Output Parameters
  /// \param result : output array
  void operator()(const int &dim, T *result, const T *vector1,
                  const Real constant1, const T *vector2, const Real constant2);
};

template <typename T, typename Device> struct dot_real_op {
  using Real = typename GetTypeReal<T>::type;
  /// @brief dot_real_op computes the dot product of the given complex
  /// arrays(treated as float arrays). And there's may have MPI communications
  /// while enabling planewave parallization strategy.
  ///
  /// Input Parameters
  /// \param dim : array size
  /// \param psi_L : input array A
  /// \param psi_R : input array B
  /// \param reduce : flag to control whether to perform the MPI communications
  ///
  /// \return
  /// FPTYPE : dot product result
  Real operator()(const int &dim, const T *psi_L,
                  const T *psi_R, const bool reduce = true);
};


//---------------------------------------------------------------------------------
//-----------------------------2. Matrix Operations--------------------------------
//---------------------------------------------------------------------------------

// compute y = alpha * op(A) * x + beta * y
template <typename T, typename Device> struct gemv_op {
  /// @brief y = alpha * op(A) * x + beta * y
  ///
  /// Input Parameters
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
  void operator()(const char &trans, const int &m,
                  const int &n, const T *alpha, const T *A, const int &lda,
                  const T *X, const int &incx, const T *beta, T *Y,
                  const int &incy);
};

// compute C = alpha * op(A) * op(B) + beta * C
template <typename T, typename Device> struct gemm_op {
  /// @brief C = alpha * op(A) * op(B) + beta * C
  ///
  /// Input Parameters
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
  void operator()(const char &transa, const char &transb,
                  const int &m, const int &n, const int &k, const T *alpha,
                  const T *a, const int &lda, const T *b, const int &ldb,
                  const T *beta, T *c, const int &ldc);
};

#ifdef __DSP
// compute C = alpha * op(A) * op(B) + beta * C on DSP Hardware
template <typename T, typename Device> struct gemm_op_mt {
  /// @brief C = alpha * op(A) * op(B) + beta * C
  ///
  /// Input Parameters
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
  void operator()(const char &transa, const char &transb,
                  const int &m, const int &n, const int &k, const T *alpha,
                  const T *a, const int &lda, const T *b, const int &ldb,
                  const T *beta, T *c, const int &ldc);
};
#endif

template <typename T, typename Device> struct matrixTranspose_op {
  /// @brief transpose the input matrix
  ///
  /// Input Parameters
  /// \param row : first dimension of matrix
  /// \param col : second dimension of matrix
  /// \param input_matrix : input matrix
  ///
  /// Output Parameters
  /// \param output_matrix : output matrix
  void operator()(const int &row, const int &col,
                  const T *input_matrix, T *output_matrix);
};

template <typename T, typename Device> struct matrixCopy {
  /// @brief copy matrix A to B, they can have different leading dimensions
  ///
  /// Input Parameters
  /// \param n1 : first dimension of matrix
  /// \param n2 : second dimension of matrix
  /// \param A : input matrix A
  /// \param LDA : leading dimension of A
  /// \param LDB : leading dimension of B
  ///
  /// Output Parameters
  /// \param B : output matrix B
  void operator()(const int& n1, const int& n2, const T* A, const int& LDA, T* B, const int& LDB);
};

template <typename T, typename Device>
struct apply_eigenvalues_op {
    using Real = typename GetTypeReal<T>::type;

    void operator()(const Device *d, const int &nbase, const int &nbase_x, const int &notconv,
                    T *result, const T *vectors, const Real *eigenvalues);
};

template <typename T, typename Device>
struct precondition_op {
    using Real = typename GetTypeReal<T>::type;
    void operator()(const Device* d,
                   const int& dim,
                   T* psi_iter,
                   const int& nbase,
                   const int& notconv,
                   const Real* precondition,  
                   const Real* eigenvalues);
};

template <typename T, typename Device>
struct normalize_op {
    using Real = typename GetTypeReal<T>::type;
    void operator()(const Device* d,
                   const int& dim,
                   T* psi_iter,
                   const int& nbase,
                   const int& notconv,
                   Real* psi_norm = nullptr);
};

template <typename T>
struct normalize_op<T, base_device::DEVICE_GPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(const base_device::DEVICE_GPU* d,
                   const int& dim,
                   T* psi_iter,
                   const int& nbase,
                   const int& notconv,
                   Real* psi_norm);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
// Partially specialize functor for base_device::GpuDevice.
template <typename T> struct dot_real_op<T, base_device::DEVICE_GPU> {
  using Real = typename GetTypeReal<T>::type;
  Real operator()(const int &dim,
                  const T *psi_L, const T *psi_R, const bool reduce = true);
};

// vector operator: result[i] = vector[i] / constant
template <typename T>
struct vector_mul_real_op<T, base_device::DEVICE_GPU>
{
  using Real = typename GetTypeReal<T>::type;
  void operator()(const int dim, T* result, const T* vector, const Real constant);
};

// vector operator: result[i] = vector1[i](complex) * vector2[i](not complex)
template <typename T> struct vector_mul_vector_op<T, base_device::DEVICE_GPU> {
  using Real = typename GetTypeReal<T>::type;
  void operator()(const int& dim, T* result, const T* vector1, const Real* vector2, const bool& add = false);
};

// vector operator: result[i] = vector[i] / constant
template <typename T> struct vector_div_constant_op<T, base_device::DEVICE_GPU> {
  using Real = typename GetTypeReal<T>::type;
  void operator()(const int& dim, T* result, const T* vector, const Real constant);
};

// vector operator: result[i] = vector1[i](complex) / vector2[i](not complex)
template <typename T> struct vector_div_vector_op<T, base_device::DEVICE_GPU> {
  using Real = typename GetTypeReal<T>::type;
  void operator()(const int &dim, T *result,
                  const T *vector1, const Real *vector2);
};

// vector operator: result[i] = vector1[i] * constant1 + vector2[i] * constant2
template <typename T>
struct vector_add_vector_op<T, base_device::DEVICE_GPU> {
  using Real = typename GetTypeReal<T>::type;
  void operator()(const int &dim, T *result,
                  const T *vector1, const Real constant1, const T *vector2,
                  const Real constant2);
};

template <typename T> struct matrixCopy<T, base_device::DEVICE_GPU> {
    void operator()(const int& n1,
                    const int& n2,
                    const T* A, // input
                    const int& LDA,
                    T* B, // output
                    const int& LDB);
};

void createGpuBlasHandle();
void destoryBLAShandle();

// vector operator: result[i] = -lambda[i] * vector[i]
template <typename T> struct apply_eigenvalues_op<T, base_device::DEVICE_GPU> {
    using Real = typename GetTypeReal<T>::type;

    void operator()(const base_device::DEVICE_GPU *d, const int &nbase, const int &nbase_x, const int &notconv,
                    T *result, const T *vectors, const Real *eigenvalues);
};

template <typename T>
struct precondition_op<T, base_device::DEVICE_GPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(const base_device::DEVICE_GPU* d,
                   const int& dim,
                   T* psi_iter,
                   const int& nbase,
                   const int& notconv,
                   const Real* precondition,
                   const Real* eigenvalues);
};

#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
} // namespace hsolver

#endif // MODULE_HSOLVER_MATH_KERNEL_H