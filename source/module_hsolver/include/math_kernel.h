// TODO: This is a temperary location for these functions.
// And will be moved to a global module(module base) later.
#ifndef MODULE_HSOLVER_MATH_KERNEL_H
#define MODULE_HSOLVER_MATH_KERNEL_H

#include "module_base/blas_connector.h"
#include "module_psi/psi.h"
#include "src_parallel/parallel_reduce.h"

namespace hsolver
{

template <typename FPTYPE, typename Device> 
struct zdot_real_op {
  FPTYPE operator() (
      const Device* d,
      const int& dim,
      const std::complex<FPTYPE>* psi_L,
      const std::complex<FPTYPE>* psi_R,
      const bool reduce = true);
};

// vector operator: result[i] = vector[i] / constant
template <typename FPTYPE, typename Device> struct vector_div_constant_op
{
    void operator()(const Device* d,
                    const int dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector,
                    const FPTYPE constant);
};

// replace vector_div_constant_op : x = alpha * x
template <typename FPTYPE, typename Device> struct scal_op
{
    void operator()(const Device* d,
                    const int& N,
                    const std::complex<FPTYPE>* alpha,
                    std::complex<FPTYPE>* X,
                    const int& incx);
};

// vector operator: result[i] = vector1[i](complex) * vector2[i](not complex)
template <typename FPTYPE, typename Device> struct vector_mul_vector_op
{
    void operator()(const Device* d,
                    const int& dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector1,
                    const FPTYPE* vector2);
};

// vector operator: result[i] = vector1[i](complex) / vector2[i](not complex)
template <typename FPTYPE, typename Device> struct vector_div_vector_op
{
    void operator()(const Device* d,
                    const int& dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector1,
                    const FPTYPE* vector2);
};

// vector operator: result[i] = vector1[i] * constant1 + vector2[i] * constant2
template <typename FPTYPE, typename Device> struct constantvector_addORsub_constantVector_op
{
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
    void operator()(const Device* d,
                    const int& N,
                    const std::complex<FPTYPE>* alpha,
                    const std::complex<FPTYPE>* X,
                    const int& incX,
                    std::complex<FPTYPE>* Y,
                    const int& incY);
};

// compute y = alpha * op(A) * x + beta * y
template <typename FPTYPE, typename Device> struct gemv_op
{
    void operator()(const Device* d,
                    const char& trans,
                    const int& m,
                    const int& n,
                    const std::complex<FPTYPE>* alpha,
                    const std::complex<FPTYPE>* A,
                    const int& lda,
                    const std::complex<FPTYPE>* X,
                    const int& incx,
                    const std::complex<FPTYPE>* beta,
                    std::complex<FPTYPE>* Y,
                    const int& incy);
};


// compute C = alpha * op(A) * op(B) + beta * C
template <typename FPTYPE, typename Device> struct gemm_op
{
    void operator()(const Device* d,
                    const char& transa, 
                    const char& transb, 
                    const int& m, 
                    const int& n, 
                    const int& k,
		            const std::complex<FPTYPE> *alpha, 
                    const std::complex<FPTYPE> *a, 
                    const int& lda, 
                    const std::complex<FPTYPE> *b, 
                    const int& ldb,
		            const std::complex<FPTYPE> *beta, 
                    std::complex<FPTYPE> *c, 
                    const int& ldc);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM

// Partially specialize functor for psi::GpuDevice.
template <typename FPTYPE> 
struct zdot_real_op<FPTYPE, psi::DEVICE_GPU> {
  FPTYPE operator()(
    const psi::DEVICE_GPU* d,
    const int& dim,
    const std::complex<FPTYPE>* psi_L,
    const std::complex<FPTYPE>* psi_R,
    const bool reduce = true);
};

// vector operator: result[i] = vector[i] / constant
template <typename FPTYPE> struct vector_div_constant_op<FPTYPE, psi::DEVICE_GPU>
{
    void operator()(const psi::DEVICE_GPU* d,
                    const int dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector,
                    const FPTYPE constant);
};

// vector operator: result[i] = vector1[i](complex) * vector2[i](not complex)
template <typename FPTYPE> struct vector_mul_vector_op<FPTYPE, psi::DEVICE_GPU>
{
    void operator()(const psi::DEVICE_GPU* d,
                    const int& dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector1,
                    const FPTYPE* vector2);
};

// vector operator: result[i] = vector1[i](complex) / vector2[i](not complex)
template <typename FPTYPE> struct vector_div_vector_op<FPTYPE, psi::DEVICE_GPU>
{
    void operator()(const psi::DEVICE_GPU* d,
                    const int& dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector1,
                    const FPTYPE* vector2);
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

void createBLAShandle();
void destoryBLAShandle();

#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
} // namespace hsolver

#endif // MODULE_HSOLVER_MATH_KERNEL_H