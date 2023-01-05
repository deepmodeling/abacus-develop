#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_psi/psi.h"
#include "module_psi/kernels/memory_op.h"

#include <hipblas.h>
#include <thrust/complex.h>
#include <hip/hip_runtime.h>

namespace hsolver {


static hipblasHandle_t cublas_handle = nullptr;

static inline
void xdot_wrapper(const int &n, const float * x, const int &incx, const float * y, const int &incy, float &result) {
    hipblasSdot(cublas_handle, n, x, incx, y, incy, &result);
}

static inline
void xdot_wrapper(const int &n, const double * x, const int &incx, const double * y, const int &incy, double &result) {
    hipblasDdot(cublas_handle, n, x, incx, y, incy, &result);
}

void createBLAShandle(){
    if (cublas_handle == nullptr) {
        hipblasCreate(&cublas_handle);
    }
}

void destoryBLAShandle(){
    if (cublas_handle != nullptr) {
        hipblasDestroy(cublas_handle);
        cublas_handle = nullptr;
    }
}

// Define the CUDA kernel:
template <typename FPTYPE>
__launch_bounds__(1024) 
__global__ void vector_div_constant_kernel(
    const int size, 
    thrust::complex<FPTYPE>* result, 
    const thrust::complex<FPTYPE>* vector, 
    const FPTYPE constant)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) 
    {
        result[i] = vector[i] / constant;
    }
}

template <typename FPTYPE>
__launch_bounds__(1024) 
__global__ void vector_mul_vector_kernel(
    const int size, 
    thrust::complex<FPTYPE>* result, 
    const thrust::complex<FPTYPE>* vector1, 
    const FPTYPE* vector2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) 
    {
        result[i] = vector1[i] * vector2[i];
    }
}

template <typename FPTYPE>
__launch_bounds__(1024) 
__global__ void vector_div_vector_kernel(
    const int size, 
    thrust::complex<FPTYPE>* result, 
    const thrust::complex<FPTYPE>* vector1, 
    const FPTYPE* vector2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) 
    {
        result[i] = vector1[i] / vector2[i];
    }
}

template <typename FPTYPE>
__launch_bounds__(1024) 
__global__ void constantvector_addORsub_constantVector_kernel(
    const int size,
    thrust::complex<FPTYPE>* result,
    const thrust::complex<FPTYPE>* vector1,
    const FPTYPE constant1,
    const thrust::complex<FPTYPE>* vector2,
    const FPTYPE constant2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) 
    {
        result[i] = vector1[i] * constant1 + vector2[i] * constant2;
    }
}

template <typename FPTYPE>
__launch_bounds__(1024) 
__global__ void matrix_transpose_kernel(
        const int row,
        const int col,
        const thrust::complex<FPTYPE>* in,
        thrust::complex<FPTYPE>* out)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < row)
    {
        for (int j = 0; j < col; j++)
        {
            out[j * row + i] = in[i * col + j];
        }
    }
}


template <typename FPTYPE>
__launch_bounds__(1024) 
__global__ void matrix_setTo_another_kernel(
        const int n,
        const int LDA,
        const int LDB,
        const thrust::complex<FPTYPE>* matrix_A,
        thrust::complex<FPTYPE>* matrix_B)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j < LDA && j < LDB)
    {
        for (int i = 0; i < n; i++)
        {
            matrix_B[i * LDB + j] = matrix_A[i * LDA + j];
        }
    }
}


// for this implementation, please check
// https://thrust.github.io/doc/group__transformed__reductions_ga321192d85c5f510e52300ae762c7e995.html denghui modify
// 2022-10-03 Note that ddot_(2*dim,a,1,           b,1) = REAL( zdotc_(dim,a,1,b,1) ) GPU specialization of actual computation.
template <typename FPTYPE>
FPTYPE zdot_real_op<FPTYPE, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& dim,
    const std::complex<FPTYPE>* psi_L,
    const std::complex<FPTYPE>* psi_R,
    const bool reduce)
{
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  // denghui modify 2022-10-07
  // Note that  ddot_(2*dim,a,1,b,1) = REAL( zdotc_(dim,a,1,b,1) )
  const FPTYPE* pL = reinterpret_cast<const FPTYPE*>(psi_L);
  const FPTYPE* pR = reinterpret_cast<const FPTYPE*>(psi_R);
  FPTYPE result = 0.0;
  xdot_wrapper(dim * 2, pL, 1, pR, 1, result);
  if (reduce) {
      Parallel_Reduce::reduce_double_pool(result);
  }
  return result;
}

// vector operator: result[i] = vector[i] / constant
template <typename FPTYPE>
void vector_div_constant_op<FPTYPE, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int dim,
    std::complex<FPTYPE>* result,
    const std::complex<FPTYPE>* vector,
    const FPTYPE constant)
{
    thrust::complex<FPTYPE>* result_tmp = reinterpret_cast<thrust::complex<FPTYPE>*>(result);
    const thrust::complex<FPTYPE>* vector_tmp = reinterpret_cast<const thrust::complex<FPTYPE>*>(vector);

    int thread = 1024;
    int block = (dim + thread - 1) / thread;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(vector_div_constant_kernel<FPTYPE>), dim3(block), dim3(thread), 0, 0, dim, result_tmp, vector_tmp, constant);
}

// vector operator: result[i] = vector1[i](complex) * vector2[i](not complex)
template <typename FPTYPE> 
void vector_mul_vector_op<FPTYPE, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& dim,
    std::complex<FPTYPE>* result,
    const std::complex<FPTYPE>* vector1,
    const FPTYPE* vector2)
{
    thrust::complex<FPTYPE>* result_tmp = reinterpret_cast<thrust::complex<FPTYPE>*>(result);
    const thrust::complex<FPTYPE>* vector1_tmp = reinterpret_cast<const thrust::complex<FPTYPE>*>(vector1);

    int thread = 1024;
    int block = (dim + thread - 1) / thread;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(vector_mul_vector_kernel<FPTYPE>), dim3(block), dim3(thread), 0, 0, dim, result_tmp, vector1_tmp, vector2);
}


// vector operator: result[i] = vector1[i](complex) / vector2[i](not complex)
template <typename FPTYPE> 
void vector_div_vector_op<FPTYPE, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& dim,
    std::complex<FPTYPE>* result,
    const std::complex<FPTYPE>* vector1,
    const FPTYPE* vector2)
{
    thrust::complex<FPTYPE>* result_tmp = reinterpret_cast<thrust::complex<FPTYPE>*>(result);
    const thrust::complex<FPTYPE>* vector1_tmp = reinterpret_cast<const thrust::complex<FPTYPE>*>(vector1);

    int thread = 1024;
    int block = (dim + thread - 1) / thread;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(vector_div_vector_kernel<FPTYPE>), dim3(block), dim3(thread), 0, 0, dim, result_tmp, vector1_tmp, vector2);
}

// vector operator: result[i] = vector1[i] * constant1 + vector2[i] * constant2
template <typename FPTYPE> 
void constantvector_addORsub_constantVector_op<FPTYPE, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& dim,
    std::complex<FPTYPE>* result,
    const std::complex<FPTYPE>* vector1,
    const FPTYPE constant1,
    const std::complex<FPTYPE>* vector2,
    const FPTYPE constant2)
{
    thrust::complex<FPTYPE>* result_tmp = reinterpret_cast<thrust::complex<FPTYPE>*>(result);
    const thrust::complex<FPTYPE>* vector1_tmp = reinterpret_cast<const thrust::complex<FPTYPE>*>(vector1);
    const thrust::complex<FPTYPE>* vector2_tmp = reinterpret_cast<const thrust::complex<FPTYPE>*>(vector2);

    int thread = 1024;
    int block = (dim + thread - 1) / thread;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(constantvector_addORsub_constantVector_kernel<FPTYPE>), dim3(block), dim3(thread), 0, 0, dim, result_tmp, vector1_tmp,constant1, vector2_tmp, constant2);
}

template <> 
void axpy_op<float, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& N,
    const std::complex<float> *alpha,
    const std::complex<float> *X,
    const int& incX,
    std::complex<float> *Y,
    const int& incY)
{
    hipblasCaxpy(cublas_handle, N, (hipblasComplex*)alpha, (hipblasComplex*)X, incX, (hipblasComplex*)Y, incY);
}

template <> 
void axpy_op<double, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& N,
    const std::complex<double> *alpha,
    const std::complex<double> *X,
    const int& incX,
    std::complex<double> *Y,
    const int& incY)
{
    hipblasZaxpy(cublas_handle, N, (hipblasDoubleComplex*)alpha, (hipblasDoubleComplex*)X, incX, (hipblasDoubleComplex*)Y, incY);
}

template <> 
void gemv_op<float, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const char& trans,
    const int& m,
    const int& n,
    const std::complex<float> *alpha,
    const std::complex<float> *A, 
    const int& lda,
    const std::complex<float> *X, 
    const int& incx,
    const std::complex<float> *beta,
    std::complex<float> *Y, 
    const int& incy)
{
    hipblasOperation_t cutrans = {};
    if (trans == 'N'){
        cutrans = HIPBLAS_OP_N;
    } 
    else if (trans == 'T'){
        cutrans = HIPBLAS_OP_T;
    } 
    else if (trans == 'C'){
        cutrans = HIPBLAS_OP_C;
    }
    hipblasCgemv(cublas_handle, cutrans, m, n, (hipblasComplex*)alpha, (hipblasComplex*)A, lda, (hipblasComplex*)X, incx, (hipblasComplex*)beta, (hipblasComplex*)Y, incx);
}

template <> 
void gemv_op<double, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const char& trans,
    const int& m,
    const int& n,
    const std::complex<double> *alpha,
    const std::complex<double> *A, 
    const int& lda,
    const std::complex<double> *X, 
    const int& incx,
    const std::complex<double> *beta,
    std::complex<double> *Y, 
    const int& incy)
{
    hipblasOperation_t cutrans = {};
    if (trans == 'N'){
        cutrans = HIPBLAS_OP_N;
    } 
    else if (trans == 'T'){
        cutrans = HIPBLAS_OP_T;
    } 
    else if (trans == 'C'){
        cutrans = HIPBLAS_OP_C;
    }
    hipblasZgemv(cublas_handle, cutrans, m, n, (hipblasDoubleComplex*)alpha, (hipblasDoubleComplex*)A, lda, (hipblasDoubleComplex*)X, incx, (hipblasDoubleComplex*)beta, (hipblasDoubleComplex*)Y, incx);
}

template <>
void scal_op<float, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
                                                 const int& N,
                                                 const std::complex<float>* alpha,
                                                 std::complex<float>* X,
                                                 const int& incx)
{
    hipblasCscal(cublas_handle, N, (hipblasComplex*)alpha, (hipblasComplex*)X, incx);
}

template <>
void scal_op<double, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
                                                  const int& N,
                                                  const std::complex<double>* alpha,
                                                  std::complex<double>* X,
                                                  const int& incx)
{
    hipblasZscal(cublas_handle, N, (hipblasDoubleComplex*)alpha, (hipblasDoubleComplex*)X, incx);
}


template <>
void gemm_op<float, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
                                                 const char& transa, 
                                                 const char& transb, 
                                                 const int& m, 
                                                 const int& n, 
                                                 const int& k,
                                                 const std::complex<float> *alpha, 
                                                 const std::complex<float> *a, 
                                                 const int& lda, 
                                                 const std::complex<float> *b, 
                                                 const int& ldb,
                                                 const std::complex<float> *beta, 
                                                 std::complex<float> *c, 
                                                 const int& ldc)
{
    hipblasOperation_t cutransA = {};
    hipblasOperation_t cutransB = {};
    // cutransA
    if (transa == 'N'){
        cutransA = HIPBLAS_OP_N;
    } 
    else if (transa == 'T'){
        cutransA = HIPBLAS_OP_T;
    } 
    else if (transa == 'C'){
        cutransA = HIPBLAS_OP_C;
    } 
    // cutransB
    if (transb == 'N'){
        cutransB = HIPBLAS_OP_N;
    } 
    else if (transb == 'T'){
        cutransB = HIPBLAS_OP_T;
    } 
    else if (transb == 'C'){
        cutransB = HIPBLAS_OP_C;
    }
    hipblasCgemm(cublas_handle, cutransA, cutransB, m, n ,k, (hipblasComplex*)alpha, (hipblasComplex*)a , lda, (hipblasComplex*)b, ldb, (hipblasComplex*)beta, (hipblasComplex*)c, ldc);
}

template <>
void gemm_op<double, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
                                                 const char& transa, 
                                                 const char& transb, 
                                                 const int& m, 
                                                 const int& n, 
                                                 const int& k,
                                                 const std::complex<double> *alpha, 
                                                 const std::complex<double> *a, 
                                                 const int& lda, 
                                                 const std::complex<double> *b, 
                                                 const int& ldb,
                                                 const std::complex<double> *beta, 
                                                 std::complex<double> *c, 
                                                 const int& ldc)
{
    hipblasOperation_t cutransA;
    hipblasOperation_t cutransB;
    // cutransA
    if (transa == 'N'){
        cutransA = HIPBLAS_OP_N;
    } 
    else if (transa == 'T'){
        cutransA = HIPBLAS_OP_T;
    } 
    else if (transa == 'C'){
        cutransA = HIPBLAS_OP_C;
    } 
    // cutransB
    if (transb == 'N'){
        cutransB = HIPBLAS_OP_N;
    } 
    else if (transb == 'T'){
        cutransB = HIPBLAS_OP_T;
    } 
    else if (transb == 'C'){
        cutransB = HIPBLAS_OP_C;
    }
    hipblasZgemm(cublas_handle, cutransA, cutransB, m, n ,k, (hipblasDoubleComplex*)alpha, (hipblasDoubleComplex*)a , lda, (hipblasDoubleComplex*)b, ldb, (hipblasDoubleComplex*)beta, (hipblasDoubleComplex*)c, ldc);
}

template <>
void matrixTranspose_op<float, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
                                                             const int& row,
                                                             const int& col,
                                                             const std::complex<float>* input_matrix,
                                                             std::complex<float>* output_matrix)
{
    std::complex<float>* device_temp = nullptr;
    psi::memory::resize_memory_op<std::complex<float>, psi::DEVICE_GPU>()(d, device_temp, row * col);

    if (row == col)
    {
        float2 ONE, ZERO;
        ONE.x = 1.0;
        ONE.y = 0.0;
        ZERO.x = ZERO.y = 0.0;

        // use 'geam' API todo transpose.
        hipblasCgeam(cublas_handle, HIPBLAS_OP_T, HIPBLAS_OP_N, col, row,
                                   reinterpret_cast<const hipblasComplex *>(&ONE), (hipblasComplex*)input_matrix, col,
                                   reinterpret_cast<const hipblasComplex *>(&ZERO), (hipblasComplex*)input_matrix, col, (hipblasComplex*)device_temp, col);
    } else
    {
        int thread = 1024;
        int block = (row + col + thread - 1) / thread;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(matrix_transpose_kernel<float>), dim3(block), dim3(thread), 0, 0, row, col, (thrust::complex<float>*)input_matrix, (thrust::complex<float>*)device_temp);
    }

    psi::memory::synchronize_memory_op<std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_GPU>()(d, d, output_matrix, device_temp, row * col);

    psi::memory::delete_memory_op<std::complex<float>, psi::DEVICE_GPU>()(d, device_temp);

}

template <>
void matrixTranspose_op<double, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
                                                             const int& row,
                                                             const int& col,
                                                             const std::complex<double>* input_matrix,
                                                             std::complex<double>* output_matrix)
{
    std::complex<double>* device_temp = nullptr;
    psi::memory::resize_memory_op<std::complex<double>, psi::DEVICE_GPU>()(d, device_temp, row * col);

    if (row == col)
    {
        hipblasDoubleComplex ONE{1.0, 0.0}, ZERO{0.0, 0.0};
        // use 'geam' API todo transpose.
        hipblasZgeam(cublas_handle, HIPBLAS_OP_T, HIPBLAS_OP_N, col, row, &ONE, (hipblasDoubleComplex*)input_matrix, col, &ZERO, (hipblasDoubleComplex*)input_matrix, col, (hipblasDoubleComplex*)device_temp, col);
    } else
    {
        int thread = 1024;
        int block = (row + col + thread - 1) / thread;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(matrix_transpose_kernel<double>), dim3(block), dim3(thread), 0, 0, row, col, (thrust::complex<double>*)input_matrix, (thrust::complex<double>*)device_temp);
    }
    
    psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_GPU>()(d, d, output_matrix, device_temp, row * col);

    psi::memory::delete_memory_op<std::complex<double>, psi::DEVICE_GPU>()(d, device_temp);
    
}


template <typename FPTYPE> 
void matrixSetToAnother<FPTYPE, psi::DEVICE_GPU>::operator()(
            const psi::DEVICE_GPU* d,
            const int& n,
            const std::complex<FPTYPE>* A,
            const int& LDA,
            std::complex<FPTYPE>* B,
            const int& LDB)
{
    int thread = 1024;
    int block = (LDA + thread - 1) / thread;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(matrix_setTo_another_kernel<FPTYPE>), dim3(block), dim3(thread), 0, 0, n, LDA, LDB, reinterpret_cast<const thrust::complex<FPTYPE>*>(A), reinterpret_cast<thrust::complex<FPTYPE>*>(B));
}



// Explicitly instantiate functors for the types of functor registered.
template struct zdot_real_op<float, psi::DEVICE_GPU>;
template struct vector_div_constant_op<float, psi::DEVICE_GPU>;
template struct vector_mul_vector_op<float, psi::DEVICE_GPU>;
template struct vector_div_vector_op<float, psi::DEVICE_GPU>;
template struct constantvector_addORsub_constantVector_op<float, psi::DEVICE_GPU>;
template struct matrixSetToAnother<float, psi::DEVICE_GPU>;

template struct zdot_real_op<double, psi::DEVICE_GPU>;
template struct vector_div_constant_op<double, psi::DEVICE_GPU>;
template struct vector_mul_vector_op<double, psi::DEVICE_GPU>;
template struct vector_div_vector_op<double, psi::DEVICE_GPU>;
template struct constantvector_addORsub_constantVector_op<double, psi::DEVICE_GPU>;
template struct matrixSetToAnother<double, psi::DEVICE_GPU>;

}  // namespace hsolver
