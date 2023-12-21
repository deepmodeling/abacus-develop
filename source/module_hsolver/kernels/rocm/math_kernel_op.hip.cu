#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_psi/psi.h"
#include "module_psi/kernels/memory_op.h"

#include <thrust/complex.h>
#include <hipblas/hipblas.h>
#include <hip/hip_runtime.h>

#define WARP_SIZE 32
#define FULL_MASK 0xffffffff
#define THREAD_PER_BLOCK 256
template <>
struct GetTypeReal<thrust::complex<float>> {
    using type = float; /**< The return type specialization for std::complex<double>. */
};
template <>
struct GetTypeReal<thrust::complex<double>> {
    using type = double; /**< The return type specialization for std::complex<double>. */
};

namespace hsolver {

template <typename T>
struct GetTypeThrust {
    using type = T;
};

template <>
struct GetTypeThrust<std::complex<float>> {
    using type = thrust::complex<float>; /**< The return type specialization for std::complex<float>. */
};

template <>
struct GetTypeThrust<std::complex<double>> {
    using type = thrust::complex<double>; /**< The return type specialization for std::complex<float>. */
};

static hipblasHandle_t cublas_handle = nullptr;

static inline
void xdot_wrapper(const int &n, const float * x, const int &incx, const float * y, const int &incy, float &result) {
    hipblasSdot(cublas_handle, n, x, incx, y, incy, &result);
}

static inline
void xdot_wrapper(const int &n, const double * x, const int &incx, const double * y, const int &incy, double &result) {
    hipblasDdot(cublas_handle, n, x, incx, y, incy, &result);
}

void createGpuBlasHandle(){
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

template <typename Real>
__global__ void line_minimize_with_block(
        thrust::complex<Real>* grad,
        thrust::complex<Real>* hgrad,
        thrust::complex<Real>* psi,
        thrust::complex<Real>* hpsi,
        const int n_basis,
        const int n_basis_max)
{
    int band_idx = blockIdx.x; // band_idx
    int tid = threadIdx.x; // basis_idx
    int item = 0;
    Real epsilo_0 = 0.0, epsilo_1 = 0.0, epsilo_2 = 0.0;
    Real theta = 0.0, cos_theta = 0.0, sin_theta = 0.0;
    __shared__ Real data[THREAD_PER_BLOCK * 3];

    data[tid] = 0;

    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += THREAD_PER_BLOCK) {
        item = band_idx * n_basis_max + basis_idx;
        data[tid] += (grad[item] * thrust::conj(grad[item])).real();
    }
    __syncthreads();
    // just do some parallel reduction in shared memory
    for (int ii = THREAD_PER_BLOCK >> 1; ii > 0; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
        }
        __syncthreads();
    }

    Real norm = 1.0 / sqrt(data[0]);
    __syncthreads();

    data[tid] = 0;
    data[THREAD_PER_BLOCK + tid] = 0;
    data[2 * THREAD_PER_BLOCK + tid] = 0;
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += THREAD_PER_BLOCK) {
        item = band_idx * n_basis_max + basis_idx;
        grad[item] *= norm;
        hgrad[item] *= norm;
        data[tid] += (hpsi[item] * thrust::conj(psi[item])).real();
        data[THREAD_PER_BLOCK + tid] += (grad[item] * thrust::conj(hpsi[item])).real();
        data[2 * THREAD_PER_BLOCK + tid] += (grad[item] * thrust::conj(hgrad[item])).real();
    }
    __syncthreads();

    // just do some parallel reduction in shared memory
    for (int ii = THREAD_PER_BLOCK >> 1; ii > 0; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
            data[THREAD_PER_BLOCK + tid] += data[THREAD_PER_BLOCK + tid + ii];
            data[2 * THREAD_PER_BLOCK + tid] += data[2 * THREAD_PER_BLOCK + tid + ii];
        }
        __syncthreads();
    }
    epsilo_0 = data[0];
    epsilo_1 = data[THREAD_PER_BLOCK];
    epsilo_2 = data[2 * THREAD_PER_BLOCK];

    theta = 0.5 * abs(atan(2 * epsilo_1/(epsilo_0 - epsilo_2)));
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += THREAD_PER_BLOCK) {
        item = band_idx * n_basis_max + basis_idx;
        psi [item] = psi [item] * cos_theta + grad [item] * sin_theta;
        hpsi[item] = hpsi[item] * cos_theta + hgrad[item] * sin_theta;
    }
}

template <typename Real>
__global__ void calc_grad_with_block(
        const Real* prec,
        Real* err,
        Real* beta,
        thrust::complex<Real>* psi,
        thrust::complex<Real>* hpsi,
        thrust::complex<Real>* grad,
        thrust::complex<Real>* grad_old,
        const int n_basis,
        const int n_basis_max)
{
    int band_idx = blockIdx.x; // band_idx
    int tid = threadIdx.x; // basis_idx
    int item = 0;
    Real err_st = 0.0;
    Real beta_st = 0.0;
    Real epsilo = 0.0;
    Real grad_2 = 0.0;
    thrust::complex<Real> grad_1 = {0, 0};
    __shared__ Real data[THREAD_PER_BLOCK * 2];

    // Init shared memory
    data[tid] = 0;

    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += THREAD_PER_BLOCK) {
        item = band_idx * n_basis_max + basis_idx;
        data[tid] += (psi[item] * thrust::conj(psi[item])).real();
    }
    __syncthreads();
    // just do some parallel reduction in shared memory
    for (int ii = THREAD_PER_BLOCK >> 1; ii > 0; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
        }
        __syncthreads();
    }

    Real norm = 1.0 / sqrt(data[0]);
    __syncthreads();

    data[tid] = 0;
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += THREAD_PER_BLOCK) {
        item = band_idx * n_basis_max + basis_idx;
        psi[item] *= norm;
        hpsi[item] *= norm;
        data[tid] += (hpsi[item] * thrust::conj(psi[item])).real();
    }
    __syncthreads();

    // just do some parallel reduction in shared memory
    for (int ii = THREAD_PER_BLOCK >> 1; ii > 0; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
        }
        __syncthreads();
    }
    epsilo = data[0];
    __syncthreads();

    data[tid] = 0;
    data[THREAD_PER_BLOCK + tid] = 0;
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += THREAD_PER_BLOCK) {
        item = band_idx * n_basis_max + basis_idx;
        grad_1 = hpsi[item] - epsilo * psi[item];
        grad_2 = thrust::norm(grad_1);
        data[tid] += grad_2;
        data[THREAD_PER_BLOCK + tid] += grad_2 / prec[basis_idx];
    }
    __syncthreads();

    // just do some parallel reduction in shared memory
    for (int ii = THREAD_PER_BLOCK >> 1; ii > 0; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
            data[THREAD_PER_BLOCK + tid] += data[THREAD_PER_BLOCK + tid + ii];
        }
        __syncthreads();
    }
    err_st = data[0];
    beta_st = data[THREAD_PER_BLOCK];
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += THREAD_PER_BLOCK) {
        item = band_idx * n_basis_max + basis_idx;
        grad_1 = hpsi[item] - epsilo * psi[item];
        grad[item] = -grad_1 / prec[basis_idx] + beta_st / beta[band_idx] * grad_old[item];
    }

    __syncthreads();
    if (tid == 0) {
        beta[band_idx] = beta_st;
        err[band_idx] = sqrt(err_st);
    }
}

// Define the CUDA kernel:
template <typename T>
__launch_bounds__(1024) 
__global__ void vector_div_constant_kernel(
    const int size, 
    T* result,
    const T* vector,
    const typename GetTypeReal<T>::type constant)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) 
    {
        result[i] = vector[i] / constant;
    }
}

template <typename T>
__launch_bounds__(1024) 
__global__ void vector_mul_vector_kernel(
    const int size, 
    T* result,
    const T* vector1,
    const typename GetTypeReal<T>::type* vector2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) 
    {
        result[i] = vector1[i] * vector2[i];
    }
}

template <typename T>
__launch_bounds__(1024) 
__global__ void vector_div_vector_kernel(
    const int size, 
    T* result,
    const T* vector1,
    const typename GetTypeReal<T>::type* vector2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) 
    {
        result[i] = vector1[i] / vector2[i];
    }
}

template <typename T, typename Real>
__launch_bounds__(1024) 
__global__ void constantvector_addORsub_constantVector_kernel(
    const int size,
    T* result,
    const T* vector1,
    const Real constant1,
    const T* vector2,
    const Real constant2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) 
    {
        result[i] = vector1[i] * constant1 + vector2[i] * constant2;
    }
}

template <typename T>
__launch_bounds__(1024) 
__global__ void matrix_transpose_kernel(
        const int row,
        const int col,
    const T* in,
    T* out)
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


template <typename T>
__launch_bounds__(1024) 
__global__ void matrix_setTo_another_kernel(
        const int n,
        const int LDA,
        const int LDB,
    const T* matrix_A,
    T* matrix_B)
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

template <typename T>
void line_minimize_with_block_op<T, psi::DEVICE_GPU>::operator()(
        T* grad_out,
        T* hgrad_out,
        T* psi_out,
        T* hpsi_out,
        const int &n_basis,
        const int &n_basis_max,
        const int &n_band)
{
    auto A = reinterpret_cast<thrust::complex<Real>*>(grad_out);
    auto B = reinterpret_cast<thrust::complex<Real>*>(hgrad_out);
    auto C = reinterpret_cast<thrust::complex<Real>*>(psi_out);
    auto D = reinterpret_cast<thrust::complex<Real>*>(hpsi_out);

    line_minimize_with_block<Real><<<n_band, THREAD_PER_BLOCK>>>(
            A, B, C, D,
            n_basis, n_basis_max);
}

template <typename T>
void calc_grad_with_block_op<T, psi::DEVICE_GPU>::operator()(
        const Real* prec_in,
        Real* err_out,
        Real* beta_out,
        T* psi_out,
        T* hpsi_out,
        T* grad_out,
        T* grad_old_out,
        const int &n_basis,
        const int &n_basis_max,
        const int &n_band)
{
    auto A = reinterpret_cast<thrust::complex<Real>*>(psi_out);
    auto B = reinterpret_cast<thrust::complex<Real>*>(hpsi_out);
    auto C = reinterpret_cast<thrust::complex<Real>*>(grad_out);
    auto D = reinterpret_cast<thrust::complex<Real>*>(grad_old_out);

    calc_grad_with_block<Real><<<n_band, THREAD_PER_BLOCK>>>(
            prec_in, err_out, beta_out,
            A, B, C, D,
            n_basis, n_basis_max);
}

template <>
double dot_real_op<double, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& dim,
    const double* psi_L,
    const double* psi_R,
    const bool reduce)
{
    double result = 0.0;
    xdot_wrapper(dim, psi_L, 1, psi_R, 1, result);
    if (reduce) {
        Parallel_Reduce::reduce_pool(result);
    }
    return result;
}

// for this implementation, please check
// https://thrust.github.io/doc/group__transformed__reductions_ga321192d85c5f510e52300ae762c7e995.html denghui modify
// 2022-10-03 Note that ddot_(2*dim,a,1,           b,1) = REAL( zdotc_(dim,a,1,b,1) ) GPU specialization of actual computation.
template <typename FPTYPE>
inline FPTYPE dot_complex_wrapper(
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
        Parallel_Reduce::reduce_pool(result);
    }
    return result;
}
template <>
float dot_real_op<std::complex<float>, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& dim,
    const std::complex<float>* psi_L,
    const std::complex<float>* psi_R,
    const bool reduce)
{
    return dot_complex_wrapper(d, dim, psi_L, psi_R, reduce);
}
template <>
double dot_real_op<std::complex<double>, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& dim,
    const std::complex<double>* psi_L,
    const std::complex<double>* psi_R,
    const bool reduce)
{
    return dot_complex_wrapper(d, dim, psi_L, psi_R, reduce);
}

template <>
void vector_div_constant_op<double, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int dim,
    double* result,
    const double* vector,
    const double constant)
{
    int thread = 1024;
    int block = (dim + thread - 1) / thread;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(vector_div_constant_kernel<double>), dim3(block), dim3(thread), 0, 0, dim, result, vector, constant);
}
// vector operator: result[i] = vector[i] / constant
template <typename FPTYPE>
inline void vector_div_constant_complex_wrapper(
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
    hipLaunchKernelGGL(HIP_KERNEL_NAME(vector_div_constant_kernel<thrust::complex<FPTYPE>>), dim3(block), dim3(thread), 0, 0, dim, result_tmp, vector_tmp, constant);
}
template <>
void vector_div_constant_op<std::complex<float>, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int dim,
    std::complex<float>* result,
    const std::complex<float>* vector,
    const float constant)
{
    vector_div_constant_complex_wrapper(d, dim, result, vector, constant);
}
template <>
void vector_div_constant_op<std::complex<double>, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int dim,
    std::complex<double>* result,
    const std::complex<double>* vector,
    const double constant)
{
    vector_div_constant_complex_wrapper(d, dim, result, vector, constant);
}
// vector operator: result[i] = vector1[i](not complex) * vector2[i](not complex)
template <>
void vector_mul_vector_op<double, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& dim,
    double* result,
    const double* vector1,
    const double* vector2)
{
    int thread = 1024;
    int block = (dim + thread - 1) / thread;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(vector_mul_vector_kernel<double>), dim3(block), dim3(thread), 0, 0, dim, result, vector1, vector2);
}

// vector operator: result[i] = vector1[i](complex) * vector2[i](not complex)
template <typename FPTYPE>
inline void vector_mul_vector_complex_wrapper(
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
    hipLaunchKernelGGL(HIP_KERNEL_NAME(vector_mul_vector_kernel<thrust::complex<FPTYPE>>), dim3(block), dim3(thread), 0, 0, dim, result_tmp, vector1_tmp, vector2);
}
template <>
void vector_mul_vector_op<std::complex<float>, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& dim,
    std::complex<float>* result,
    const std::complex<float>* vector1,
    const float* vector2)
{
    vector_mul_vector_complex_wrapper(d, dim, result, vector1, vector2);
}
template <>
void vector_mul_vector_op<std::complex<double>, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& dim,
    std::complex<double>* result,
    const std::complex<double>* vector1,
    const double* vector2)
{
    vector_mul_vector_complex_wrapper(d, dim, result, vector1, vector2);
}
// vector operator: result[i] = vector1[i](complex) / vector2[i](not complex)
template <>
void vector_div_vector_op<double, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& dim,
    double* result,
    const double* vector1,
    const double* vector2)
{
    int thread = 1024;
    int block = (dim + thread - 1) / thread;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(vector_div_vector_kernel<double>), dim3(block), dim3(thread), 0, 0, dim, result, vector1, vector2);
}
// vector operator: result[i] = vector1[i](complex) / vector2[i](not complex)
template <typename FPTYPE> 
inline void vector_div_vector_op_complex_wrapper(
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
    hipLaunchKernelGGL(HIP_KERNEL_NAME(vector_div_vector_kernel<thrust::complex<FPTYPE>>), dim3(block), dim3(thread), 0, 0, dim, result_tmp, vector1_tmp, vector2);
}
template <>
void vector_div_vector_op<std::complex<float>, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& dim,
    std::complex<float>* result,
    const std::complex<float>* vector1,
    const float* vector2)
{
    vector_div_vector_op_complex_wrapper(d, dim, result, vector1, vector2);
}
template <>
void vector_div_vector_op<std::complex<double>, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& dim,
    std::complex<double>* result,
    const std::complex<double>* vector1,
    const double* vector2)
{
    vector_div_vector_op_complex_wrapper(d, dim, result, vector1, vector2);
}

// vector operator: result[i] = vector1[i] * constant1 + vector2[i] * constant2
template <typename T>
void constantvector_addORsub_constantVector_op<T, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& dim,
    T* result,
    const T* vector1,
    const Real constant1,
    const T* vector2,
    const Real constant2) 
{
    using Type = typename GetTypeThrust<T>::type;
    using Real = typename GetTypeReal<T>::type;
    
    auto result_tmp = reinterpret_cast<Type*>(result);
    auto vector1_tmp = reinterpret_cast<const Type*>(vector1);
    auto vector2_tmp = reinterpret_cast<const Type*>(vector2);

    int thread = 1024;
    int block = (dim + thread - 1) / thread;
    constantvector_addORsub_constantVector_kernel<Type, Real> <<<block, thread >>>(dim, result_tmp, vector1_tmp, constant1, vector2_tmp, constant2);
}

template <>
void axpy_op<double, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const int& N,
    const double* alpha,
    const double* X,
    const int& incX,
    double* Y,
    const int& incY)
{
    hipblasDaxpy(cublas_handle, N, alpha, X, incX, Y, incY);
}

template <>
void axpy_op<std::complex<float>, psi::DEVICE_GPU>::operator()(
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
void axpy_op<std::complex<double>, psi::DEVICE_GPU>::operator()(
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
void gemv_op<double, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* d,
    const char& trans,
    const int& m,
    const int& n,
    const double* alpha,
    const double* A,
    const int& lda,
    const double* X,
    const int& incx,
    const double* beta,
    double* Y,
    const int& incy)
{
    hipblasOperation_t cutrans = {};
    if (trans == 'N') {
        cutrans = HIPBLAS_OP_N;
    }
    else if (trans == 'T') {
        cutrans = HIPBLAS_OP_T;
    }
    else if (trans == 'C') {
        cutrans = HIPBLAS_OP_C;
    }
    hipblasDgemv(cublas_handle, cutrans, m, n, alpha, A, lda, X, incx, beta, Y, incx);
}

template <>
void gemv_op<std::complex<float>, psi::DEVICE_GPU>::operator()(
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
    hipblasCgemv(cublas_handle, cutrans, m, n, (hipblasComplex*)alpha, (hipblasComplex*)A, lda, (hipblasComplex*)X, incx, (hipblasComplex*)beta, (hipblasComplex*)Y, incx);
}

template <> 
void gemv_op<std::complex<double>, psi::DEVICE_GPU>::operator()(
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
void gemm_op<double, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
    const char& transa,
    const char& transb,
    const int& m,
    const int& n,
    const int& k,
    const double* alpha,
    const double* a,
    const int& lda,
    const double* b,
    const int& ldb,
    const double* beta,
    double* c,
    const int& ldc)
{
    hipblasOperation_t cutransA;
    hipblasOperation_t cutransB;
    // cutransA
    if (transa == 'N') {
        cutransA = HIPBLAS_OP_N;
    }
    else if (transa == 'T') {
        cutransA = HIPBLAS_OP_T;
    }
    // cutransB
    if (transb == 'N') {
        cutransB = HIPBLAS_OP_N;
    }
    else if (transb == 'T') {
        cutransB = HIPBLAS_OP_T;
    }
    hipblasDgemm(cublas_handle, cutransA, cutransB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

template <>
void gemm_op<std::complex<float>, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
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
void gemm_op<std::complex<double>, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
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
void matrixTranspose_op<double, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
    const int& row,
    const int& col,
    const double* input_matrix,
    double* output_matrix)
{
    double* device_temp = nullptr;
    psi::memory::resize_memory_op<double, psi::DEVICE_GPU>()(d, device_temp, row * col);

    if (row == col)
    {
        double ONE = 1.0, ZERO = 0.0;
        // use 'geam' API todo transpose.
        hipblasDgeam(cublas_handle, HIPBLAS_OP_T, HIPBLAS_OP_N, col, row, &ONE, input_matrix, col, &ZERO, input_matrix, col, device_temp, col);
    }
    else
    {
        int thread = 1024;
        int block = (row + col + thread - 1) / thread;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(matrix_transpose_kernel<double>), dim3(block), dim3(thread), 0, 0, row, col, input_matrix, device_temp);
    }

    psi::memory::synchronize_memory_op<double, psi::DEVICE_GPU, psi::DEVICE_GPU>()(d, d, output_matrix, device_temp, row * col);

    psi::memory::delete_memory_op<double, psi::DEVICE_GPU>()(d, device_temp);

}

template <>
void matrixTranspose_op<std::complex<float>, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
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
        hipLaunchKernelGGL(HIP_KERNEL_NAME(matrix_transpose_kernel<thrust::complex<float>>), dim3(block), dim3(thread), 0, 0, row, col, (thrust::complex<float>*)input_matrix, (thrust::complex<float>*)device_temp);
    }

    psi::memory::synchronize_memory_op<std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_GPU>()(d, d, output_matrix, device_temp, row * col);

    psi::memory::delete_memory_op<std::complex<float>, psi::DEVICE_GPU>()(d, device_temp);

}

template <>
void matrixTranspose_op<std::complex<double>, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
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
        hipLaunchKernelGGL(HIP_KERNEL_NAME(matrix_transpose_kernel<thrust::complex<double>>), dim3(block), dim3(thread), 0, 0, row, col, (thrust::complex<double>*)input_matrix, (thrust::complex<double>*)device_temp);
    }
    
    psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_GPU>()(d, d, output_matrix, device_temp, row * col);

    psi::memory::delete_memory_op<std::complex<double>, psi::DEVICE_GPU>()(d, device_temp);
    
}

template <>
void matrixSetToAnother < double, psi::DEVICE_GPU > ::operator()(
    const psi::DEVICE_GPU* d,
    const int& n,
    const double* A,
    const int& LDA,
    double* B,
    const int& LDB)
{
    int thread = 1024;
    int block = (LDA + thread - 1) / thread;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(matrix_setTo_another_kernel<double>), dim3(block), dim3(thread), 0, 0, n, LDA, LDB, A, B);
}
template <>
void matrixSetToAnother < std::complex<float>, psi::DEVICE_GPU > ::operator()(
            const psi::DEVICE_GPU* d,
            const int& n,
    const std::complex<float>* A,
            const int& LDA,
    std::complex<float>* B,
            const int& LDB)
{
    int thread = 1024;
    int block = (LDA + thread - 1) / thread;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(matrix_setTo_another_kernel<thrust::complex<float>>), dim3(block), dim3(thread), 0, 0, n, LDA, LDB, reinterpret_cast<const thrust::complex<float>*>(A), reinterpret_cast<thrust::complex<float>*>(B));
}
template <>
void matrixSetToAnother < std::complex<double>, psi::DEVICE_GPU > ::operator()(
    const psi::DEVICE_GPU* d,
    const int& n,
    const std::complex<double>* A,
    const int& LDA,
    std::complex<double>* B,
    const int& LDB)
{
    int thread = 1024;
    int block = (LDA + thread - 1) / thread;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(matrix_setTo_another_kernel<thrust::complex<double>>), dim3(block), dim3(thread), 0, 0, n, LDA, LDB, reinterpret_cast<const thrust::complex<double>*>(A), reinterpret_cast<thrust::complex<double>*>(B));
}



// Explicitly instantiate functors for the types of functor registered.
template struct dot_real_op<std::complex<float>, psi::DEVICE_GPU>;
template struct calc_grad_with_block_op<std::complex<float>, psi::DEVICE_GPU>;
template struct line_minimize_with_block_op<std::complex<float>, psi::DEVICE_GPU>;
template struct vector_div_constant_op<std::complex<float>, psi::DEVICE_GPU>;
template struct vector_mul_vector_op<std::complex<float>, psi::DEVICE_GPU>;
template struct vector_div_vector_op<std::complex<float>, psi::DEVICE_GPU>;
template struct constantvector_addORsub_constantVector_op<std::complex<float>, psi::DEVICE_GPU>;
template struct matrixSetToAnother<std::complex<float>, psi::DEVICE_GPU>;

template struct dot_real_op<std::complex<double>, psi::DEVICE_GPU>;
template struct calc_grad_with_block_op<std::complex<double>, psi::DEVICE_GPU>;
template struct line_minimize_with_block_op<std::complex<double>, psi::DEVICE_GPU>;
template struct vector_div_constant_op<std::complex<double>, psi::DEVICE_GPU>;
template struct vector_mul_vector_op<std::complex<double>, psi::DEVICE_GPU>;
template struct vector_div_vector_op<std::complex<double>, psi::DEVICE_GPU>;
template struct constantvector_addORsub_constantVector_op<std::complex<double>, psi::DEVICE_GPU>;
template struct matrixSetToAnother<std::complex<double>, psi::DEVICE_GPU>;

#ifdef __LCAO
template struct dot_real_op<double, psi::DEVICE_GPU>;
template struct vector_div_constant_op<double, psi::DEVICE_GPU>;
template struct vector_mul_vector_op<double, psi::DEVICE_GPU>;
template struct vector_div_vector_op<double, psi::DEVICE_GPU>;
template struct matrixSetToAnother<double, psi::DEVICE_GPU>;
template struct constantvector_addORsub_constantVector_op<double, psi::DEVICE_GPU>;
#endif
}  // namespace hsolver
