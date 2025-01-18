#include "module_base/module_device/memory_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_psi/psi.h"
#include "module_base/tool_quit.h"

#include <base/macros/macros.h>
#include <cuda_runtime.h>
#include <thrust/complex.h>
#include <thrust/execution_policy.h>
#include <thrust/inner_product.h>

namespace hsolver
{
const int warp_size = 32;
// const unsigned int full_mask = 0xffffffff;
const int thread_per_block = 256;
}

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

static cublasHandle_t cublas_handle = nullptr;

static inline
void xdot_wrapper(const int &n, const float * x, const int &incx, const float * y, const int &incy, float &result) {
    cublasErrcheck(cublasSdot(cublas_handle, n, x, incx, y, incy, &result));
}

static inline
void xdot_wrapper(const int &n, const double * x, const int &incx, const double * y, const int &incy, double &result) {
    cublasErrcheck(cublasDdot(cublas_handle, n, x, incx, y, incy, &result));
}

void createGpuBlasHandle(){
    if (cublas_handle == nullptr) {
        cublasErrcheck(cublasCreate(&cublas_handle));
    }
}

void destoryBLAShandle(){
    if (cublas_handle != nullptr) {
        cublasErrcheck(cublasDestroy(cublas_handle));
        cublas_handle = nullptr;
    }
}

// template <typename FPTYPE>
// __forceinline__ __device__ void warp_reduce(FPTYPE& val) {
//     for (int offset = 16; offset > 0; offset >>= 1)
//         val += __shfl_down_sync(full_mask, val, offset);
// }

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
    __shared__ Real data[thread_per_block * 3];

    data[tid] = 0;

    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += thread_per_block) {
        item = band_idx * n_basis_max + basis_idx;
        data[tid] += (grad[item] * thrust::conj(grad[item])).real();
    }
    __syncthreads();
    // just do some parallel reduction in shared memory
    for (int ii = thread_per_block >> 1; ii > warp_size; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
        }
        __syncthreads();
    }

    // For threads in the same warp, it is better that they process the same work
    // Also, __syncwarp() should be used instead of __syncthreads()
    // Therefore we unroll the loop and ensure that the threads does the same work
    if (tid < warp_size) {
        data[tid] += data[tid + 32]; __syncwarp();
        data[tid] += data[tid + 16]; __syncwarp();
        data[tid] += data[tid + 8]; __syncwarp();
        data[tid] += data[tid + 4]; __syncwarp();
        data[tid] += data[tid + 2]; __syncwarp();
        data[tid] += data[tid + 1]; __syncwarp();
    }

    __syncthreads();

    Real norm = 1.0 / sqrt(data[0]);
    __syncthreads();

    data[tid] = 0;
    data[thread_per_block + tid] = 0;
    data[2 * thread_per_block + tid] = 0;
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += thread_per_block) {
        item = band_idx * n_basis_max + basis_idx;
        grad[item] *= norm;
        hgrad[item] *= norm;
        data[tid] += (hpsi[item] * thrust::conj(psi[item])).real();
        data[thread_per_block + tid] += (grad[item] * thrust::conj(hpsi[item])).real();
        data[2 * thread_per_block + tid] += (grad[item] * thrust::conj(hgrad[item])).real();
    }
    __syncthreads();

    // just do some parallel reduction in shared memory
    for (int ii = thread_per_block >> 1; ii > warp_size; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
            data[thread_per_block + tid] += data[thread_per_block + tid + ii];
            data[2 * thread_per_block + tid] += data[2 * thread_per_block + tid + ii];
        }
        __syncthreads();
    }
    if (tid < warp_size) {
        data[tid] += data[tid + 32]; __syncwarp();
        data[tid] += data[tid + 16]; __syncwarp();
        data[tid] += data[tid + 8]; __syncwarp();
        data[tid] += data[tid + 4]; __syncwarp();
        data[tid] += data[tid + 2]; __syncwarp();
        data[tid] += data[tid + 1]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 32]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 16]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 8]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 4]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 2]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 1]; __syncwarp();
        data[2 * thread_per_block + tid] += data[2 * thread_per_block + tid + 32]; __syncwarp();
        data[2 * thread_per_block + tid] += data[2 * thread_per_block + tid + 16]; __syncwarp();
        data[2 * thread_per_block + tid] += data[2 * thread_per_block + tid + 8]; __syncwarp();
        data[2 * thread_per_block + tid] += data[2 * thread_per_block + tid + 4]; __syncwarp();
        data[2 * thread_per_block + tid] += data[2 * thread_per_block + tid + 2]; __syncwarp();
        data[2 * thread_per_block + tid] += data[2 * thread_per_block + tid + 1]; __syncwarp();
    }
    __syncthreads();
    epsilo_0 = data[0];
    epsilo_1 = data[thread_per_block];
    epsilo_2 = data[2 * thread_per_block];

    theta = 0.5 * abs(atan(2 * epsilo_1/(epsilo_0 - epsilo_2)));
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += thread_per_block) {
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
    __shared__ Real data[thread_per_block * 2];

    // Init shared memory
    data[tid] = 0;

    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += thread_per_block) {
        item = band_idx * n_basis_max + basis_idx;
        data[tid] += (psi[item] * thrust::conj(psi[item])).real();
    }
    __syncthreads();
    // just do some parallel reduction in shared memory
    for (int ii = thread_per_block >> 1; ii > warp_size; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
        }
        __syncthreads();
    }

    if (tid < warp_size) {
        data[tid] += data[tid + 32]; __syncwarp();
        data[tid] += data[tid + 16]; __syncwarp();
        data[tid] += data[tid + 8]; __syncwarp();
        data[tid] += data[tid + 4]; __syncwarp();
        data[tid] += data[tid + 2]; __syncwarp();
        data[tid] += data[tid + 1]; __syncwarp();
    }

    __syncthreads();

    Real norm = 1.0 / sqrt(data[0]);
    __syncthreads();

    data[tid] = 0;
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += thread_per_block) {
        item = band_idx * n_basis_max + basis_idx;
        psi[item] *= norm;
        hpsi[item] *= norm;
        data[tid] += (hpsi[item] * thrust::conj(psi[item])).real();
    }
    __syncthreads();

    // just do some parallel reduction in shared memory
    for (int ii = thread_per_block >> 1; ii > warp_size; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
        }
        __syncthreads();
    }

    if (tid < warp_size) {
        data[tid] += data[tid + 32]; __syncwarp();
        data[tid] += data[tid + 16]; __syncwarp();
        data[tid] += data[tid + 8]; __syncwarp();
        data[tid] += data[tid + 4]; __syncwarp();
        data[tid] += data[tid + 2]; __syncwarp();
        data[tid] += data[tid + 1]; __syncwarp();
    }

    __syncthreads();
    epsilo = data[0];
    __syncthreads();

    data[tid] = 0;
    data[thread_per_block + tid] = 0;
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += thread_per_block) {
        item = band_idx * n_basis_max + basis_idx;
        grad_1 = hpsi[item] - epsilo * psi[item];
        grad_2 = thrust::norm(grad_1);
        data[tid] += grad_2;
        data[thread_per_block + tid] += grad_2 / prec[basis_idx];
    }
    __syncthreads();

    // just do some parallel reduction in shared memory
    for (int ii = thread_per_block >> 1; ii > warp_size; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
            data[thread_per_block + tid] += data[thread_per_block + tid + ii];
        }
        __syncthreads();
    }

    if (tid < warp_size) {
        data[tid] += data[tid + 32]; __syncwarp();
        data[tid] += data[tid + 16]; __syncwarp();
        data[tid] += data[tid + 8]; __syncwarp();
        data[tid] += data[tid + 4]; __syncwarp();
        data[tid] += data[tid + 2]; __syncwarp();
        data[tid] += data[tid + 1]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 32]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 16]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 8]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 4]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 2]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 1]; __syncwarp();
    }

    __syncthreads();
    err_st = data[0];
    beta_st = data[thread_per_block];
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += thread_per_block) {
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
void line_minimize_with_block_op<T, base_device::DEVICE_GPU>::operator()(T* grad_out,
                                                                         T* hgrad_out,
                                                                         T* psi_out,
                                                                         T* hpsi_out,
                                                                         const int& n_basis,
                                                                         const int& n_basis_max,
                                                                         const int& n_band)
{
    auto A = reinterpret_cast<thrust::complex<Real>*>(grad_out);
    auto B = reinterpret_cast<thrust::complex<Real>*>(hgrad_out);
    auto C = reinterpret_cast<thrust::complex<Real>*>(psi_out);
    auto D = reinterpret_cast<thrust::complex<Real>*>(hpsi_out);

    line_minimize_with_block<Real><<<n_band, thread_per_block>>>(
            A, B, C, D,
            n_basis, n_basis_max);

    cudaCheckOnDebug();
}

template <typename T>
void calc_grad_with_block_op<T, base_device::DEVICE_GPU>::operator()(const Real* prec_in,
                                                                     Real* err_out,
                                                                     Real* beta_out,
                                                                     T* psi_out,
                                                                     T* hpsi_out,
                                                                     T* grad_out,
                                                                     T* grad_old_out,
                                                                     const int& n_basis,
                                                                     const int& n_basis_max,
                                                                     const int& n_band)
{
    auto A = reinterpret_cast<thrust::complex<Real>*>(psi_out);
    auto B = reinterpret_cast<thrust::complex<Real>*>(hpsi_out);
    auto C = reinterpret_cast<thrust::complex<Real>*>(grad_out);
    auto D = reinterpret_cast<thrust::complex<Real>*>(grad_old_out);

    calc_grad_with_block<Real><<<n_band, thread_per_block>>>(
            prec_in, err_out, beta_out,
            A, B, C, D,
            n_basis, n_basis_max);

    cudaCheckOnDebug();
}

template <>
double dot_real_op<double, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
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
// 2022-10-03 Note that ddot_(2*dim,a,1,b,1) = REAL( zdotc_(dim,a,1,b,1) ) GPU specialization of actual computation.
template <typename FPTYPE>
inline FPTYPE dot_complex_wrapper(const base_device::DEVICE_GPU* d,
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
float dot_real_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                            const int& dim,
                                                                            const std::complex<float>* psi_L,
                                                                            const std::complex<float>* psi_R,
                                                                            const bool reduce)
{
    return dot_complex_wrapper(d, dim, psi_L, psi_R, reduce);
}
template <>
double dot_real_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                              const int& dim,
                                                                              const std::complex<double>* psi_L,
                                                                              const std::complex<double>* psi_R,
                                                                              const bool reduce)
{
    return dot_complex_wrapper(d, dim, psi_L, psi_R, reduce);
}

// vector operator: result[i] = vector[i] / constant
template <>
void vector_div_constant_op<double, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                         const int dim,
                                                                         double* result,
                                                                         const double* vector,
                                                                         const double constant)
{
    // In small cases, 1024 threads per block will only utilize 17 blocks, much less than 40
    int thread = thread_per_block;
    int block = (dim + thread - 1) / thread;
    vector_div_constant_kernel<double> <<<block, thread >>> (dim, result, vector, constant);

    cudaCheckOnDebug();
}

// vector operator: result[i] = vector[i] / constant
template <typename FPTYPE>
inline void vector_div_constant_complex_wrapper(const base_device::DEVICE_GPU* d,
                                                const int dim,
                                                std::complex<FPTYPE>* result,
                                                const std::complex<FPTYPE>* vector,
                                                const FPTYPE constant)
{
    thrust::complex<FPTYPE>* result_tmp = reinterpret_cast<thrust::complex<FPTYPE>*>(result);
    const thrust::complex<FPTYPE>* vector_tmp = reinterpret_cast<const thrust::complex<FPTYPE>*>(vector);

    int thread = thread_per_block;
    int block = (dim + thread - 1) / thread;
    vector_div_constant_kernel<thrust::complex<FPTYPE>> <<<block, thread >>> (dim, result_tmp, vector_tmp, constant);

    cudaCheckOnDebug();
}
template <>
void vector_div_constant_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                                      const int dim,
                                                                                      std::complex<float>* result,
                                                                                      const std::complex<float>* vector,
                                                                                      const float constant)
{
    vector_div_constant_complex_wrapper(d, dim, result, vector, constant);
}
template <>
void vector_div_constant_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(
    const base_device::DEVICE_GPU* d,
    const int dim,
    std::complex<double>* result,
    const std::complex<double>* vector,
    const double constant)
{
    vector_div_constant_complex_wrapper(d, dim, result, vector, constant);
}
// vector operator: result[i] = vector1[i](not complex) * vector2[i](not complex)
template <>
void vector_mul_vector_op<double, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                       const int& dim,
                                                                       double* result,
                                                                       const double* vector1,
                                                                       const double* vector2)
{
    int thread = thread_per_block;
    int block = (dim + thread - 1) / thread;
    vector_mul_vector_kernel<double> <<<block, thread >>> (dim, result, vector1, vector2);

    cudaCheckOnDebug();
}
// vector operator: result[i] = vector1[i](complex) * vector2[i](not complex)
template <typename FPTYPE>
inline void vector_mul_vector_complex_wrapper(const base_device::DEVICE_GPU* d,
                                              const int& dim,
                                              std::complex<FPTYPE>* result,
                                              const std::complex<FPTYPE>* vector1,
                                              const FPTYPE* vector2)
{
    thrust::complex<FPTYPE>* result_tmp = reinterpret_cast<thrust::complex<FPTYPE>*>(result);
    const thrust::complex<FPTYPE>* vector1_tmp = reinterpret_cast<const thrust::complex<FPTYPE>*>(vector1);
    int thread = thread_per_block;
    int block = (dim + thread - 1) / thread;
    vector_mul_vector_kernel<thrust::complex<FPTYPE>> <<<block, thread >>> (dim, result_tmp, vector1_tmp, vector2);

    cudaCheckOnDebug();
}
template <>
void vector_mul_vector_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                                    const int& dim,
                                                                                    std::complex<float>* result,
                                                                                    const std::complex<float>* vector1,
                                                                                    const float* vector2)
{
    vector_mul_vector_complex_wrapper(d, dim, result, vector1, vector2);
}
template <>
void vector_mul_vector_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(
    const base_device::DEVICE_GPU* d,
    const int& dim,
    std::complex<double>* result,
    const std::complex<double>* vector1,
    const double* vector2)
{
    vector_mul_vector_complex_wrapper(d, dim, result, vector1, vector2);
}

// vector operator: result[i] = vector1[i](not complex) / vector2[i](not complex)
template <>
void vector_div_vector_op<double, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                       const int& dim,
                                                                       double* result,
                                                                       const double* vector1,
                                                                       const double* vector2)
{
    int thread = thread_per_block;
    int block = (dim + thread - 1) / thread;
    vector_div_vector_kernel<double> <<<block, thread >>> (dim, result, vector1, vector2);

    cudaCheckOnDebug();
}
// vector operator: result[i] = vector1[i](complex) / vector2[i](not complex)
template <typename FPTYPE>
inline void vector_div_vector_complex_wrapper(const base_device::DEVICE_GPU* d,
                                              const int& dim,
                                              std::complex<FPTYPE>* result,
                                              const std::complex<FPTYPE>* vector1,
                                              const FPTYPE* vector2)
{
    thrust::complex<FPTYPE>* result_tmp = reinterpret_cast<thrust::complex<FPTYPE>*>(result);
    const thrust::complex<FPTYPE>* vector1_tmp = reinterpret_cast<const thrust::complex<FPTYPE>*>(vector1);
    int thread = thread_per_block;
    int block = (dim + thread - 1) / thread;
    vector_div_vector_kernel<thrust::complex<FPTYPE>> <<<block, thread >>> (dim, result_tmp, vector1_tmp, vector2);

    cudaCheckOnDebug();
}
template <>
void vector_div_vector_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                                    const int& dim,
                                                                                    std::complex<float>* result,
                                                                                    const std::complex<float>* vector1,
                                                                                    const float* vector2)
{
    vector_div_vector_complex_wrapper(d, dim, result, vector1, vector2);
}
template <>
void vector_div_vector_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(
    const base_device::DEVICE_GPU* d,
    const int& dim,
    std::complex<double>* result,
    const std::complex<double>* vector1,
    const double* vector2)
{
    vector_div_vector_complex_wrapper(d, dim, result, vector1, vector2);
}
// vector operator: result[i] = vector1[i] * constant1 + vector2[i] * constant2
template <typename T>
void constantvector_addORsub_constantVector_op<T, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
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

    int thread = thread_per_block;
    int block = (dim + thread - 1) / thread;
    constantvector_addORsub_constantVector_kernel<Type, Real> <<<block, thread >>>(dim, result_tmp, vector1_tmp, constant1, vector2_tmp, constant2);

    cudaCheckOnDebug();
}

template <>
void axpy_op<double, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                          const int& N,
                                                          const double* alpha,
                                                          const double* X,
                                                          const int& incX,
                                                          double* Y,
                                                          const int& incY)
{
    cublasErrcheck(cublasDaxpy(cublas_handle, N, alpha, X, incX, Y, incY));
}

template <>
void axpy_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                       const int& N,
                                                                       const std::complex<float>* alpha,
                                                                       const std::complex<float>* X,
                                                                       const int& incX,
                                                                       std::complex<float>* Y,
                                                                       const int& incY)
{
    cublasErrcheck(cublasCaxpy(cublas_handle, N, (float2*)alpha, (float2*)X, incX, (float2*)Y, incY));
}

template <>
void axpy_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                        const int& N,
                                                                        const std::complex<double>* alpha,
                                                                        const std::complex<double>* X,
                                                                        const int& incX,
                                                                        std::complex<double>* Y,
                                                                        const int& incY)
{
    cublasErrcheck(cublasZaxpy(cublas_handle, N, (double2*)alpha, (double2*)X, incX, (double2*)Y, incY));
}

cublasOperation_t judge_trans_op(bool is_complex, const char& trans, const char* name)
{
    if (trans == 'N')
    {
        return CUBLAS_OP_N;
    }
    else if(trans == 'T')
    {
        return CUBLAS_OP_T;
    }
    else if(is_complex && trans == 'C')
    {
        return CUBLAS_OP_C;
    }
    else 
    {
        ModuleBase::WARNING_QUIT(name, std::string("Unknown trans type ") + trans + std::string(" !"));
    }
}

template <>
void gemv_op<double, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
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
    cublasOperation_t cutrans = judge_trans_op(false, trans, "gemv_op");
    cublasErrcheck(cublasDgemv(cublas_handle, cutrans, m, n, alpha, A, lda, X, incx, beta, Y, incx));
}

template <>
void gemv_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                       const char& trans,
                                                                       const int& m,
                                                                       const int& n,
                                                                       const std::complex<float>* alpha_in,
                                                                       const std::complex<float>* A,
                                                                       const int& lda,
                                                                       const std::complex<float>* X,
                                                                       const int& incx,
                                                                       const std::complex<float>* beta_in,
                                                                       std::complex<float>* Y,
                                                                       const int& incy)
{
    cublasOperation_t cutrans = judge_trans_op(true, trans, "gemv_op");
    cuFloatComplex alpha = make_cuFloatComplex(alpha_in->real(), alpha_in->imag());
    cuFloatComplex beta = make_cuFloatComplex(beta_in->real(), beta_in->imag());
    cublasErrcheck(cublasCgemv(cublas_handle, cutrans, m, n, &alpha, (cuFloatComplex*)A, lda, (cuFloatComplex*)X, incx, &beta, (cuFloatComplex*)Y, incx));
}

template <>
void gemv_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                        const char& trans,
                                                                        const int& m,
                                                                        const int& n,
                                                                        const std::complex<double>* alpha_in,
                                                                        const std::complex<double>* A,
                                                                        const int& lda,
                                                                        const std::complex<double>* X,
                                                                        const int& incx,
                                                                        const std::complex<double>* beta_in,
                                                                        std::complex<double>* Y,
                                                                        const int& incy)
{
    cublasOperation_t cutrans = judge_trans_op(true, trans, "gemv_op");
    cuDoubleComplex alpha = make_cuDoubleComplex(alpha_in->real(), alpha_in->imag());
    cuDoubleComplex beta = make_cuDoubleComplex(beta_in->real(), beta_in->imag());
    // icpc and nvcc have some compatible problems
    // We must use cuDoubleComplex instead of converting std::complex<double>* to cuDoubleComplex*
    cublasErrcheck(cublasZgemv(cublas_handle, cutrans, m, n, &alpha, (cuDoubleComplex*)A, lda, (cuDoubleComplex*)X, incx, &beta, (cuDoubleComplex*)Y, incx));
}

template <>
void scal_op<float, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                         const int& N,
                                                         const std::complex<float>* alpha,
                                                         std::complex<float>* X,
                                                         const int& incx)
{
    cublasErrcheck(cublasCscal(cublas_handle, N, (float2*)alpha, (float2*)X, incx));
}

template <>
void scal_op<double, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                          const int& N,
                                                          const std::complex<double>* alpha,
                                                          std::complex<double>* X,
                                                          const int& incx)
{
    cublasErrcheck(cublasZscal(cublas_handle, N, (double2*)alpha, (double2*)X, incx));
}

template <>
void gemm_op<double, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
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
    cublasOperation_t cutransA = judge_trans_op(false, transa, "gemm_op");
    cublasOperation_t cutransB = judge_trans_op(false, transb, "gemm_op");
    cublasErrcheck(cublasDgemm(cublas_handle, cutransA, cutransB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc));
}
template <>
void gemm_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                       const char& transa,
                                                                       const char& transb,
                                                                       const int& m,
                                                                       const int& n,
                                                                       const int& k,
                                                                       const std::complex<float>* alpha,
                                                                       const std::complex<float>* a,
                                                                       const int& lda,
                                                                       const std::complex<float>* b,
                                                                       const int& ldb,
                                                                       const std::complex<float>* beta,
                                                                       std::complex<float>* c,
                                                                       const int& ldc)
{
    cublasOperation_t cutransA = judge_trans_op(true, transa, "gemm_op");
    cublasOperation_t cutransB = judge_trans_op(true, transb, "gemm_op");
    cublasErrcheck(cublasCgemm(cublas_handle, cutransA, cutransB, m, n ,k, (float2*)alpha, (float2*)a , lda, (float2*)b, ldb, (float2*)beta, (float2*)c, ldc));
}

template <>
void gemm_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                        const char& transa,
                                                                        const char& transb,
                                                                        const int& m,
                                                                        const int& n,
                                                                        const int& k,
                                                                        const std::complex<double>* alpha,
                                                                        const std::complex<double>* a,
                                                                        const int& lda,
                                                                        const std::complex<double>* b,
                                                                        const int& ldb,
                                                                        const std::complex<double>* beta,
                                                                        std::complex<double>* c,
                                                                        const int& ldc)
{
    cublasOperation_t cutransA = judge_trans_op(true, transa, "gemm_op");
    cublasOperation_t cutransB = judge_trans_op(true, transb, "gemm_op");
    cublasErrcheck(cublasZgemm(cublas_handle, cutransA, cutransB, m, n ,k, (double2*)alpha, (double2*)a , lda, (double2*)b, ldb, (double2*)beta, (double2*)c, ldc));
}

template <>
void matrixTranspose_op<double, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                     const int& row,
                                                                     const int& col,
                                                                     const double* input_matrix,
                                                                     double* output_matrix)
{
    double* device_temp = nullptr;
    base_device::memory::resize_memory_op<double, base_device::DEVICE_GPU>()(device_temp, row * col);

    if (row == col)
    {
        double ONE = 1.0, ZERO = 0.0;

        // use 'geam' API todo transpose.
        cublasErrcheck(cublasDgeam(cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, col, row, &ONE, input_matrix, col, &ZERO, input_matrix, col, device_temp, col));
    }
    else
    {
        int thread = 1024;
        int block = (row + col + thread - 1) / thread;
        matrix_transpose_kernel<double> <<<block, thread >>> (row, col, input_matrix, device_temp);

        cudaCheckOnDebug();
    }

    base_device::memory::synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(
        output_matrix,
        device_temp,
        row * col);

    base_device::memory::delete_memory_op<double, base_device::DEVICE_GPU>()(device_temp);
}

template <>
void matrixTranspose_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(
    const base_device::DEVICE_GPU* d,
    const int& row,
    const int& col,
    const std::complex<float>* input_matrix,
    std::complex<float>* output_matrix)
{
    std::complex<float>* device_temp = nullptr;
    base_device::memory::resize_memory_op<std::complex<float>, base_device::DEVICE_GPU>()(device_temp, row * col);

    if (row == col)
    {
        double2 ONE, ZERO;
        ONE.x = 1.0;
        ONE.y = 0.0;
        ZERO.x = ZERO.y = 0.0;

        // use 'geam' API todo transpose.
        cublasErrcheck(cublasCgeam(cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, col, row,
                                   reinterpret_cast<const float2 *>(&ONE), (float2*)input_matrix, col,
                                   reinterpret_cast<const float2 *>(&ZERO), (float2*)input_matrix, col, (float2*)device_temp, col));
    } else
    {
        int thread = 1024;
        int block = (row + col + thread - 1) / thread;
        matrix_transpose_kernel<thrust::complex<float>> <<<block, thread >>> (row, col, (thrust::complex<float>*)input_matrix, (thrust::complex<float>*)device_temp);

        cudaCheckOnDebug();
    }

    base_device::memory::synchronize_memory_op<std::complex<float>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(
        output_matrix,
        device_temp,
        row * col);

    base_device::memory::delete_memory_op<std::complex<float>, base_device::DEVICE_GPU>()(device_temp);

    cudaCheckOnDebug();

}

template <>
void matrixTranspose_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(
    const base_device::DEVICE_GPU* d,
    const int& row,
    const int& col,
    const std::complex<double>* input_matrix,
    std::complex<double>* output_matrix)
{
    std::complex<double>* device_temp = nullptr;
    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(device_temp, row * col);

    if (row == col)
    {
        double2 ONE, ZERO;
        ONE.x = 1.0;
        ONE.y = 0.0;
        ZERO.x = ZERO.y = 0.0;

        // use 'geam' API todo transpose.
        cublasErrcheck(cublasZgeam(cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, col, row, &ONE, (double2*)input_matrix, col, &ZERO, (double2*)input_matrix, col, (double2*)device_temp, col));
    } else
    {
        int thread = 1024;
        int block = (row + col + thread - 1) / thread;
        matrix_transpose_kernel<thrust::complex<double>> <<<block, thread >>> (row, col, (thrust::complex<double>*)input_matrix, (thrust::complex<double>*)device_temp);
        cudaCheckOnDebug();
    }

    base_device::memory::synchronize_memory_op<std::complex<double>,
                                               base_device::DEVICE_GPU,
                                               base_device::DEVICE_GPU>()(output_matrix, device_temp, row * col);

    base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(device_temp);
}

template <>
void matrixSetToAnother<double, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                     const int& n,
                                                                     const double* A,
                                                                     const int& LDA,
                                                                     double* B,
                                                                     const int& LDB)
{
    int thread = 1024;
    int block = (LDA + thread - 1) / thread;
    matrix_setTo_another_kernel<double> <<<block, thread >>> (n, LDA, LDB, A, B);
    cudaCheckOnDebug();
}
template <>
void matrixSetToAnother<std::complex<float>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                                  const int& n,
                                                                                  const std::complex<float>* A,
                                                                                  const int& LDA,
                                                                                  std::complex<float>* B,
                                                                                  const int& LDB)
{
    int thread = 1024;
    int block = (LDA + thread - 1) / thread;
    matrix_setTo_another_kernel<thrust::complex<float>> <<<block, thread >>> (n, LDA, LDB, reinterpret_cast<const thrust::complex<float>*>(A), reinterpret_cast<thrust::complex<float>*>(B));
    cudaCheckOnDebug();
}
template <>
void matrixSetToAnother<std::complex<double>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                                   const int& n,
                                                                                   const std::complex<double>* A,
                                                                                   const int& LDA,
                                                                                   std::complex<double>* B,
                                                                                   const int& LDB)
{
    int thread = 1024;
    int block = (LDA + thread - 1) / thread;
    matrix_setTo_another_kernel<thrust::complex<double>> <<<block, thread >>> (n, LDA, LDB, reinterpret_cast<const thrust::complex<double>*>(A), reinterpret_cast<thrust::complex<double>*>(B));

    cudaCheckOnDebug();
}


// Explicitly instantiate functors for the types of functor registered.
template struct dot_real_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct calc_grad_with_block_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct line_minimize_with_block_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct vector_div_constant_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct vector_mul_vector_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct vector_div_vector_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct constantvector_addORsub_constantVector_op<float, base_device::DEVICE_GPU>;
template struct constantvector_addORsub_constantVector_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct matrixSetToAnother<std::complex<float>, base_device::DEVICE_GPU>;

template struct dot_real_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct calc_grad_with_block_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct line_minimize_with_block_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct vector_div_constant_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct vector_mul_vector_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct vector_div_vector_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct constantvector_addORsub_constantVector_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct matrixSetToAnother<std::complex<double>, base_device::DEVICE_GPU>;

#ifdef __LCAO
template struct dot_real_op<double, base_device::DEVICE_GPU>;
template struct vector_div_constant_op<double, base_device::DEVICE_GPU>;
template struct vector_mul_vector_op<double, base_device::DEVICE_GPU>;
template struct vector_div_vector_op<double, base_device::DEVICE_GPU>;
template struct matrixSetToAnother<double, base_device::DEVICE_GPU>;
template struct constantvector_addORsub_constantVector_op<double, base_device::DEVICE_GPU>;
#endif
}  // namespace hsolver
