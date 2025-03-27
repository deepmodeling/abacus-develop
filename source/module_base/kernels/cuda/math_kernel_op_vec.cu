#include "module_base/kernels/math_kernel_op.h"

#include <base/macros/macros.h>
#include <thrust/complex.h>

template <>
struct GetTypeReal<thrust::complex<float>> {
    using type = float; /**< The return type specialization for std::complex<double>. */
};
template <>
struct GetTypeReal<thrust::complex<double>> {
    using type = double; /**< The return type specialization for std::complex<double>. */
};
namespace ModuleBase
{
const int thread_per_block = 256;
void xdot_wrapper(const int &n, const float * x, const int &incx, const float * y, const int &incy, float &result);
void xdot_wrapper(const int &n, const double * x, const int &incx, const double * y, const int &incy, double &result);

// Define the CUDA kernel:
template <typename T>
__global__ void vector_mul_real_kernel(const int size,
                                       T* result,
                                       const T* vector,
                                       const typename GetTypeReal<T>::type constant)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size)
    {
        result[i] = vector[i] * constant;
    }
}

template <typename T>
__global__ void vector_mul_vector_kernel(const int size,
                                         T* result,
                                         const T* vector1,
                                         const typename GetTypeReal<T>::type* vector2,
                                         const bool add)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size)
    {
        if (add)
        {
            result[i] += vector1[i] * vector2[i];
        }
        else
        {
            result[i] = vector1[i] * vector2[i];
        }
    }
}

template <typename T>
__global__ void vector_div_vector_kernel(const int size,
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
__global__ void constantvector_addORsub_constantVector_kernel(const int size,
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

// vector operator: result[i] = vector[i] * constant
template <>
void vector_mul_real_op<double, base_device::DEVICE_GPU>::operator()(const int dim,
                                                                     double* result,
                                                                     const double* vector,
                                                                     const double constant)
{
    // In small cases, 1024 threads per block will only utilize 17 blocks, much less than 40
    int thread = thread_per_block;
    int block = (dim + thread - 1) / thread;
    vector_mul_real_kernel<double><<<block, thread>>>(dim, result, vector, constant);

    cudaCheckOnDebug();
}

template <typename FPTYPE>
inline void vector_mul_real_wrapper(const int dim,
                                    std::complex<FPTYPE>* result,
                                    const std::complex<FPTYPE>* vector,
                                    const FPTYPE constant)
{
    thrust::complex<FPTYPE>* result_tmp = reinterpret_cast<thrust::complex<FPTYPE>*>(result);
    const thrust::complex<FPTYPE>* vector_tmp = reinterpret_cast<const thrust::complex<FPTYPE>*>(vector);

    int thread = thread_per_block;
    int block = (dim + thread - 1) / thread;
    vector_mul_real_kernel<thrust::complex<FPTYPE>><<<block, thread>>>(dim, result_tmp, vector_tmp, constant);

    cudaCheckOnDebug();
}
template <>
void vector_mul_real_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const int dim,
                                                                                  std::complex<float>* result,
                                                                                  const std::complex<float>* vector,
                                                                                  const float constant)
{
    vector_mul_real_wrapper(dim, result, vector, constant);
}
template <>
void vector_mul_real_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(const int dim,
                                                                                   std::complex<double>* result,
                                                                                   const std::complex<double>* vector,
                                                                                   const double constant)
{
    vector_mul_real_wrapper(dim, result, vector, constant);
}

// vector operator: result[i] = vector1[i](not complex) * vector2[i](not complex)
template <>
void vector_mul_vector_op<double, base_device::DEVICE_GPU>::operator()(const int& dim,
                                                                       double* result,
                                                                       const double* vector1,
                                                                       const double* vector2,
                                                                       const bool& add)
{
    int thread = thread_per_block;
    int block = (dim + thread - 1) / thread;
    vector_mul_vector_kernel<double><<<block, thread>>>(dim, result, vector1, vector2, add);

    cudaCheckOnDebug();
}
// vector operator: result[i] = vector1[i](complex) * vector2[i](not complex)
template <typename FPTYPE>
inline void vector_mul_vector_complex_wrapper(const int& dim,
                                              std::complex<FPTYPE>* result,
                                              const std::complex<FPTYPE>* vector1,
                                              const FPTYPE* vector2,
                                              const bool& add)
{
    thrust::complex<FPTYPE>* result_tmp = reinterpret_cast<thrust::complex<FPTYPE>*>(result);
    const thrust::complex<FPTYPE>* vector1_tmp = reinterpret_cast<const thrust::complex<FPTYPE>*>(vector1);
    int thread = thread_per_block;
    int block = (dim + thread - 1) / thread;
    vector_mul_vector_kernel<thrust::complex<FPTYPE>><<<block, thread>>>(dim, result_tmp, vector1_tmp, vector2, add);

    cudaCheckOnDebug();
}
template <>
void vector_mul_vector_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const int& dim,
                                                                                    std::complex<float>* result,
                                                                                    const std::complex<float>* vector1,
                                                                                    const float* vector2,
                                                                                    const bool& add)
{
    vector_mul_vector_complex_wrapper(dim, result, vector1, vector2, add);
}
template <>
void vector_mul_vector_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(
    const int& dim,
    std::complex<double>* result,
    const std::complex<double>* vector1,
    const double* vector2,
    const bool& add)
{
    vector_mul_vector_complex_wrapper(dim, result, vector1, vector2, add);
}

// vector operator: result[i] = vector1[i](not complex) / vector2[i](not complex)
template <>
void vector_div_vector_op<double, base_device::DEVICE_GPU>::operator()(const int& dim,
                                                                       double* result,
                                                                       const double* vector1,
                                                                       const double* vector2)
{
    int thread = thread_per_block;
    int block = (dim + thread - 1) / thread;
    vector_div_vector_kernel<double><<<block, thread>>>(dim, result, vector1, vector2);

    cudaCheckOnDebug();
}
// vector operator: result[i] = vector1[i](complex) / vector2[i](not complex)
template <typename FPTYPE>
inline void vector_div_vector_complex_wrapper(const int& dim,
                                              std::complex<FPTYPE>* result,
                                              const std::complex<FPTYPE>* vector1,
                                              const FPTYPE* vector2)
{
    thrust::complex<FPTYPE>* result_tmp = reinterpret_cast<thrust::complex<FPTYPE>*>(result);
    const thrust::complex<FPTYPE>* vector1_tmp = reinterpret_cast<const thrust::complex<FPTYPE>*>(vector1);
    int thread = thread_per_block;
    int block = (dim + thread - 1) / thread;
    vector_div_vector_kernel<thrust::complex<FPTYPE>><<<block, thread>>>(dim, result_tmp, vector1_tmp, vector2);

    cudaCheckOnDebug();
}
template <>
void vector_div_vector_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const int& dim,
                                                                                    std::complex<float>* result,
                                                                                    const std::complex<float>* vector1,
                                                                                    const float* vector2)
{
    vector_div_vector_complex_wrapper(dim, result, vector1, vector2);
}
template <>
void vector_div_vector_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(
    const int& dim,
    std::complex<double>* result,
    const std::complex<double>* vector1,
    const double* vector2)
{
    vector_div_vector_complex_wrapper(dim, result, vector1, vector2);
}

// vector operator: result[i] = vector1[i] * constant1 + vector2[i] * constant2
template <typename T>
void vector_add_vector_op<T, base_device::DEVICE_GPU>::operator()(const int& dim,
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
    constantvector_addORsub_constantVector_kernel<Type, Real>
        <<<block, thread>>>(dim, result_tmp, vector1_tmp, constant1, vector2_tmp, constant2);

    cudaCheckOnDebug();
}

template <>
double dot_real_op<double, base_device::DEVICE_GPU>::operator()(const int& dim,
                                                                const double* psi_L,
                                                                const double* psi_R,
                                                                const bool reduce)
{
    double result = 0.0;
    xdot_wrapper(dim, psi_L, 1, psi_R, 1, result);
    if (reduce)
    {
        Parallel_Reduce::reduce_pool(result);
    }
    return result;
}
// for this implementation, please check
// https://thrust.github.io/doc/group__transformed__reductions_ga321192d85c5f510e52300ae762c7e995.html denghui modify
// 2022-10-03 Note that ddot_(2*dim,a,1,b,1) = REAL( zdotc_(dim,a,1,b,1) ) GPU specialization of actual computation.
template <typename FPTYPE>
inline FPTYPE dot_complex_wrapper(const int& dim,
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
    if (reduce)
    {
        Parallel_Reduce::reduce_pool(result);
    }
    return result;
}

template <>
float dot_real_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const int& dim,
                                                                            const std::complex<float>* psi_L,
                                                                            const std::complex<float>* psi_R,
                                                                            const bool reduce)
{
    return dot_complex_wrapper(dim, psi_L, psi_R, reduce);
}
template <>
double dot_real_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(const int& dim,
                                                                              const std::complex<double>* psi_L,
                                                                              const std::complex<double>* psi_R,
                                                                              const bool reduce)
{
    return dot_complex_wrapper(dim, psi_L, psi_R, reduce);
}

// Explicitly instantiate functors for the types of functor registered.
template struct vector_mul_real_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct vector_mul_real_op<double, base_device::DEVICE_GPU>;
template struct vector_mul_real_op<std::complex<double>, base_device::DEVICE_GPU>;

template struct vector_mul_vector_op<float, base_device::DEVICE_GPU>;
template struct vector_mul_vector_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct vector_mul_vector_op<double, base_device::DEVICE_GPU>;
template struct vector_mul_vector_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct vector_div_vector_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct vector_div_vector_op<double, base_device::DEVICE_GPU>;
template struct vector_div_vector_op<std::complex<double>, base_device::DEVICE_GPU>;

template struct vector_add_vector_op<float, base_device::DEVICE_GPU>;
template struct vector_add_vector_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct vector_add_vector_op<double, base_device::DEVICE_GPU>;
template struct vector_add_vector_op<std::complex<double>, base_device::DEVICE_GPU>;

template struct dot_real_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct dot_real_op<double, base_device::DEVICE_GPU>;
template struct dot_real_op<std::complex<double>, base_device::DEVICE_GPU>;
} // namespace ModuleBase