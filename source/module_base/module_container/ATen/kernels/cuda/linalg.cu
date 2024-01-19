#include <ATen/kernels/linalg.h>
#include <ATen/core/tensor.h>
#include <ATen/core/tensor_types.h>
#include <base/macros/cuda.h>

#include <cuda_runtime.h>
#include <thrust/complex.h>

namespace container {
namespace kernels {

template <typename T>
__device__ static inline 
T conj(T& in) {
    return in;
}

template <typename T>
__device__ static inline 
thrust::complex<T> conj(thrust::complex<T>& in) {
    return thrust::conj(in);
}

template <typename T>
__global__ void do_add_kernel(
    const int num_element,
    const T alpha,
    const T* x,
    const T beta,
    const T* y,
    T* z)
{
    // Perform add operation for the specified range [begin, end) in the output Tensor.
    for (int o_idx = threadIdx.x; o_idx < num_element; o_idx += blockDim.x) {
        // Assign the sum of the input Tensor elements at index 'o_idx' to the output Tensor element at index 'o_idx'.
        z[o_idx] = alpha * x[o_idx] + beta * y[o_idx];
    }
}

template <typename T>
__global__ void do_mul_kernel(
    const int num_element,
    const T alpha,
    const T* x,
    T* y)
{
    for (int o_idx = threadIdx.x; o_idx < num_element; o_idx += blockDim.x) {
        // Assign the sum of the input Tensor elements at index 'o_idx' to the output Tensor element at index 'o_idx'.
        y[o_idx] = alpha * x[o_idx];
    }
}

template <typename T>
__global__ void do_mul_kernel(
    const int num_element,
    const T alpha,
    const T* x,
    const T* y,
    T* z)
{
    for (int o_idx = threadIdx.x; o_idx < num_element; o_idx += blockDim.x) {
        // Assign the sum of the input Tensor elements at index 'o_idx' to the output Tensor element at index 'o_idx'.
        z[o_idx] = alpha * x[o_idx] * y[o_idx];
    }
}

template <typename T>
__global__ void do_div_kernel(
    const int num_element,
    const T alpha,
    const T* x,
    const T* y,
    T* z)
{
    for (int o_idx = threadIdx.x; o_idx < num_element; o_idx += blockDim.x) {
        // Assign the sum of the input Tensor elements at index 'o_idx' to the output Tensor element at index 'o_idx'.
        z[o_idx] = alpha * x[o_idx] / y[o_idx];
    }
}

template <typename T>
__global__ void do_fma_kernel(
    const int num_element,
    const T alpha,
    const T* x,
    const T* y,
    const T beta,
    const T* z,
    T* out)
{
    for (int o_idx = threadIdx.x; o_idx < num_element; o_idx += blockDim.x) {
        // Assign the sum of the input Tensor elements at index 'o_idx' to the output Tensor element at index 'o_idx'.
        out[o_idx] = alpha * x[o_idx] * y[o_idx] + beta * z[o_idx];
    }
}

template <typename T, bool Conjugate>
__global__ void do_transpose_kernel(
    int ndim,
    int64_t num_elements,
    const T* p,
    const int* perm,
    const int64_t* in_strides,
    const int64_t* out_strides,
    T* q)
{
    for (int64_t o_idx = 0; o_idx < num_elements; o_idx++) {
        int64_t i_idx = 0; // Initialize the index for the input Tensor element.
        int64_t current_o_idx = o_idx; // Calculate the index for the output Tensor element.

        // Iterate over each dimension of the output Tensor.
        for (int ii = 0; ii < ndim; ++ii) {
            // Calculate the ratio of the current output Tensor index 'current_o_idx' in the current dimension.
            const int64_t ratio = current_o_idx / out_strides[ii];
            // Update the output Tensor index 'current_o_idx' by removing the offset in the current dimension.
            current_o_idx -= ratio * out_strides[ii];
            // Calculate the offset for the corresponding index position in the input Tensor and accumulate it in 'i_idx'.
            i_idx += ratio * in_strides[perm[ii]];
        }
        // Check if conjugation is needed.
        if (Conjugate) {
            // Assign the conjugate value of the input Tensor element at index 'i_idx' to the output Tensor element at index 'o_idx'.
            q[o_idx] = kernels::conj(p[i_idx]);
        } else {
            // Assign the input Tensor element at index 'i_idx' to the output Tensor element at index 'o_idx'.
            q[o_idx] = p[i_idx];
        }
    }
}

template <typename T>
__global__ void do_stride_kernel(
    int ndim,
    int64_t size,
    const T* p,
    const int64_t* stride,
    const int64_t* in_strides,
    const int64_t* out_strides,
    T* q)
{
    // Perform stride operation for the specified range [begin, end) in the output Tensor.
    for (int64_t o_idx = threadIdx.x; o_idx < size; o_idx += blockDim.x) {
        int64_t i_idx = 0; // Initialize the index for the input Tensor element.
        int64_t current_o_idx = o_idx; // Calculate the index for the output Tensor element.
        // Iterate over each dimension of the output Tensor.
        for (int ii = 0; ii < ndim; ++ii) {
            // Calculate the index in the current dimension.
            // It is natural to view a tensor as a multi-dimentional array.
            const int64_t current_dim_idx = current_o_idx / out_strides[ii];
            // Update the output Tensor index 'current_o_idx' by removing the offset in the current dimension.
            current_o_idx -= current_dim_idx * out_strides[ii];
            // Calculate the offset for the corresponding index position in the input Tensor and accumulate it in 'i_idx'.
            i_idx += current_dim_idx * stride[ii] * in_strides[ii];
        }
        // Assign the input Tensor element at index 'i_idx' to the output Tensor element at index 'o_idx'.
        q[o_idx] = p[i_idx];
    }
}

template <typename T>
__global__ void do_inflate_kernel(
    int ndim,
    int64_t size,
    const T* p,
    const int64_t* stride,
    const int64_t* in_strides,
    const int64_t* out_strides,
    T* q)
{
    // Perform stride operation for the specified range [begin, end) in the output Tensor.
    for (int64_t o_idx = threadIdx.x; o_idx < size; o_idx += blockDim.x) {
        int64_t i_idx = 0; // Initialize the index for the input Tensor element.
        int64_t current_o_idx = o_idx; // Calculate the index for the output Tensor element.
        bool valid = true;
        // Iterate over each dimension of the output Tensor.
        for (int ii = 0; ii < ndim; ++ii) {
            // Calculte the ratio of the current output Tensor index 'current_o_idx' in the current dimension.
            const int64_t current_dim_idx = current_o_idx / out_strides[ii];
            // Update the output Tensor index 'current_o_idx' by removing the offset in the current dimension.
            current_o_idx -= current_dim_idx * out_strides[ii];
            // Calculate the offset for the corresponding index position in the input Tensor and accumulate it in 'i_idx'.
            if (current_dim_idx % stride[ii] == 0) {
                i_idx += (current_dim_idx / stride[ii]) * in_strides[ii];
            }
            else {
                valid = false;
                break;
            }
        }
        // Assign the input Tensor element at index 'i_idx' to the output Tensor element at index 'o_idx'.
        q[o_idx] = p[i_idx] * static_cast<T>(valid ? 1.0 : 0.0);
    }
}

template <typename T>
__global__ void do_reduce_kernel(
    int64_t size,
    int64_t inner_most_dim,
    const T* p,
    T* q)
{
    for (int64_t o_idx = threadIdx.x; o_idx < size; o_idx += blockDim.x) {
        T sum = 0;
        for (int64_t i_idx = o_idx * inner_most_dim; i_idx < inner_most_dim + o_idx * inner_most_dim; i_idx++) {
            sum += p[i_idx];
        }
        q[o_idx] = sum;
    }
}

template <typename T>
static std::vector<T> compute_stride(const std::vector<T>& shape) {
    int ndims = shape.size();
    std::vector<T> strides(ndims);
    T stride = 1;

    auto it = shape.end(); // Start from the last element
    for (int ii = ndims - 1; ii >= 0; ii--) {
        it--;
        strides[ii] = stride;
        stride *= static_cast<T>(*it);
    }
    return std::move(strides);
}

template<typename T, typename Device>
void add<T, Device>::operator()(const int& num_element, const T& alpha, const T* x, const T& beta, const T* y, T* z) {
    using Type = typename GetTypeThrust<T>::type;
    auto alpha_ = *reinterpret_cast<const Type*>(&alpha);
    auto beta_ = *reinterpret_cast<const Type*>(&beta);
    const int block = (num_element + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    do_add_kernel<Type><<<block, THREADS_PER_BLOCK>>> (
        num_element, alpha_, reinterpret_cast<const Type*>(x),
        beta_, reinterpret_cast<const Type*>(y), reinterpret_cast<Type*>(z));
}

template<typename T, typename Device>
void mul<T, Device>::operator()(const int& num_element, const T& alpha, const T* x, T* y) {
    using Type = typename GetTypeThrust<T>::type;
    auto alpha_ = *reinterpret_cast<const Type*>(&alpha);
    const int block = (num_element + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    do_mul_kernel<Type><<<block, THREADS_PER_BLOCK>>> (
        num_element, alpha_,
        reinterpret_cast<const Type*>(x),  reinterpret_cast<Type*>(y));
}

template<typename T, typename Device>
void mul<T, Device>::operator()(const int& num_element, const T& alpha, const T* x, const T* y, T* z) {
    using Type = typename GetTypeThrust<T>::type;
    auto alpha_ = *reinterpret_cast<const Type*>(&alpha);
    const int block = (num_element + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    do_mul_kernel<Type><<<block, THREADS_PER_BLOCK>>> (
        num_element, alpha_,
        reinterpret_cast<const Type*>(x), reinterpret_cast<const Type*>(y), reinterpret_cast<Type*>(z));
}

template<typename T, typename Device>
void div<T, Device>::operator()(const int& num_element, const T& alpha, const T* x, const T* y, T* z) {
    using Type = typename GetTypeThrust<T>::type;
    auto alpha_ = *reinterpret_cast<const Type*>(&alpha);
    const int block = (num_element + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    do_div_kernel<Type><<<block, THREADS_PER_BLOCK>>> (
        num_element, alpha_, reinterpret_cast<const Type*>(x), reinterpret_cast<const Type*>(y), reinterpret_cast<Type*>(z));
}

template<typename T, typename Device>
void fma<T, Device>::operator()(const int& num_element, const T& alpha, const T* x, const T* y, const T& beta, const T* z, T* out) {
    using Type = typename GetTypeThrust<T>::type;
    auto alpha_ = *reinterpret_cast<const Type*>(&alpha);
    auto beta_ = *reinterpret_cast<const Type*>(&beta);
    const int block = (num_element + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    do_fma_kernel<Type><<<block, THREADS_PER_BLOCK>>> (
        num_element, alpha_, reinterpret_cast<const Type*>(x), reinterpret_cast<const Type*>(y),
        beta_, reinterpret_cast<const Type*>(z), reinterpret_cast<Type*>(out));
}

template<typename T, typename Device, bool Conjugate>
void transpose<T, Device, Conjugate>::operator()(
    const std::vector<int> &perm,
    const std::vector<int64_t> &p_shape,
    const std::vector<int64_t> &q_shape,
    const T *p,
    T *q)
{
    using Type = typename GetTypeThrust<T>::type;

    REQUIRES_OK(p_shape.size() == q_shape.size(),
                "transpose: p and q must have the same number of dimensions");
    const int ndim = static_cast<int>(p_shape.size());
    auto in_strides = compute_stride(p_shape);
    auto out_strides = compute_stride(q_shape);

    int num_elements = 1;
    for (int ii = 0; ii < ndim; ++ii) {
        num_elements *= static_cast<int>(q_shape[ii]);
    }
    num_elements = ndim ? num_elements : 0;

    Tensor t_perm(DataType::DT_INT, DeviceType::GpuDevice, {ndim});
    Tensor t_in_strides(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});
    Tensor t_out_strides(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});

    kernels::synchronize_memory<int, DEVICE_GPU, DEVICE_CPU>()(
            t_perm.data<int>(), perm.data(), perm.size());
    kernels::synchronize_memory<int64_t, DEVICE_GPU, DEVICE_CPU>()(
            t_in_strides.data<int64_t>(), in_strides.data(), in_strides.size());
    kernels::synchronize_memory<int64_t, DEVICE_GPU, DEVICE_CPU>()(
            t_out_strides.data<int64_t>(), out_strides.data(), out_strides.size());

    const Type* p_ = reinterpret_cast<const Type*>(p);
    Type* q_ = reinterpret_cast<Type*>((q));

    const int block = (num_elements + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    do_transpose_kernel<Type, Conjugate><<<block, THREADS_PER_BLOCK>>> (
        ndim, num_elements, p_, t_perm.data<int>(),
        t_in_strides.data<int64_t>(), t_out_strides.data<int64_t>(), q_);
}

template<typename T, typename Device>
void stride<T, Device>::operator()(
    const std::vector<int64_t> &stride,
    const std::vector<int64_t> &p_shape,
    const std::vector<int64_t> &q_shape,
    const T *p,
    T *q)
{
    using Type = typename GetTypeThrust<T>::type;

    REQUIRES_OK(p_shape.size() == q_shape.size(),
                "transpose: p and q must have the same number of dimensions");
    const int ndim = static_cast<int>(p_shape.size());
    auto in_strides = compute_stride(p_shape);
    auto out_strides = compute_stride(q_shape);

    int num_elements = 1;
    for (int ii = 0; ii < ndim; ++ii) {
        num_elements *= static_cast<int>(q_shape[ii]);
    }
    num_elements = ndim ? num_elements : 0;

    Tensor t_stride(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});
    Tensor t_in_strides(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});
    Tensor t_out_strides(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});

    kernels::synchronize_memory<int64_t, DEVICE_GPU, DEVICE_CPU>()(
        t_stride.data<int64_t>(), stride.data(), stride.size());
    kernels::synchronize_memory<int64_t, DEVICE_GPU, DEVICE_CPU>()(
        t_in_strides.data<int64_t>(), in_strides.data(), in_strides.size());
    kernels::synchronize_memory<int64_t, DEVICE_GPU, DEVICE_CPU>()(
        t_out_strides.data<int64_t>(), out_strides.data(), out_strides.size());

    const Type* p_ = reinterpret_cast<const Type*>(p);
    Type* q_ = reinterpret_cast<Type*>((q));

    const int block = (num_elements + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    do_stride_kernel<Type><<<block, THREADS_PER_BLOCK>>> (
        ndim, num_elements, p_, t_stride.data<int64_t>(), t_in_strides.data<int64_t>(), t_out_strides.data<int64_t>(), q_);
}


template<typename T, typename Device>
void inflate<T, Device>::operator()(
    const std::vector<int64_t> &inflate,
    const std::vector<int64_t> &p_shape,
    const std::vector<int64_t> &q_shape,
    const T *p,
    T *q)
{
    using Type = typename GetTypeThrust<T>::type;

    REQUIRES_OK(p_shape.size() == q_shape.size(),
                "transpose: p and q must have the same number of dimensions");
    const int ndim = static_cast<int>(p_shape.size());
    auto in_strides = compute_stride(p_shape);
    auto out_strides = compute_stride(q_shape);

    int num_elements = 1;
    for (int ii = 0; ii < ndim; ++ii) {
        num_elements *= static_cast<int>(q_shape[ii]);
    }
    num_elements = ndim ? num_elements : 0;

    Tensor t_stride(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});
    Tensor t_in_strides(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});
    Tensor t_out_strides(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});

    kernels::synchronize_memory<int64_t, DEVICE_GPU, DEVICE_CPU>()(
        t_stride.data<int64_t>(), inflate.data(), inflate.size());
    kernels::synchronize_memory<int64_t, DEVICE_GPU, DEVICE_CPU>()(
        t_in_strides.data<int64_t>(), in_strides.data(), in_strides.size());
    kernels::synchronize_memory<int64_t, DEVICE_GPU, DEVICE_CPU>()(
        t_out_strides.data<int64_t>(), out_strides.data(), out_strides.size());

    const Type* p_ = reinterpret_cast<const Type*>(p);
    Type* q_ = reinterpret_cast<Type*>((q));

    const int block = (num_elements + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    do_inflate_kernel<Type><<<block, THREADS_PER_BLOCK>>> (
        ndim, num_elements, p_, t_stride.data<int64_t>(), t_in_strides.data<int64_t>(), t_out_strides.data<int64_t>(), q_);
}


template<typename T, typename Device>
void reduce<T, Device>::operator()(
    const int64_t &num_element,
    const int64_t &inner_most_dim,
    const T *p,
    T *q)
{
    using Type = typename GetTypeThrust<T>::type;

    const Type* p_ = reinterpret_cast<const Type*>(p);
    Type* q_ = reinterpret_cast<Type*>((q));

    const int block = (static_cast<int>(num_element) + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    do_reduce_kernel<Type><<<block, THREADS_PER_BLOCK>>> (
        num_element, inner_most_dim, p_, q_);
}

template struct add<int, DEVICE_GPU>;
template struct add<int64_t, DEVICE_GPU>;
template struct add<float, DEVICE_GPU>;
template struct add<double, DEVICE_GPU>;
template struct add<std::complex<float>, DEVICE_GPU>;
template struct add<std::complex<double>, DEVICE_GPU>;

template struct mul<int, DEVICE_GPU>;
template struct mul<int64_t, DEVICE_GPU>;
template struct mul<float, DEVICE_GPU>;
template struct mul<double, DEVICE_GPU>;
template struct mul<std::complex<float>, DEVICE_GPU>;
template struct mul<std::complex<double>, DEVICE_GPU>;

template struct div<int, DEVICE_GPU>;
template struct div<int64_t, DEVICE_GPU>;
template struct div<float, DEVICE_GPU>;
template struct div<double, DEVICE_GPU>;
template struct div<std::complex<float>, DEVICE_GPU>;
template struct div<std::complex<double>, DEVICE_GPU>;

template struct fma<int, DEVICE_GPU>;
template struct fma<int64_t, DEVICE_GPU>;
template struct fma<float, DEVICE_GPU>;
template struct fma<double, DEVICE_GPU>;
template struct fma<std::complex<float>, DEVICE_GPU>;
template struct fma<std::complex<double>, DEVICE_GPU>;

template struct transpose<int, DEVICE_GPU>;
template struct transpose<int64_t, DEVICE_GPU>;
template struct transpose<float, DEVICE_GPU>;
template struct transpose<double, DEVICE_GPU>;
template struct transpose<std::complex<float>, DEVICE_GPU>;
template struct transpose<std::complex<double>, DEVICE_GPU>;

template struct stride<int, DEVICE_GPU>;
template struct stride<int64_t, DEVICE_GPU>;
template struct stride<float, DEVICE_GPU>;
template struct stride<double, DEVICE_GPU>;
template struct stride<std::complex<float>, DEVICE_GPU>;
template struct stride<std::complex<double>, DEVICE_GPU>;

template struct inflate<int, DEVICE_GPU>;
template struct inflate<int64_t, DEVICE_GPU>;
template struct inflate<float, DEVICE_GPU>;
template struct inflate<double, DEVICE_GPU>;
template struct inflate<std::complex<float>, DEVICE_GPU>;
template struct inflate<std::complex<double>, DEVICE_GPU>;

template struct reduce<int, DEVICE_GPU>;
template struct reduce<int64_t, DEVICE_GPU>;
template struct reduce<float, DEVICE_GPU>;
template struct reduce<double, DEVICE_GPU>;
template struct reduce<std::complex<float>, DEVICE_GPU>;
template struct reduce<std::complex<double>, DEVICE_GPU>;

} // namespace kernels
} // namespace container