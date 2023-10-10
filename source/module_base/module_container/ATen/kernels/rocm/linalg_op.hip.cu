#include <ATen/kernels/linalg_op.h>
#include <ATen/core/tensor.h>
#include <ATen/core/tensor_types.h>
#include <base/macros/rocm.h>

#include <hip/hip_runtime.h>
#include <thrust/complex.h>

namespace container {
namespace op {

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

template <typename T, bool Conjugate>
__global__ void do_transpose_kernel(
    int ndim,
    int64_t size,
    const T* p,
    const int* perm,
    const int64_t* in_strides,
    const int64_t* out_strides,
    T* q)
{
    for (int64_t o_idx = threadIdx.x; o_idx < size; o_idx += blockDim.x) {
        int64_t i_idx = 0; // Initialize the index for the input Tensor element.
        int64_t current_o_idx = o_idx; // Calculate the index for the output Tensor element.

        // Iterate over each dimension of the output Tensor.
        for (int ii = 0; ii < ndim; ++ii) {
            // Calculate the current_dim_idx of the current output Tensor index 'current_o_idx' in the current dimension.
            const int64_t current_dim_idx = current_o_idx / out_strides[ii];
            // Update the output Tensor index 'current_o_idx' by removing the offset in the current dimension.
            current_o_idx -= current_dim_idx * out_strides[ii];
            // Calculate the offset for the corresponding index position in the input Tensor and accumulate it in 'i_idx'.
            i_idx += current_dim_idx * in_strides[perm[ii]];
        }
        // Check if conjugation is needed.
        if (Conjugate) {
            // Assign the conjugate value of the input Tensor element at index 'i_idx' to the output Tensor element at index 'o_idx'.
            q[o_idx] = op::conj(p[i_idx]);
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
    int ndim,
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

template <typename T, bool Conjugate>
struct transpose_op<T, DEVICE_GPU, Conjugate> {
    using Type = typename GetTypeThrust<T>::type;
    void operator()(
        const Tensor& input,
        const std::vector<int>& perm,
        Tensor& output)
    {
        const int ndim = input.shape().ndim();
        const int64_t output_size = output.NumElements();
        auto in_strides = compute_stride<int64_t>(input.shape().dims());
        auto out_strides = compute_stride<int64_t>(output.shape().dims());

        Tensor t_perm(DataType::DT_INT, DeviceType::GpuDevice, {ndim});
        Tensor t_in_strides(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});
        Tensor t_out_strides(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});

        op::synchronize_memory_op<int, DEVICE_GPU, DEVICE_CPU>()(
                t_perm.data<int>(), perm.data(), perm.size());
        op::synchronize_memory_op<int64_t, DEVICE_GPU, DEVICE_CPU>()(
                t_in_strides.data<int64_t>(), in_strides.data(), in_strides.size());
        op::synchronize_memory_op<int64_t, DEVICE_GPU, DEVICE_CPU>()(
                t_out_strides.data<int64_t>(), out_strides.data(), out_strides.size());

        const Type* p = reinterpret_cast<const Type*>(input.data());
        Type* q = reinterpret_cast<Type*>((output.data()));

        const int block = (output_size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        do_transpose_kernel<Type, Conjugate><<<block, THREADS_PER_BLOCK>>> (
            ndim, output_size, p, t_perm.data<int>(), 
            t_in_strides.data<int64_t>(), t_out_strides.data<int64_t>(), q);
    }
};

template <typename T>
struct stride_op<T, DEVICE_GPU> {
    using Type = typename GetTypeThrust<T>::type;
    void operator()(
        const Tensor& input,
        const TensorShape& stride,
        Tensor& output)
    {
        const int ndim = input.shape().ndim();
        const int64_t output_size = output.NumElements();
        auto in_strides = compute_stride<int64_t>(input.shape().dims());
        auto out_strides = compute_stride<int64_t>(output.shape().dims());

        Tensor t_stride(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});
        Tensor t_in_strides(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});
        Tensor t_out_strides(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});

        op::synchronize_memory_op<int64_t, DEVICE_GPU, DEVICE_CPU>()(
                t_stride.data<int64_t>(), stride.dims().data(), t_stride.NumElements());
        op::synchronize_memory_op<int64_t, DEVICE_GPU, DEVICE_CPU>()(
                t_in_strides.data<int64_t>(), in_strides.data(), in_strides.size());
        op::synchronize_memory_op<int64_t, DEVICE_GPU, DEVICE_CPU>()(
                t_out_strides.data<int64_t>(), out_strides.data(), out_strides.size());

        const Type* p = reinterpret_cast<const Type*>(input.data());
        Type* q = reinterpret_cast<Type*>((output.data()));

        const int block = (output_size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        do_stride_kernel<Type><<<block, THREADS_PER_BLOCK>>> (
            ndim, output_size, p, t_stride.data<int64_t>(), t_in_strides.data<int64_t>(), t_out_strides.data<int64_t>(), q);
    }
};

template <typename T>
struct inflate_op<T, DEVICE_GPU> {
    using Type = typename GetTypeThrust<T>::type;
    void operator()(
        const Tensor& input,
        const TensorShape& stride,
        Tensor& output)
    {
        const int ndim = input.shape().ndim();
        const int64_t output_size = output.NumElements();
        auto in_strides = compute_stride<int64_t>(input.shape().dims());
        auto out_strides = compute_stride<int64_t>(output.shape().dims());

        Tensor t_stride(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});
        Tensor t_in_strides(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});
        Tensor t_out_strides(DataType::DT_INT64, DeviceType::GpuDevice, {ndim});

        op::synchronize_memory_op<int64_t, DEVICE_GPU, DEVICE_CPU>()(
                t_stride.data<int64_t>(), stride.dims().data(), t_stride.NumElements());
        op::synchronize_memory_op<int64_t, DEVICE_GPU, DEVICE_CPU>()(
                t_in_strides.data<int64_t>(), in_strides.data(), in_strides.size());
        op::synchronize_memory_op<int64_t, DEVICE_GPU, DEVICE_CPU>()(
                t_out_strides.data<int64_t>(), out_strides.data(), out_strides.size());

        const Type* p = reinterpret_cast<const Type*>(input.data());
        Type* q = reinterpret_cast<Type*>((output.data()));

        const int block = (output_size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        do_inflate_kernel<Type><<<block, THREADS_PER_BLOCK>>> (
            ndim, output_size, p, t_stride.data<int64_t>(), t_in_strides.data<int64_t>(), t_out_strides.data<int64_t>(), q);
    }
};

template <typename T>
struct reduce_op<T, DEVICE_GPU> {
    using Type = typename GetTypeThrust<T>::type;
    void operator()(
        const Tensor& input,
        const int64_t& inner_most_dim,
        Tensor& output)
    {
        const int ndim = input.shape().ndim();
        const int64_t output_size = output.NumElements();

        const Type* p = reinterpret_cast<const Type*>(input.data());
        Type* q = reinterpret_cast<Type*>((output.data()));

        const int block = (output_size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        do_reduce_kernel<Type><<<block, THREADS_PER_BLOCK>>> (
            ndim, output_size, inner_most_dim, p, q);
    }
};

template struct transpose_op<int, DEVICE_GPU>;
template struct transpose_op<int64_t, DEVICE_GPU>;
template struct transpose_op<float, DEVICE_GPU>;
template struct transpose_op<double, DEVICE_GPU>;
template struct transpose_op<std::complex<float>, DEVICE_GPU>;
template struct transpose_op<std::complex<double>, DEVICE_GPU>;

template struct stride_op<int, DEVICE_GPU>;
template struct stride_op<int64_t, DEVICE_GPU>;
template struct stride_op<float, DEVICE_GPU>;
template struct stride_op<double, DEVICE_GPU>;
template struct stride_op<std::complex<float>, DEVICE_GPU>;
template struct stride_op<std::complex<double>, DEVICE_GPU>;

template struct inflate_op<int, DEVICE_GPU>;
template struct inflate_op<int64_t, DEVICE_GPU>;
template struct inflate_op<float, DEVICE_GPU>;
template struct inflate_op<double, DEVICE_GPU>;
template struct inflate_op<std::complex<float>, DEVICE_GPU>;
template struct inflate_op<std::complex<double>, DEVICE_GPU>;

template struct reduce_op<int, DEVICE_GPU>;
template struct reduce_op<int64_t, DEVICE_GPU>;
template struct reduce_op<float, DEVICE_GPU>;
template struct reduce_op<double, DEVICE_GPU>;
template struct reduce_op<std::complex<float>, DEVICE_GPU>;
template struct reduce_op<std::complex<double>, DEVICE_GPU>;

} // namespace op
} // namespace container