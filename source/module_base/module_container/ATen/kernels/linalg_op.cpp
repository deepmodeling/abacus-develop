#include <ATen/kernels/linalg_op.h>

namespace container {
namespace op {

template <typename T>
static inline T conj(T& in) {
    return in;
}

template <typename T>
static inline std::complex<T> conj(std::complex<T>& in) {
    return std::conj(in);
}

template <typename T>
static std::vector<T> ComputeStride(const std::vector<T>& shape) {
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

template <typename T, typename Device, bool Conjugate>
void transpose_op<T, Device, Conjugate>::operator()(
    const Tensor& input,
    const std::vector<int>& perm,
    Tensor& output)
{
    const int ndim = input.shape().ndim();
    auto in_strides = ComputeStride(input.shape().dims());
    auto out_strides = ComputeStride(output.shape().dims());

    const T* p = reinterpret_cast<const T*>(input.data());
    T* q = reinterpret_cast<T*>((output.data()));

    // Define a lambda expression 'transpose_fn' to implement transpose operation.
    // Perform transpose operation for the specified range [begin, end) in the output Tensor.
    for (int64_t o_idx = 0; o_idx < output.shape().NumElements(); o_idx++) {
        int64_t i_idx = 0; // Initialize the index for the input Tensor element.
        int64_t t = o_idx; // Calculate the index for the output Tensor element.

        // Iterate over each dimension of the output Tensor.
        for (int ii = 0; ii < ndim; ++ii) {
            // Calculate the ratio of the current output Tensor index 't' in the current dimension.
            const int64_t ratio = t / out_strides[ii];
            // Update the output Tensor index 't' by removing the offset in the current dimension.
            t -= ratio * out_strides[ii];
            // Calculate the offset for the corresponding index position in the input Tensor and accumulate it in 'i_idx'.
            i_idx += ratio * in_strides[perm[ii]];
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

template <typename T, typename Device>
void stride_op<T, Device>::operator()(
    const Tensor& input,
    const TensorShape& stride,
    Tensor& output)
{
    const int ndim = input.shape().ndim();
    auto in_strides = ComputeStride(input.shape().dims());
    auto out_strides = ComputeStride(output.shape().dims());

    const T* p = input.data<T>();
    T* q = output.data<T>();

    // Perform stride operation for the specified range [begin, end) in the output Tensor.
    for (int64_t o_idx = 0; o_idx < output.NumElements(); o_idx++) {
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
            i_idx += (current_dim_idx * stride.dim_size(ii)) * in_strides[ii];
        }
        // Assign the input Tensor element at index 'i_idx' to the output Tensor element at index 'o_idx'.
        q[o_idx] = p[i_idx];
    }
}


template <typename T, typename Device>
void inflate_op<T, Device>::operator()(
    const Tensor& input,
    const TensorShape& stride,
    Tensor& output)
{
    const int ndim = input.shape().ndim();
    auto in_strides = ComputeStride(input.shape().dims());
    auto out_strides = ComputeStride(output.shape().dims());

    const T* p = input.data<T>();
    T* q = output.data<T>();

    // Perform stride operation for the specified range [begin, end) in the output Tensor.
    for (int64_t o_idx = 0; o_idx < output.NumElements(); o_idx++) {
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
            if (current_dim_idx % stride.dim_size(ii) == 0) {
                i_idx += (current_dim_idx / stride.dim_size(ii)) * in_strides[ii];
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

template <typename T, typename Device>
void reduce_op<T, Device>::operator()(
    const Tensor& input,
    const int64_t& inner_most_dim,
    Tensor& output)
{
    const T* p = reinterpret_cast<const T*>(input.data());
    T* q = reinterpret_cast<T*>((output.data()));

    // It's just so simple to implement the reduce operation.
    for (int64_t o_idx = 0; o_idx < output.NumElements(); o_idx++) {
        T sum = 0;
        for (int64_t i_idx = o_idx * inner_most_dim; i_idx < inner_most_dim + o_idx * inner_most_dim; i_idx++) {
            sum += p[i_idx];
        }
        q[o_idx] = sum;
    }
}

template struct transpose_op<int, DEVICE_CPU>;
template struct transpose_op<int64_t, DEVICE_CPU>;
template struct transpose_op<float, DEVICE_CPU>;
template struct transpose_op<double, DEVICE_CPU>;
template struct transpose_op<std::complex<float>, DEVICE_CPU>;
template struct transpose_op<std::complex<double>, DEVICE_CPU>;

template struct stride_op<int, DEVICE_CPU>;
template struct stride_op<int64_t, DEVICE_CPU>;
template struct stride_op<float, DEVICE_CPU>;
template struct stride_op<double, DEVICE_CPU>;
template struct stride_op<std::complex<float>, DEVICE_CPU>;
template struct stride_op<std::complex<double>, DEVICE_CPU>;

template struct inflate_op<int, DEVICE_CPU>;
template struct inflate_op<int64_t, DEVICE_CPU>;
template struct inflate_op<float, DEVICE_CPU>;
template struct inflate_op<double, DEVICE_CPU>;
template struct inflate_op<std::complex<float>, DEVICE_CPU>;
template struct inflate_op<std::complex<double>, DEVICE_CPU>;

template struct reduce_op<int, DEVICE_CPU>;
template struct reduce_op<int64_t, DEVICE_CPU>;
template struct reduce_op<float, DEVICE_CPU>;
template struct reduce_op<double, DEVICE_CPU>;
template struct reduce_op<std::complex<float>, DEVICE_CPU>;
template struct reduce_op<std::complex<double>, DEVICE_CPU>;

} // namespace op
} // namespace container