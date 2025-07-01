#include <ATen/kernels/linalg.h>

namespace container {
namespace kernels {

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

template<typename T, typename Device>
void add<T, Device>::operator()(const int& num_element, const T& alpha, const T* x, const T& beta, const T* y, T* z) {
    // Define a lambda expression 'add_fn' to implement add operation.
    // Perform add operation for the specified range [begin, end) in the output Tensor.
    for (int o_idx = 0; o_idx < num_element; o_idx++) {
        // Assign the sum of the input Tensor elements at index 'o_idx' to the output Tensor element at index 'o_idx'.
        z[o_idx] = alpha * x[o_idx] + beta * y[o_idx];
    }
}

template<typename T, typename Device>
void mul<T, Device>::operator()(const int& num_element, const T& alpha, const T* x, T* y) {
    // Define a lambda expression 'mul_fn' to implement mul operation.
    // Perform mul operation for the specified range [begin, end) in the output Tensor.
    for (int o_idx = 0; o_idx < num_element; o_idx++) {
        // Assign the product of the input Tensor elements at index 'o_idx' to the output Tensor element at index 'o_idx'.
        y[o_idx] = alpha * x[o_idx];
    }
}

template<typename T, typename Device>
void mul<T, Device>::operator()(const int& num_element, const T& alpha, const T* x, const T* y, T* z) {
    // Define a lambda expression 'mul_fn' to implement mul operation.
    // Perform mul operation for the specified range [begin, end) in the output Tensor.
    for (int o_idx = 0; o_idx < num_element; o_idx++) {
        // Assign the product of the input Tensor elements at index 'o_idx' to the output Tensor element at index 'o_idx'.
        z[o_idx] = alpha * x[o_idx] * y[o_idx];
    }
}

template<typename T, typename Device>
void div<T, Device>::operator()(const int& num_element, const T& alpha, const T* x, const T* y, T* z) {
    // Define a lambda expression 'div_fn' to implement div operation.
    // Perform div operation for the specified range [begin, end) in the output Tensor.
    for (int o_idx = 0; o_idx < num_element; o_idx++) {
        // Assign the quotient of the input Tensor elements at index 'o_idx' to the output Tensor element at index 'o_idx'.
        z[o_idx] = alpha * x[o_idx] / y[o_idx];
    }
}

template<typename T, typename Device>
void fma<T, Device>::operator()(const int& num_element, const T& alpha, const T* x, const T* y, const T& beta, const T* z, T* out) {
    // Define a lambda expression 'fma_fn' to implement fma operation.
    // Perform fma operation for the specified range [begin, end) in the output Tensor.
    for (int o_idx = 0; o_idx < num_element; o_idx++) {
        // Assign the sum of the product of the input Tensor elements at index 'o_idx' and the corresponding coefficients to the output Tensor element at index 'o_idx'.
        out[o_idx] = alpha * x[o_idx] * y[o_idx] + beta * z[o_idx];
    }
}

template<typename T, typename Device, bool Conjugate>
void transpose<T, Device, Conjugate>::operator()(
    const std::vector<int> &perm,
    const std::vector<int64_t> &p_shape,
    const std::vector<int64_t> &q_shape,
    const T *p,
    T *q)
{
    REQUIRES_OK(p_shape.size() == q_shape.size(),
        "transpose: p and q must have the same number of dimensions");
    const int ndim = static_cast<int>(p_shape.size());
    auto in_strides = ComputeStride(p_shape);
    auto out_strides = ComputeStride(q_shape);

    int64_t num_elements = 1;
    for (int ii = 0; ii < ndim; ++ii) {
        num_elements *= q_shape[ii];
    }
    num_elements = ndim ? num_elements : 0;
    // Define a lambda expression 'transpose_fn' to implement transpose operation.
    // Perform transpose operation for the specified range [begin, end) in the output Tensor.
    for (int64_t o_idx = 0; o_idx < num_elements; o_idx++) {
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
            q[o_idx] = kernels::conj(p[i_idx]);
        } else {
            // Assign the input Tensor element at index 'i_idx' to the output Tensor element at index 'o_idx'.
            q[o_idx] = p[i_idx];
        }
    }
}


template<typename T, typename Device>
void stride<T, Device>::operator()(
    const std::vector<int64_t> &stride,
    const std::vector<int64_t> &p_shape,
    const std::vector<int64_t> &q_shape,
    const T *p,
    T *q)
{
    REQUIRES_OK(p_shape.size() == q_shape.size() ,
                "stride: p and q must match the number of dimensions");
    const int ndim = static_cast<int>(p_shape.size());
    auto in_strides = ComputeStride(p_shape);
    auto out_strides = ComputeStride(q_shape);

    int64_t num_elements = 1;
    for (int ii = 0; ii < ndim; ++ii) {
        num_elements *= q_shape[ii];
    }
    num_elements = ndim ? num_elements : 0;
    // Define a lambda expression 'stride_fn' to implement stride operation.
    // Perform stride operation for the specified range [begin, end) in the output Tensor.
    // Perform stride operation for the specified range [begin, end) in the output Tensor.
    for (int64_t o_idx = 0; o_idx < num_elements; o_idx++) {
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
            i_idx += (current_dim_idx * stride[ii]) * in_strides[ii];
        }
        // Assign the input Tensor element at index 'i_idx' to the output Tensor element at index 'o_idx'.
        q[o_idx] = p[i_idx];
    }
}


template<typename T, typename Device>
void inflate<T, Device>::operator()(
    const std::vector<int64_t> &inflate,
    const std::vector<int64_t> &p_shape,
    const std::vector<int64_t> &q_shape,
    const T *p,
    T *q)
{
    REQUIRES_OK(p_shape.size() == q_shape.size(),
                "transpose: p and q must have the same number of dimensions");
    const int ndim = static_cast<int>(p_shape.size());
    auto in_strides = ComputeStride(p_shape);
    auto out_strides = ComputeStride(q_shape);

    int64_t num_elements = 1;
    for (int ii = 0; ii < ndim; ++ii) {
        num_elements *= q_shape[ii];
    }
    num_elements = ndim ? num_elements : 0;
    // Define a lambda expression 'inflate_fn' to implement inflate operation.
    // Perform inflate operation for the specified range [begin, end) in the output Tensor.
    for (int64_t o_idx = 0; o_idx < num_elements; o_idx++) {
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
            if (current_dim_idx % inflate[ii] == 0) {
                i_idx += (current_dim_idx / inflate[ii]) * in_strides[ii];
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


template<typename T, typename Device>
void reduce<T, Device>::operator()(
    const int64_t &num_element,
    const int64_t &inner_most_dim,
    const T *p,
    T *q)
{
    // It's just so simple to implement the reduce operation.
    for (int64_t o_idx = 0; o_idx < num_element; o_idx++) {
        T sum = 0;
        for (int64_t i_idx = o_idx * inner_most_dim; i_idx < inner_most_dim + o_idx * inner_most_dim; i_idx++) {
            sum += p[i_idx];
        }
        q[o_idx] = sum;
    }
}


template struct add<int, DEVICE_CPU>;
template struct add<int64_t, DEVICE_CPU>;
template struct add<float, DEVICE_CPU>;
template struct add<double, DEVICE_CPU>;
template struct add<std::complex<float>, DEVICE_CPU>;
template struct add<std::complex<double>, DEVICE_CPU>;

template struct mul<int, DEVICE_CPU>;
template struct mul<int64_t, DEVICE_CPU>;
template struct mul<float, DEVICE_CPU>;
template struct mul<double, DEVICE_CPU>;
template struct mul<std::complex<float>, DEVICE_CPU>;
template struct mul<std::complex<double>, DEVICE_CPU>;

template struct div<int, DEVICE_CPU>;
template struct div<int64_t, DEVICE_CPU>;
template struct div<float, DEVICE_CPU>;
template struct div<double, DEVICE_CPU>;
template struct div<std::complex<float>, DEVICE_CPU>;
template struct div<std::complex<double>, DEVICE_CPU>;

template struct fma<int, DEVICE_CPU>;
template struct fma<int64_t, DEVICE_CPU>;
template struct fma<float, DEVICE_CPU>;
template struct fma<double, DEVICE_CPU>;
template struct fma<std::complex<float>, DEVICE_CPU>;
template struct fma<std::complex<double>, DEVICE_CPU>;

template struct transpose<int, DEVICE_CPU>;
template struct transpose<int64_t, DEVICE_CPU>;
template struct transpose<float, DEVICE_CPU>;
template struct transpose<double, DEVICE_CPU>;
template struct transpose<std::complex<float>, DEVICE_CPU>;
template struct transpose<std::complex<double>, DEVICE_CPU>;

template struct stride<int, DEVICE_CPU>;
template struct stride<int64_t, DEVICE_CPU>;
template struct stride<float, DEVICE_CPU>;
template struct stride<double, DEVICE_CPU>;
template struct stride<std::complex<float>, DEVICE_CPU>;
template struct stride<std::complex<double>, DEVICE_CPU>;

template struct inflate<int, DEVICE_CPU>;
template struct inflate<int64_t, DEVICE_CPU>;
template struct inflate<float, DEVICE_CPU>;
template struct inflate<double, DEVICE_CPU>;
template struct inflate<std::complex<float>, DEVICE_CPU>;
template struct inflate<std::complex<double>, DEVICE_CPU>;

template struct reduce<int, DEVICE_CPU>;
template struct reduce<int64_t, DEVICE_CPU>;
template struct reduce<float, DEVICE_CPU>;
template struct reduce<double, DEVICE_CPU>;
template struct reduce<std::complex<float>, DEVICE_CPU>;
template struct reduce<std::complex<double>, DEVICE_CPU>;

} // namespace kernels
} // namespace container