#ifndef ATEN_OPS_LINALG_H_
#define ATEN_OPS_LINALG_H_

#include <ATen/core/tensor.h>

#include <ATen/ops/einsum_op.h>

namespace container {
namespace op {

/**
 * @brief A functor to perform add operation on a Tensor.
 *
 * This functor adds two Tensors element-wise, resulting in a new Tensor with the same
 * shape as the input Tensors.
 */
struct add_op {
    /**
     * @brief Perform add operation on the input Tensors.
     *
     * This function adds two Tensors element-wise, resulting in a new Tensor with the same
     * shape as the input Tensors.
     *
     * @param x The first input Tensor.
     * @param y The second input Tensor.
     * @param z The output Tensor that will hold the result of the add operation.
     *          It must have the same shape as the input Tensors.
     */
    void operator()(
        const Tensor& x,
        const Tensor& y,
        Tensor& z);

    template <typename T>
    void operator()(
        const T& alpha,
        const Tensor& x,
        const T& beta,
        const Tensor& y,
        Tensor& z);
};

struct mul_op {
    // z = x * y
    void operator()(
        const Tensor& x,
        const Tensor& y,
        Tensor& z);

    // y = alpha * x
    template <typename T>
    void operator()(
        const T& alpha,
        const Tensor& x,
        Tensor& y);
};

struct div_op {
    // z = x / y
    void operator()(
        const Tensor& x,
        const Tensor& y,
        Tensor& z);
};

template <bool Conjugate = false>
struct transpose_op {
    /**
     * @brief Perform the transpose operation on the input tensor.
     *
     * @param input The input tensor to be transposed.
     * @param permutation A vector representing the desired permutation of the
     *                    dimensions. The elements of the vector indicate the new
     *                    order of the dimensions, and the size of the vector must
     *                    be equal to the rank of the input tensor.
     * @param output The output tensor that will store the result of the transpose
     *               operation. The output tensor should be pre-allocated and have
     *               the correct shape for the transpose result.
     *
     * @return True if the transpose operation is successful; otherwise, false.
     *
     * @note The size and shape of the output tensor must match the new dimensions
     *       specified by the permutation vector. If the permutation is not valid or
     *       the output tensor is not pre-allocated with the correct shape, the
     *       function will return false.
     */
    void operator()(
        const Tensor& input,
        const std::vector<int>& permutation,
        Tensor& output);
};

/**
 * @brief A functor to perform stride operation on a Tensor.
 *
 * This functor applies a stride operation on the input Tensor, resulting in a new
 * Tensor with data elements selected based on the specified stride. The stride
 * specifies how elements are skipped in each dimension during the selection process.
 *
 * @tparam T The data type of the Tensor.
 * @tparam Device The execution device (e.g., CPU or GPU).
 */
struct stride_op {
    /**
     * @brief Perform stride operation on the input Tensor.
     *
     * This function applies the stride operation on the input Tensor, generating a new
     * Tensor with data elements selected based on the specified stride.
     *
     * @param input The input Tensor on which the stride operation is performed.
     * @param stride The stride specification, represented as a TensorShape.
     *               The stride specifies how elements are skipped in each dimension during
     *               the selection process.
     * @param output The output Tensor that will hold the result of the stride operation.
     *               It must have the appropriate size to store the selected elements.
     */
    void operator()(
        const Tensor& input,
        const std::vector<int64_t>& stride,
        Tensor& output);
};

/**
 * @brief A functor for inflating a tensor.
 *
 * This struct defines a functor that can be used to inflate a tensor using the specified stride.
 */
struct inflate_op {
    /**
     * @brief Inflate the input tensor.
     *
     * This function inflates the input tensor using the given stride and stores the result in the output tensor.
     *
     * @param input The input tensor to be inflated.
     * @param stride The stride to use for inflation.
     * @param output The output tensor where the inflated data will be stored.
     */
    void operator()(
        const Tensor& input,
        const std::vector<int64_t>& stride,
        Tensor& output);
};


struct reduce_op {
    void operator()(
        const Tensor& input,
        const int64_t& inner_most_dim,
        Tensor& output);
};

} // namespace op
} // namespace container

ct::Tensor   operator+(const ct::Tensor& self, const ct::Tensor& other);
ct::Tensor   operator-(const ct::Tensor& self, const ct::Tensor& other);
ct::Tensor   operator*(const ct::Tensor& self, const ct::Tensor& other);
ct::Tensor   operator/(const ct::Tensor& self, const ct::Tensor& other);
ct::Tensor& operator+=(ct::Tensor& self, const ct::Tensor& other);
ct::Tensor& operator-=(ct::Tensor& self, const ct::Tensor& other);
ct::Tensor& operator*=(ct::Tensor& self, const ct::Tensor& other);
ct::Tensor& operator/=(ct::Tensor& self, const ct::Tensor& other);

#endif // ATEN_OPS_LINALG_H_