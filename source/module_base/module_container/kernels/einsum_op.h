#ifndef CONTAINER_KERNELS_EINSUM_OP_H_
#define CONTAINER_KERNELS_EINSUM_OP_H_

#include <tuple>

#include "../tensor.h"

namespace container {
namespace op {

// TODO: implement this method this week!

template <typename T, typename Device, typename... Tensors>
struct einsum_op {
    typename std::enable_if<std::is_same<
            typename std::common_type<Tensors...>::type, Tensor>::value, Tensor>::type
    einsum(const std::string &expr, Tensors *... tensors) {
        // Check the input dimension
        auto _tensors = std::make_tuple(tensors...);
        constexpr int num_tensors = sizeof...(Tensors);
        if (num_tensors > 2) {
            throw std::invalid_argument("Einstein notation only support two or less tensors");
        }

        return Tensor(DataType::DT_INT, {});
    }
};

/**
 *
 * @brief Computes the Einstein summation convention on the given tensors
 * based on the expression passed as a string.
 *
 * @tparam Tensors Variadic template parameter pack for the input tensors.
 *
 * @param expr The Einstein summation convention expression.
 * @param tensors The input tensors for the summation operation.
 * @return The resulting tensor after applying the Einstein summation convention on the input tensors.
 *
 * @throws std::invalid_argument if the expression or the input tensors are invalid.
 * @throws std::runtime_error if an error occurs while performing the summation operation.
 */
template <typename... Tensors>
    typename std::enable_if<std::is_same<
            typename std::common_type<Tensors...>::type, Tensor>::value, Tensor>::type
einsum(const std::string &expr, Tensors *... tensors) {
    // Check the input dimension
    auto _tensors = std::make_tuple(tensors...);
    constexpr int num_tensors = sizeof...(Tensors);
    if (num_tensors > 2) {
        throw std::invalid_argument("Einstein notation only support two or less tensors");
    }

    return Tensor(DataType::DT_INT, {});
}



} // namespace op
} // namespace container

#endif // CONTAINER_KERNELS_EINSUM_OP_H_