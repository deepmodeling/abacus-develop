#ifndef MODULE_BASE_MACROS_H
#define MODULE_BASE_MACROS_H

#include <complex>

template <typename T>
struct GetTypeReal {
    using type = T; /**< The return type based on the input type. */
};

/**
 * @brief Specialization of GetTypeReal for std::complex<float>.
 *
 * This specialization sets the return type to be float when the input type is std::complex<float>.
 */
template <>
struct GetTypeReal<std::complex<float>> {
    using type = float; /**< The return type specialization for std::complex<float>. */
};

/**
 * @brief Specialization of GetTypeReal for std::complex<double>.
 *
 * This specialization sets the return type to be double when the input type is std::complex<double>.
 */
template <>
struct GetTypeReal<std::complex<double>> {
    using type = double; /**< The return type specialization for std::complex<double>. */
};

#endif