#ifndef MODULE_BASE_MACROS_H
#define MODULE_BASE_MACROS_H

#include <complex>

template <typename T>
struct PossibleComplexToReal {
    using type = T; /**< The return type based on the input type. */
};

/**
 * @brief Specialization of PossibleComplexToReal for std::complex<float>.
 *
 * This specialization sets the return type to be float when the input type is std::complex<float>.
 */
template <>
struct PossibleComplexToReal<std::complex<float>> {
    using type = float; /**< The return type specialization for std::complex<float>. */
};

/**
 * @brief Specialization of PossibleComplexToReal for std::complex<double>.
 *
 * This specialization sets the return type to be double when the input type is std::complex<double>.
 */
template <>
struct PossibleComplexToReal<std::complex<double>> {
    using type = double; /**< The return type specialization for std::complex<double>. */
};

#endif