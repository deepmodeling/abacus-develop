/**
 * @file tensor_types.h
 * @brief This file contains the definition of the DataType enum class.
 */
#ifndef ATEN_CORE_TENSOR_TYPES_H_
#define ATEN_CORE_TENSOR_TYPES_H_

#include <string>
#include <complex>
#include <iostream>
#include <stdexcept>

#include "source_base/module_device/types.h"
#include "ATen/core/tensor_enums.h"

// NOTE: Previously this file included <base/macros/cuda.h> or <base/macros/rocm.h>,
// which transitively provided GetTypeThrust, GetTypeCuda, GetTypeRocm, etc.
// Now replaced with <thrust/complex.h> which only provides thrust::complex type.
// If you need GetTypeThrust/GetTypeCuda/GetTypeRocm, include <base/macros/cuda.h>
// or <base/macros/rocm.h> directly in your .cpp file.
// mohan add 20260605
#if defined(__CUDACC__) || defined(__HIPCC__)
#include <thrust/complex.h>
#endif

namespace container {

template <typename T, int Accuracy>
static inline bool element_compare(T& a, T& b) {
    if (Accuracy <= 4) {
        return (a == b) || (std::norm(a - b) < 1e-7);
    } 
    else if (Accuracy <= 8) {
        return (a == b) || (std::norm(a - b) < 1e-15);
    } 
    else {
        return (a == b);
    }
}

/**
 * @brief Template struct to determine the return type based on the input type.
 *
 * This struct defines a template class that is used to determine the appropriate return type
 * based on the input type. By default, the return type is the same as the input type.
 *
 * @tparam T The input type for which the return type needs to be determined.
 */
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

template <typename T> 
struct PsiToContainer {
    using type = T; /**< The return type based on the input type. */
};

template <>
struct PsiToContainer<base_device::DEVICE_CPU>
{
    using type = container::DEVICE_CPU; /**< The return type specialization for std::complex<float>. */
};

template <>
struct PsiToContainer<base_device::DEVICE_GPU>
{
    using type = container::DEVICE_GPU; /**< The return type specialization for std::complex<double>. */
};

template <typename T> 
struct ContainerToPsi {
    using type = T; /**< The return type based on the input type. */
};

template <>
struct ContainerToPsi<container::DEVICE_CPU> {
    using type = base_device::DEVICE_CPU; /**< The return type specialization for std::complex<float>. */
};

template <>
struct ContainerToPsi<container::DEVICE_GPU> {
    using type = base_device::DEVICE_GPU; /**< The return type specialization for std::complex<double>. */
};

#if defined(__CUDACC__) || defined(__HIPCC__)
template <>
struct DataTypeToEnum<thrust::complex<float>> {
    static constexpr DataType value = DataType::DT_COMPLEX;
};

template <>
struct DataTypeToEnum<thrust::complex<double>> {
    static constexpr DataType value = DataType::DT_COMPLEX_DOUBLE;
};
#endif // defined(__CUDACC__) || defined(__HIPCC__)

/**
 * @brief Overloaded operator<< for the Tensor class.
 *
 * Prints the data type of the enum type DataType.
 *
 * @param os The output stream to write to.
 * @param tensor The Tensor object to print.
 *
 * @return The output stream.
 */
std::ostream& operator<<(std::ostream& os, const DataType& data_type);

/**
 * @brief Overloaded operator<< for the Tensor class.
 *
 * Prints the memory type of the enum type MemoryType.
 *
 * @param os The output stream to write to.
 * @param tensor The Tensor object to print.
 *
 * @return The output stream.
 */
std::ostream& operator<<(std::ostream& os, const DeviceType& memory_type);

} // namespace container
#endif // ATEN_CORE_TENSOR_TYPES_H_