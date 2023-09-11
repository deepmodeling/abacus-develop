/**
 * @file tensor_types.h
 * @brief This file contains the definition of the DataType enum class.
 */
#ifndef ATEN_CORE_TENSOR_TYPES_H_
#define ATEN_CORE_TENSOR_TYPES_H_

#include <regex>
#include <string>
#include <vector>
#include <numeric>
#include <complex>
#include <utility>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <unordered_map>
#include <initializer_list>

#include <module_psi/kernels/types.h>

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
@brief Enumeration of data types for tensors.
The DataType enum lists the supported data types for tensors. Each data type
is identified by a unique value. The DT_INVALID value is reserved for invalid
data types.
*/
enum class DataType {
    DT_INVALID = 0, ///< Invalid data type */
    DT_FLOAT = 1, ///< Single-precision floating point */
    DT_DOUBLE = 2, ///< Double-precision floating point */
    DT_INT = 3, ///< 32-bit integer */
    DT_INT64 = 4, ///< 64-bit integer */
    DT_COMPLEX = 5, ///< 32-bit complex */
    DT_COMPLEX_DOUBLE = 6, /**< 64-bit complex */
// ... other data types
};

/**
 *@struct DEVICE_CPU, DEVICE_GPU
 *@brief A tag type for identifying CPU and GPU devices.
*/
struct DEVICE_CPU;
struct DEVICE_GPU;

struct DEVICE_CPU {};
struct DEVICE_GPU {};
/**
 * @brief The type of memory used by an allocator.
 */
enum class DeviceType {
    UnKnown = 0,  ///< Memory type is unknown.
    CpuDevice = 1,     ///< Memory type is CPU.
    GpuDevice = 2,     ///< Memory type is GPU(CUDA or ROCm).
};

/**
 * @brief Template struct to determine the return type based on the input type.
 *
 * This struct defines a template class that is used to determine the appropriate return type
 * based on the input type. By default, the return type is the same as the input type.
 *
 * @tparam T The input type for which the return type needs to be determined.
 */
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

template <typename T> 
struct PsiToContainer {
    using type = T; /**< The return type based on the input type. */
};

template <>
struct PsiToContainer<psi::DEVICE_CPU> {
    using type = container::DEVICE_CPU; /**< The return type specialization for std::complex<float>. */
};

template <>
struct PsiToContainer<psi::DEVICE_GPU> {
    using type = container::DEVICE_GPU; /**< The return type specialization for std::complex<double>. */
};



/**
 * @brief Template struct for mapping a Device Type to its corresponding enum value.
 *
 * @param T The DataType to map to its enum value.
 *
 * @return The enumeration value corresponding to the data type.
 *  This method uses template specialization to map each supported data type to its
 *  corresponding enumeration value. If the template argument T is not a supported
 *  data type, this method will cause a compile-time error.
 *  Example usage:
 *      DataTypeToEnum<float>::value; // Returns DataType::DT_FLOAT
 */
template <typename T>
struct DeviceTypeToEnum {
    static constexpr DeviceType value = {};
};
// Specializations of DeviceTypeToEnum for supported devices.
template <>
struct DeviceTypeToEnum<DEVICE_CPU> {
    static constexpr DeviceType value = DeviceType::CpuDevice;
};
template <>
struct DeviceTypeToEnum<DEVICE_GPU> {
    static constexpr DeviceType value = DeviceType::GpuDevice;
};
template <>
struct DeviceTypeToEnum<psi::DEVICE_CPU> {
    static constexpr DeviceType value = DeviceType::CpuDevice;
};
template <>
struct DeviceTypeToEnum<psi::DEVICE_GPU> {
    static constexpr DeviceType value = DeviceType::GpuDevice;
};

/**
 * @brief Template struct for mapping a DataType to its corresponding enum value.
 *
 * @param T The DataType to map to its enum value.
 *
 * @return The enumeration value corresponding to the data type.
 *  This method uses template specialization to map each supported data type to its
 *  corresponding enumeration value. If the template argument T is not a supported
 *  data type, this method will cause a compile-time error.
 *  Example usage:
 *      DataTypeToEnum<float>::value; // Returns DataType::DT_FLOAT
 */
template <typename T>
struct DataTypeToEnum {
    static constexpr DataType value = {};
};
// Specializations of DataTypeToEnum for supported types.
template <>
struct DataTypeToEnum<int> {
    static constexpr DataType value = DataType::DT_INT;
};
template <>
struct DataTypeToEnum<float> {
    static constexpr DataType value = DataType::DT_FLOAT;
};
template <>
struct DataTypeToEnum<double> {
    static constexpr DataType value = DataType::DT_DOUBLE;
};
template <>
struct DataTypeToEnum<int64_t> {
    static constexpr DataType value = DataType::DT_INT64;
};
template <>
struct DataTypeToEnum<std::complex<float>> {
    static constexpr DataType value = DataType::DT_COMPLEX;
};
template <>
struct DataTypeToEnum<std::complex<double>> {
    static constexpr DataType value = DataType::DT_COMPLEX_DOUBLE;
};

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