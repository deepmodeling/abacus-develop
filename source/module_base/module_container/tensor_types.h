/**
 * @file tensor_types.h
 * @brief This file contains the definition of the DataType enum class.
 */
#ifndef CONTAINER_TENSOR_TYPES_H_
#define CONTAINER_TENSOR_TYPES_H_

#include <complex>
#include <iostream>

namespace container {

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

/**
 * @brief The type of memory used by an allocator.
 */
enum class DeviceType {
    UnKnown = 0,  ///< Memory type is unknown.
    CpuDevice = 1,      ///< Memory type is CPU.
    GpuDevice = 2,     ///< Memory type is GPU(CUDA or ROCm).
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

// The macro TEMPLATE_1() expands to a switch statement conditioned on
// TYPE_ENUM. Each case expands the STMTS after a typedef for T.
#define SINGLE_ARG(...) __VA_ARGS__

#define CASE_2(TYPE, DEVICE, STMTS)               \
  case (int(DataTypeToEnum<TYPE>::value) * 10 +   \
        int(DeviceTypeToEnum<DEVICE>::value)): {  \
    typedef TYPE T_;                              \
    typedef DEVICE DEVICE_;                       \
    STMTS;                                        \
    break;                                        \
  }

#define CASES_ALL_WITH_DEFAULT_2(TYPE_ENUM, DEVICE_ENUM, STMTS, DEFAULT) \
  switch (int(TYPE_ENUM) * 10 + int(DEVICE_ENUM)) {                      \
    CASE_2(float, DEVICE_CPU, SINGLE_ARG(STMTS))                         \
    CASE_2(float, DEVICE_GPU, SINGLE_ARG(STMTS))                         \
    CASE_2(double, DEVICE_CPU, SINGLE_ARG(STMTS))                        \
    CASE_2(double, DEVICE_GPU, SINGLE_ARG(STMTS))                        \
    CASE_2(int, DEVICE_CPU, SINGLE_ARG(STMTS))                           \
    CASE_2(int, DEVICE_GPU, SINGLE_ARG(STMTS))                           \
    CASE_2(int64_t, DEVICE_CPU, SINGLE_ARG(STMTS))                       \
    CASE_2(int64_t, DEVICE_GPU, SINGLE_ARG(STMTS))                       \
    CASE_2(std::complex<float>, DEVICE_CPU, SINGLE_ARG(STMTS))           \
    CASE_2(std::complex<float>, DEVICE_GPU, SINGLE_ARG(STMTS))           \
    CASE_2(std::complex<double>, DEVICE_CPU, SINGLE_ARG(STMTS))          \
    CASE_2(std::complex<double>, DEVICE_GPU, SINGLE_ARG(STMTS))          \
    default:                                                             \
      DEFAULT;                                                           \
      break;                                                             \
  }

#define CASES_CZ_WITH_DEFAULT_2(TYPE_ENUM, DEVICE_ENUM, STMTS, DEFAULT)  \
  switch (int(TYPE_ENUM) * 10 + int(DEVICE_ENUM)) {                      \
    CASE_2(std::complex<float>, DEVICE_CPU, SINGLE_ARG(STMTS))           \
    CASE_2(std::complex<float>, DEVICE_GPU, SINGLE_ARG(STMTS))           \
    CASE_2(std::complex<double>, DEVICE_CPU, SINGLE_ARG(STMTS))          \
    CASE_2(std::complex<double>, DEVICE_GPU, SINGLE_ARG(STMTS))          \
    default:                                                             \
      DEFAULT;                                                           \
      break;                                                             \
  }

// TODO: add a log solution to container.
#define TEMPLATE_ALL_2(TYPE_ENUM, DEVICE_ENUM, ...)                      \
  CASES_ALL_WITH_DEFAULT_2(TYPE_ENUM, DEVICE_ENUM, (__VA_ARGS__),        \
                           std::cerr << "Unexpected type: " << TYPE_ENUM; exit(EXIT_FAILURE));

#define TEMPLATE_CZ_2(TYPE_ENUM, DEVICE_ENUM, ...)                       \
  CASES_CZ_WITH_DEFAULT_2(TYPE_ENUM, DEVICE_ENUM, (__VA_ARGS__),         \
                           std::cerr << "Unexpected type: " << TYPE_ENUM; exit(EXIT_FAILURE));

} // namespace container
#endif // CONTAINER_TENSOR_TYPES_H_