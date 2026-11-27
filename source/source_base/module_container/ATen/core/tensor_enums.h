// Extracted from tensor_types.h to break the circular dependency:
// macros.h -> tensor_types.h -> cuda.h/rocm.h -> macros.h
// This file provides only enums and template specializations needed by macros.h,
// without any GPU library dependencies.
// mohan add 20260605
#ifndef ATEN_CORE_TENSOR_ENUMS_H_
#define ATEN_CORE_TENSOR_ENUMS_H_

#include <complex>

#include "source_base/module_device/types.h"

namespace container {

enum class DataType {
    DT_INVALID = 0,
    DT_FLOAT = 1,
    DT_DOUBLE = 2,
    DT_INT = 3,
    DT_INT64 = 4,
    DT_COMPLEX = 5,
    DT_COMPLEX_DOUBLE = 6,
};

struct DEVICE_CPU {};
struct DEVICE_GPU {};

enum class DeviceType {
    UnKnown = 0,
    CpuDevice = 1,
    GpuDevice = 2,
};

template <typename T>
struct DeviceTypeToEnum {
    static constexpr DeviceType value = {};
};

template <>
struct DeviceTypeToEnum<DEVICE_CPU> {
    static constexpr DeviceType value = DeviceType::CpuDevice;
};

template <>
struct DeviceTypeToEnum<DEVICE_GPU> {
    static constexpr DeviceType value = DeviceType::GpuDevice;
};

template <>
struct DeviceTypeToEnum<base_device::DEVICE_CPU> {
    static constexpr DeviceType value = DeviceType::CpuDevice;
};

template <>
struct DeviceTypeToEnum<base_device::DEVICE_GPU> {
    static constexpr DeviceType value = DeviceType::GpuDevice;
};

template <typename T>
struct DataTypeToEnum {
    static constexpr DataType value = {};
};

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

} // namespace container

#endif // ATEN_CORE_TENSOR_ENUMS_H_
