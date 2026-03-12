#ifndef DEVICE_HELPERS_H_
#define DEVICE_HELPERS_H_

/**
 * @file device_helpers.h
 * @brief Type trait templates for device and precision detection.
 *
 * This header provides template declarations for:
 * - get_device_type<Device>() - returns device type enum
 * - get_current_precision<T>() - returns "single" or "double"
 */

#include "types.h"
#include <complex>
#include <string>
#include <type_traits>

namespace base_device
{

// Forward declaration
class DeviceContext;

/**
 * @brief Get the device type enum from DeviceContext (runtime version).
 * @param ctx Pointer to DeviceContext
 * @return AbacusDevice_t enum value
 */
inline AbacusDevice_t get_device_type(const DeviceContext* ctx);

/**
 * @brief Get the device type enum for a given device type (compile-time version).
 * @tparam Device The device type (DEVICE_CPU or DEVICE_GPU)
 * @param dev Pointer to device (used for template deduction)
 * @return AbacusDevice_t enum value
 */
template <typename Device>
AbacusDevice_t get_device_type(const Device* dev)
{
    if (std::is_same<Device, DEVICE_CPU>::value) return CpuDevice;
    else if (std::is_same<Device, DEVICE_GPU>::value) return GpuDevice;
    else if (std::is_same<Device, DEVICE_DSP>::value) return DspDevice;
    else return UnKnown;
}

/**
 * @brief Get the precision string for a given numeric type.
 * @tparam T The numeric type (float, double, std::complex<float>, std::complex<double>)
 * @param var Pointer to variable (used for template deduction)
 * @return "single" or "double"
 */
template <typename T>
std::string get_current_precision(const T* var);

// Template specialization declarations
template <>
std::string get_current_precision<float>(const float* var);

template <>
std::string get_current_precision<double>(const double* var);

template <>
std::string get_current_precision<std::complex<float>>(const std::complex<float>* var);

template <>
std::string get_current_precision<std::complex<double>>(const std::complex<double>* var);

} // end of namespace base_device

#endif // DEVICE_HELPERS_H_
