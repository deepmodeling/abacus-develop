#include "device_helpers.h"

namespace base_device
{

// Device type specializations
template <>
AbacusDevice_t get_device_type<DEVICE_CPU>(const DEVICE_CPU* dev)
{
    return CpuDevice;
}

template <>
AbacusDevice_t get_device_type<DEVICE_GPU>(const DEVICE_GPU* dev)
{
    return GpuDevice;
}

// Precision specializations
template <>
std::string get_current_precision<float>(const float* var)
{
    return "single";
}

template <>
std::string get_current_precision<double>(const double* var)
{
    return "double";
}

template <>
std::string get_current_precision<std::complex<float>>(const std::complex<float>* var)
{
    return "single";
}

template <>
std::string get_current_precision<std::complex<double>>(const std::complex<double>* var)
{
    return "double";
}

} // end of namespace base_device
