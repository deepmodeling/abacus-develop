#include "diag_const_nums.h"

#include <complex>

template class const_nums<double>;
template class const_nums<float>;
template class const_nums<std::complex<double>>;
template class const_nums<std::complex<float>>;

// Specialize templates to support double types
template <>
const_nums<double>::const_nums()
{
    base_device::memory::resize_memory_op<double, base_device::DEVICE_CPU>()(this->zero, 1);
    this->zero[0] = 0.0;
    base_device::memory::resize_memory_op<double, base_device::DEVICE_CPU>()(this->one, 1);
    this->one[0] = 1.0;
    base_device::memory::resize_memory_op<double, base_device::DEVICE_CPU>()(this->neg_one, 1);
    this->neg_one[0] = -1.0;
}

// Specialize templates to support double types
template <>
const_nums<float>::const_nums()
{
    base_device::memory::resize_memory_op<float, base_device::DEVICE_CPU>()(this->zero, 1);
    this->zero[0] = 0.0;
    base_device::memory::resize_memory_op<float, base_device::DEVICE_CPU>()(this->one, 1);
    this->one[0] = 1.0;
    base_device::memory::resize_memory_op<float, base_device::DEVICE_CPU>()(this->neg_one, 1);
    this->neg_one[0] = -1.0;
}

// Specialized templates to support std:: complex<double>types
template <>
const_nums<std::complex<double>>::const_nums()
{
    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_CPU>()(this->zero, 1);
    this->zero[0] = std::complex<double>(0.0, 0.0);
    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_CPU>()(this->one, 1);
    this->one[0] = std::complex<double>(1.0, 0.0);
    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_CPU>()(this->neg_one, 1);
    this->neg_one[0] = std::complex<double>(-1.0, 0.0);
}

// Specialized templates to support std:: complex<float>types
template <>
const_nums<std::complex<float>>::const_nums()
{
    base_device::memory::resize_memory_op<std::complex<float>, base_device::DEVICE_CPU>()(this->zero, 1);
    this->zero[0] = std::complex<float>(0.0, 0.0);
    base_device::memory::resize_memory_op<std::complex<float>, base_device::DEVICE_CPU>()(this->one, 1);
    this->one[0] = std::complex<float>(1.0, 0.0);
    base_device::memory::resize_memory_op<std::complex<float>, base_device::DEVICE_CPU>()(this->neg_one, 1);
    this->neg_one[0] = std::complex<float>(-1.0, 0.0);
}