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
    //base_device::memory::resize_memory_op<double, base_device::DEVICE_CPU>()(
    //                    this->cpu_ctx, this->zero, 1);
    this->zero = (double*)malloc(sizeof(double));
    this->zero[0] = 0.0;
    //base_device::memory::resize_memory_op<double, base_device::DEVICE_CPU>()(
    //                    this->cpu_ctx, this->one, 1);
    this->one = (double*)malloc(sizeof(double));
    this->one[0] = 1.0;
    //base_device::memory::resize_memory_op<double, base_device::DEVICE_CPU>()(
    //                    this->cpu_ctx, this->neg_one, 1);
    this->neg_one = (double*)malloc(sizeof(double));
    this->neg_one[0] = -1.0;
}

// Specialize templates to support double types
template <>
const_nums<float>::const_nums()
{
    //base_device::memory::resize_memory_op<float, base_device::DEVICE_CPU>()(
    //                    this->cpu_ctx, this->zero, 1);
    this->zero = (float*)malloc(sizeof(float));
    this->zero[0] = 0.0;
    //base_device::memory::resize_memory_op<float, base_device::DEVICE_CPU>()(
    //                    this->cpu_ctx, this->one, 1);
    this->one = (float*)malloc(sizeof(float));
    this->one[0] = 1.0;
    //base_device::memory::resize_memory_op<float, base_device::DEVICE_CPU>()(
    //                    this->cpu_ctx, this->neg_one, 1);
    this->neg_one = (float*)malloc(sizeof(float));
    this->neg_one[0] = -1.0;
}

// Specialized templates to support std:: complex<double>types
template <>
const_nums<std::complex<double>>::const_nums()
{
    //base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_CPU>()(
    //                    this->cpu_ctx, this->zero, 1);
    this->zero = (std::complex<double>*)malloc(sizeof(std::complex<double>));
    this->zero[0] = std::complex<double>(0.0, 0.0);
    //base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_CPU>()(
    //                    this->cpu_ctx, this->one, 1);
    this->one = (std::complex<double>*)malloc(sizeof(std::complex<double>));
    this->one[0] = std::complex<double>(1.0, 0.0);
    //base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_CPU>()(
    //                    this->cpu_ctx, this->neg_one, 1);
    this->neg_one = (std::complex<double>*)malloc(sizeof(std::complex<double>));
    this->neg_one[0] = std::complex<double>(-1.0, 0.0);
}

// Specialized templates to support std:: complex<float>types
template <>
const_nums<std::complex<float>>::const_nums()
{
    //base_device::memory::resize_memory_op<std::complex<float>, base_device::DEVICE_CPU>()(
    //                    this->cpu_ctx, this->zero, 1);
    this->zero = (std::complex<float>*)malloc(sizeof(std::complex<float>));
    this->zero[0] = std::complex<float>(0.0, 0.0);
    //base_device::memory::resize_memory_op<std::complex<float>, base_device::DEVICE_CPU>()(
    //                    this->cpu_ctx, this->one, 1);
    this->one = (std::complex<float>*)malloc(sizeof(std::complex<float>));
    this->one[0] = std::complex<float>(1.0, 0.0);
    //base_device::memory::resize_memory_op<std::complex<float>, base_device::DEVICE_CPU>()(
    //                    this->cpu_ctx, this->neg_one, 1);
    this->neg_one = (std::complex<float>*)malloc(sizeof(std::complex<float>));
    this->neg_one[0] = std::complex<float>(-1.0, 0.0);
}