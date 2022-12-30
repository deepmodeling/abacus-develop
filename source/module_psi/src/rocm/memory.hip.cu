#include "module_psi/include/memory.h"

#include <vector>
#include <stdio.h>
#include <complex>
#include <assert.h>

#include <hip/hip_runtime.h>

namespace psi {
namespace memory {

template <typename FPTYPE>
void resize_memory_op<FPTYPE, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* dev, 
    FPTYPE*& arr, 
    const size_t size)
{
  if (arr != nullptr) {
    delete_memory_op<FPTYPE, psi::DEVICE_GPU>()(dev, arr);
  }
  hipMalloc((void **)&arr, sizeof(FPTYPE) * size);
}

template <typename FPTYPE>
void set_memory_op<FPTYPE, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* dev, 
    FPTYPE* arr, 
    const int var, 
    const size_t size) 
{
  hipMemset(arr, var, sizeof(FPTYPE) * size);  
}

template <typename FPTYPE> 
void synchronize_memory_op<FPTYPE, psi::DEVICE_CPU, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_CPU* dev_out, 
    const psi::DEVICE_GPU* dev_in, 
    FPTYPE* arr_out, 
    const FPTYPE* arr_in,
    const size_t size) 
{
  hipMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, hipMemcpyDeviceToHost);  
}

template <typename FPTYPE> 
void synchronize_memory_op<FPTYPE, psi::DEVICE_GPU, psi::DEVICE_CPU>::operator()(
    const psi::DEVICE_GPU* dev_out, 
    const psi::DEVICE_CPU* dev_in, 
    FPTYPE* arr_out, 
    const FPTYPE* arr_in,
    const size_t size) 
{
  hipMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, hipMemcpyHostToDevice);  
}

template <typename FPTYPE> 
void synchronize_memory_op<FPTYPE, psi::DEVICE_GPU, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* dev_out, 
    const psi::DEVICE_GPU* dev_in, 
    FPTYPE* arr_out, 
    const FPTYPE* arr_in,
    const size_t size) 
{
  hipMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, hipMemcpyDeviceToDevice);  
}

template <typename FPTYPE>
void delete_memory_op<FPTYPE, psi::DEVICE_GPU>::operator() (
    const psi::DEVICE_GPU* dev, 
    FPTYPE* arr) 
{
  hipFree(arr);
}

template struct resize_memory_op<int, psi::DEVICE_GPU>;
template struct resize_memory_op<double, psi::DEVICE_GPU>;
template struct resize_memory_op<std::complex<double>, psi::DEVICE_GPU>;

template struct set_memory_op<int, psi::DEVICE_GPU>;
template struct set_memory_op<double, psi::DEVICE_GPU>;
template struct set_memory_op<std::complex<double>, psi::DEVICE_GPU>;

template struct synchronize_memory_op<int, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<int, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<int, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<double, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<double, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<double, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_GPU>;

template struct delete_memory_op<int, psi::DEVICE_GPU>;
template struct delete_memory_op<double, psi::DEVICE_GPU>;
template struct delete_memory_op<std::complex<double>, psi::DEVICE_GPU>;

} // end of namespace gpu_cuda
} // end of namespace psi