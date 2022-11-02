#ifndef MODULE_PSI_GPU_CUDA_H_
#define MODULE_PSI_GPU_CUDA_H_

#include <vector>
#include <stdio.h>
#include <complex>
#include <assert.h>
#include <cuda_runtime.h>
#include "module_psi/include/memory.h"

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
  cudaMalloc((void **)&arr, sizeof(FPTYPE) * size);
}

template <typename FPTYPE>
void set_memory_op<FPTYPE, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* dev, 
    FPTYPE* arr, 
    const int var, 
    const size_t size) 
{
  cudaMemset(arr, var, sizeof(FPTYPE) * size);  
}

template <typename FPTYPE> 
void synchronize_memory_op<FPTYPE, psi::DEVICE_CPU, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_CPU* dev_out, 
    const psi::DEVICE_GPU* dev_in, 
    FPTYPE* arr_out, 
    const FPTYPE* arr_in,
    const size_t size) 
{
  cudaMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, cudaMemcpyDeviceToHost);  
}

template <typename FPTYPE> 
void synchronize_memory_op<FPTYPE, psi::DEVICE_GPU, psi::DEVICE_CPU>::operator()(
    const psi::DEVICE_GPU* dev_out, 
    const psi::DEVICE_CPU* dev_in, 
    FPTYPE* arr_out, 
    const FPTYPE* arr_in,
    const size_t size) 
{
  cudaMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, cudaMemcpyHostToDevice);  
}

template <typename FPTYPE> 
void synchronize_memory_op<FPTYPE, psi::DEVICE_GPU, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* dev_out, 
    const psi::DEVICE_GPU* dev_in, 
    FPTYPE* arr_out, 
    const FPTYPE* arr_in,
    const size_t size) 
{
  cudaMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, cudaMemcpyDeviceToDevice);  
}

template <typename FPTYPE>
void delete_memory_op<FPTYPE, psi::DEVICE_GPU>::operator() (
    const psi::DEVICE_GPU* dev, 
    FPTYPE* arr) 
{
  cudaFree(arr);
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

#endif  // MODULE_PSI_GPU_CUDA_H_