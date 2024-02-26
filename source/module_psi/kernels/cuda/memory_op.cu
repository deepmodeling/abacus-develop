#include "module_psi/kernels/memory_op.h"

#include <complex>
#include <iostream>

#include <cuda_runtime.h>
#include <thrust/complex.h>
#include <base/macros/macros.h>

#include "module_base/memory.h"

#define THREADS_PER_BLOCK 256

namespace psi {
namespace memory {

template <typename FPTYPE_out, typename FPTYPE_in>
__global__ void cast_memory(
        FPTYPE_out* out,
        const FPTYPE_in* in,
        const int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size) {return;}
    out[idx] = static_cast<FPTYPE_out>(in[idx]);
}

template <typename FPTYPE_out, typename FPTYPE_in>
__global__ void cast_memory(
        std::complex<FPTYPE_out>* out,
        const std::complex<FPTYPE_in>* in,
        const int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size) {return;}
    auto* _out = reinterpret_cast<thrust::complex<FPTYPE_out>*>(out);
    const auto* _in = reinterpret_cast<const thrust::complex<FPTYPE_in>*>(in);
    _out[idx] = static_cast<thrust::complex<FPTYPE_out>>(_in[idx]);
}

template <typename FPTYPE>
void resize_memory_op<FPTYPE, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* dev, 
    FPTYPE*& arr, 
    const size_t size,
    const char* record_in)
{
  if (arr != nullptr) {
    delete_memory_op<FPTYPE, psi::DEVICE_GPU>()(dev, arr);
  }
  cudaErrcheck(cudaMalloc((void **)&arr, sizeof(FPTYPE) * size));
  ModuleBase::Memory::record_gpu("memory_op","MemOp::CMem",sizeof(FPTYPE) * size);
}

template <typename FPTYPE>
void set_memory_op<FPTYPE, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* dev, 
    FPTYPE* arr, 
    const int var, 
    const size_t size) 
{
  cudaErrcheck(cudaMemset(arr, var, sizeof(FPTYPE) * size));  
}

template <typename FPTYPE> 
void synchronize_memory_op<FPTYPE, psi::DEVICE_CPU, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_CPU* dev_out, 
    const psi::DEVICE_GPU* dev_in, 
    FPTYPE* arr_out, 
    const FPTYPE* arr_in,
    const size_t size) 
{
  cudaErrcheck(cudaMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, cudaMemcpyDeviceToHost));  
}

template <typename FPTYPE> 
void synchronize_memory_op<FPTYPE, psi::DEVICE_GPU, psi::DEVICE_CPU>::operator()(
    const psi::DEVICE_GPU* dev_out, 
    const psi::DEVICE_CPU* dev_in, 
    FPTYPE* arr_out, 
    const FPTYPE* arr_in,
    const size_t size) 
{
  cudaErrcheck(cudaMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, cudaMemcpyHostToDevice));  
  //ModuleBase::Memory::record_gpu("memoryOp","synchronize",sizeof(FPTYPE) * size);
}

template <typename FPTYPE> 
void synchronize_memory_op<FPTYPE, psi::DEVICE_GPU, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* dev_out, 
    const psi::DEVICE_GPU* dev_in, 
    FPTYPE* arr_out, 
    const FPTYPE* arr_in,
    const size_t size) 
{
  cudaErrcheck(cudaMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, cudaMemcpyDeviceToDevice));  
  //ModuleBase::Memory::record_gpu("memoryOp","synchronize",sizeof(FPTYPE) * size);
}

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, psi::DEVICE_GPU, psi::DEVICE_GPU> {
    void operator()(const psi::DEVICE_GPU* dev_out,
                    const psi::DEVICE_GPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size) {
        if (size == 0) {return;}
        const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        cast_memory<<<block, THREADS_PER_BLOCK>>>(arr_out, arr_in, size);

        cudaErrcheck(cudaGetLastError());
        cudaErrcheck(cudaDeviceSynchronize());
    }
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, psi::DEVICE_GPU, psi::DEVICE_CPU> {
    void operator()(const psi::DEVICE_GPU* dev_out,
                    const psi::DEVICE_CPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size) {
        
        if (size == 0) {return;}
        FPTYPE_in * arr = nullptr;
        cudaErrcheck(cudaMalloc((void **)&arr, sizeof(FPTYPE_in) * size));
        ModuleBase::Memory::record_gpu("memory_op","MemOp::CMemOp",sizeof(FPTYPE_in) * size);
        cudaErrcheck(cudaMemcpy(arr, arr_in, sizeof(FPTYPE_in) * size, cudaMemcpyHostToDevice));
        //ModuleBase::Memory::record_gpu("memoryOp","synchronize",sizeof(FPTYPE_in) * size);
        const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        cast_memory<<<block, THREADS_PER_BLOCK>>>(arr_out, arr, size);
        cudaErrcheck(cudaGetLastError());
        cudaErrcheck(cudaDeviceSynchronize());
        cudaErrcheck(cudaFree(arr));
    }
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, psi::DEVICE_CPU, psi::DEVICE_GPU> {
    void operator()(const psi::DEVICE_CPU* dev_out,
                    const psi::DEVICE_GPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size) {
        auto * arr = (FPTYPE_in*) malloc(sizeof(FPTYPE_in) * size);
        cudaErrcheck(cudaMemcpy(arr, arr_in, sizeof(FPTYPE_in) * size, cudaMemcpyDeviceToHost));
        for (int ii = 0; ii < size; ii++) {
            arr_out[ii] = static_cast<FPTYPE_out>(arr[ii]);
        }
        free(arr);
    }
};

template <typename FPTYPE>
void delete_memory_op<FPTYPE, psi::DEVICE_GPU>::operator() (
    const psi::DEVICE_GPU* dev, 
    FPTYPE* arr) 
{
  cudaErrcheck(cudaFree(arr));
}

template struct resize_memory_op<int, psi::DEVICE_GPU>;
template struct resize_memory_op<float, psi::DEVICE_GPU>;
template struct resize_memory_op<double, psi::DEVICE_GPU>;
template struct resize_memory_op<std::complex<float>, psi::DEVICE_GPU>;
template struct resize_memory_op<std::complex<double>, psi::DEVICE_GPU>;

template struct set_memory_op<int, psi::DEVICE_GPU>;
template struct set_memory_op<float, psi::DEVICE_GPU>;
template struct set_memory_op<double, psi::DEVICE_GPU>;
template struct set_memory_op<std::complex<float>, psi::DEVICE_GPU>;
template struct set_memory_op<std::complex<double>, psi::DEVICE_GPU>;

template struct synchronize_memory_op<int, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<int, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<int, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<float, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<float, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<float, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<double, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<double, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<double, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<float>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_GPU>;

template struct cast_memory_op<float, float, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<double, double, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<float, double, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<double, float, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<float, float, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<double, double, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<float, double, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<double, float, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<float, float, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct cast_memory_op<double, double, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct cast_memory_op<float, double, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct cast_memory_op<double, float, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, psi::DEVICE_CPU, psi::DEVICE_GPU>;

template struct delete_memory_op<int, psi::DEVICE_GPU>;
template struct delete_memory_op<float, psi::DEVICE_GPU>;
template struct delete_memory_op<double, psi::DEVICE_GPU>;
template struct delete_memory_op<std::complex<float>, psi::DEVICE_GPU>;
template struct delete_memory_op<std::complex<double>, psi::DEVICE_GPU>;

} // end of namespace gpu_cuda
} // end of namespace psi