#include <ATen/kernels/memory.h>
#include <base/macros/cuda.h>

#include <cuda_runtime.h>
#include <thrust/complex.h>

#include "module_base/memory.h"

namespace container {
namespace kernels {

template <typename T>
__global__ void do_set_memory(
    T* out,
    const T var,
    const size_t size)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size) {return;}
    out[idx] = var;
}

template <typename T_out, typename T_in>
__global__ void do_cast_memory(
    T_out* out,
    const T_in* in,
    const int size)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size) {return;}
    out[idx] = static_cast<T_out>(in[idx]);
}

template <typename T_out, typename T_in>
__global__ void do_cast_memory(
    std::complex<T_out>* out,
    const std::complex<T_in>* in,
    const int size)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size) {return;}
    auto* _out = reinterpret_cast<thrust::complex<T_out>*>(out);
    const auto* _in = reinterpret_cast<const thrust::complex<T_in>*>(in);
    _out[idx] = static_cast<thrust::complex<T_out>>(_in[idx]);
}

template <typename T, typename Device>
void resize_memory<T, Device>::operator()(
    T*& arr,
    const size_t& size,
    const char* record_in)
{
    if (arr != nullptr) {
        delete_memory<T, container::DEVICE_GPU>()(arr);
        //ModuleBase::Memory::record_gpu("memory","resize_free",sizeof(arr));
    }
    cudaErrcheck(cudaMalloc((void **)&arr, sizeof(T) * size));
    ModuleBase::Memory::record_gpu("memory","Mem::resize",sizeof(T) * size);
}

template <typename T, typename Device>
void set_memory<T, Device>::operator()(
    T* arr,
    const T& var, 
    const size_t& size) 
{
    const int block = static_cast<int>((size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
    do_set_memory<<<block, THREADS_PER_BLOCK>>>(arr, var, size);
    cudaErrcheck(cudaGetLastError());
    cudaErrcheck(cudaDeviceSynchronize());
}

template <typename T>
struct synchronize_memory<T, DEVICE_CPU, DEVICE_GPU> {
    void operator()(
        T *arr_out,
        const T *arr_in,
        const size_t& size)
    {
        cudaErrcheck(cudaMemcpy(arr_out, arr_in, sizeof(T) * size, cudaMemcpyDeviceToHost));
    }
};

template <typename T>
struct synchronize_memory<T, DEVICE_GPU, DEVICE_CPU> {
    void operator()(
        T *arr_out,
        const T *arr_in,
        const size_t& size)
    {
        cudaErrcheck(cudaMemcpy(arr_out, arr_in, sizeof(T) * size, cudaMemcpyHostToDevice));
        //ModuleBase::Memory::record_gpu("memory","synchronize",sizeof(T) * size);
    }
};

template <typename T>
struct synchronize_memory<T, DEVICE_GPU, DEVICE_GPU> {
    void operator()(
        T *arr_out,
        const T *arr_in,
        const size_t& size)
    {
        cudaErrcheck(cudaMemcpy(arr_out, arr_in, sizeof(T) * size, cudaMemcpyHostToDevice));
        //ModuleBase::Memory::record_gpu("memory","synchronize",sizeof(T) * size);
    }
};


template <typename T_out, typename T_in>
struct cast_memory<T_out, T_in, container::DEVICE_GPU, container::DEVICE_GPU> {
    void operator()(
        T_out* arr_out,
        const T_in* arr_in,
        const size_t& size)
    {
        const int block = static_cast<int>((size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
        do_cast_memory<<<block, THREADS_PER_BLOCK>>>(arr_out, arr_in, size);
        cudaErrcheck(cudaGetLastError());
        cudaErrcheck(cudaDeviceSynchronize());
    }
};


template <typename T_out, typename T_in>
struct cast_memory<T_out, T_in, container::DEVICE_GPU, container::DEVICE_CPU> {
    void operator()(
        T_out* arr_out,
        const T_in* arr_in,
        const size_t& size)
    {
        T_in * arr = nullptr;
        cudaErrcheck(cudaMalloc((void **)&arr, sizeof(T_in) * size));
        ModuleBase::Memory::record_gpu("memory","Mem::cast",sizeof(T_in) * size);
        cudaErrcheck(cudaMemcpy(arr, arr_in, sizeof(T_in) * size, cudaMemcpyHostToDevice));
        //ModuleBase::Memory::record_gpu("memory","cast",sizeof(T_in) * size);
        const int block = static_cast<int>((size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
        do_cast_memory<<<block, THREADS_PER_BLOCK>>>(arr_out, arr, size);
        cudaErrcheck(cudaGetLastError());
        cudaErrcheck(cudaDeviceSynchronize());
        cudaErrcheck(cudaFree(arr));
        ModuleBase::Memory::record_gpu("memory","cast_free",sizeof(arr));
    }
};


template <typename T_out, typename T_in>
struct cast_memory<T_out, T_in, container::DEVICE_CPU, container::DEVICE_GPU> {
    void operator()(
        T_out* arr_out,
        const T_in* arr_in,
        const size_t& size)
    {
        auto * arr = (T_in*) malloc(sizeof(T_in) * size);
        cudaErrcheck(cudaMemcpy(arr, arr_in, sizeof(T_in) * size, cudaMemcpyDeviceToHost));
        for (int ii = 0; ii < size; ii++) {
            arr_out[ii] = static_cast<T_out>(arr[ii]);
        }
        free(arr);
    }
};

template <typename T, typename Device>
void delete_memory<T, Device>::operator() (
    T* arr)
{
    cudaErrcheck(cudaFree(arr));
}

template struct resize_memory<int, container::DEVICE_GPU>;
template struct resize_memory<int64_t, container::DEVICE_GPU>;
template struct resize_memory<float, container::DEVICE_GPU>;
template struct resize_memory<double, container::DEVICE_GPU>;
template struct resize_memory<std::complex<float>, container::DEVICE_GPU>;
template struct resize_memory<std::complex<double>, container::DEVICE_GPU>;

template struct set_memory<int, container::DEVICE_GPU>;
template struct set_memory<int64_t , container::DEVICE_GPU>;
template struct set_memory<float, container::DEVICE_GPU>;
template struct set_memory<double, container::DEVICE_GPU>;
template struct set_memory<std::complex<float>, container::DEVICE_GPU>;
template struct set_memory<std::complex<double>, container::DEVICE_GPU>;

template struct synchronize_memory<int, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct synchronize_memory<int, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct synchronize_memory<int, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct synchronize_memory<int64_t, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct synchronize_memory<int64_t, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct synchronize_memory<int64_t, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct synchronize_memory<float, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct synchronize_memory<float, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct synchronize_memory<float, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct synchronize_memory<double, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct synchronize_memory<double, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct synchronize_memory<double, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct synchronize_memory<std::complex<float>, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct synchronize_memory<std::complex<float>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct synchronize_memory<std::complex<float>, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct synchronize_memory<std::complex<double>, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct synchronize_memory<std::complex<double>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct synchronize_memory<std::complex<double>, container::DEVICE_GPU, container::DEVICE_GPU>;

template struct cast_memory<float, float, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<double, double, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<float, double, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<double, float, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<float>, std::complex<float>, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<double>, std::complex<double>, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<float>, std::complex<double>, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<double>, std::complex<float>, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<float, float, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<double, double, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<float, double, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<double, float, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<std::complex<float>, std::complex<float>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<std::complex<double>, std::complex<double>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<std::complex<float>, std::complex<double>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<std::complex<double>, std::complex<float>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<float, float, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory<double, double, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory<float, double, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory<double, float, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<float>, std::complex<float>, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<double>, std::complex<double>, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<float>, std::complex<double>, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<double>, std::complex<float>, container::DEVICE_CPU, container::DEVICE_GPU>;

template struct delete_memory<int, container::DEVICE_GPU>;
template struct delete_memory<int64_t, container::DEVICE_GPU>;
template struct delete_memory<float, container::DEVICE_GPU>;
template struct delete_memory<double, container::DEVICE_GPU>;
template struct delete_memory<std::complex<float>, container::DEVICE_GPU>;
template struct delete_memory<std::complex<double>, container::DEVICE_GPU>;

} // namespace kernels
} // namespace container