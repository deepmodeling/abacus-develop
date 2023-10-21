#include <ATen/kernels/memory_op.h>
#include <base/macros/cuda.h>

#include <cuda_runtime.h>
#include <thrust/complex.h>

namespace container {
namespace op {

template <typename T>
__global__ void set_memory(
    T* out,
    const T var,
    const size_t size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size) {return;}
    out[idx] = var;
}

template <typename T_out, typename T_in>
__global__ void cast_memory(
    T_out* out,
    const T_in* in,
    const int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size) {return;}
    out[idx] = static_cast<T_out>(in[idx]);
}

template <typename T_out, typename T_in>
__global__ void cast_memory(
    std::complex<T_out>* out,
    const std::complex<T_in>* in,
    const int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size) {return;}
    auto* _out = reinterpret_cast<thrust::complex<T_out>*>(out);
    const auto* _in = reinterpret_cast<const thrust::complex<T_in>*>(in);
    _out[idx] = static_cast<thrust::complex<T_out>>(_in[idx]);
}

template <typename T>
void resize_memory_op<T, container::DEVICE_GPU>::operator()(
    T*& arr,
    const size_t size,
    const char* record_in)
{
    if (arr != nullptr) {
        delete_memory_op<T, container::DEVICE_GPU>()(arr);
    }
    cudaMalloc((void **)&arr, sizeof(T) * size);
}

template <typename T>
void set_memory_op<T, container::DEVICE_GPU>::operator()(
    T* arr,
    const T& var, 
    const size_t& size) 
{
    const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_memory<<<block, THREADS_PER_BLOCK>>>(arr, var, size);
}

template <typename T>
void synchronize_memory_op<T, container::DEVICE_CPU, container::DEVICE_GPU>::operator()(
    T* arr_out,
    const T* arr_in,
    const size_t size) 
{
    cudaMemcpy(arr_out, arr_in, sizeof(T) * size, cudaMemcpyDeviceToHost);
}

template <typename T>
void synchronize_memory_op<T, container::DEVICE_GPU, container::DEVICE_CPU>::operator()(
    T* arr_out,
    const T* arr_in,
    const size_t size) 
{
    cudaMemcpy(arr_out, arr_in, sizeof(T) * size, cudaMemcpyHostToDevice);
}

template <typename T>
void synchronize_memory_op<T, container::DEVICE_GPU, container::DEVICE_GPU>::operator()(
    T* arr_out,
    const T* arr_in,
    const size_t size) 
{
    cudaMemcpy(arr_out, arr_in, sizeof(T) * size, cudaMemcpyDeviceToDevice);
}

template <typename T_out, typename T_in>
struct cast_memory_op<T_out, T_in, container::DEVICE_GPU, container::DEVICE_GPU> {
    void operator()(
        T_out* arr_out,
        const T_in* arr_in,
        const size_t size) 
    {
        const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        cast_memory<<<block, THREADS_PER_BLOCK>>>(arr_out, arr_in, size);
    }
};

template <typename T_out, typename T_in>
struct cast_memory_op<T_out, T_in, container::DEVICE_GPU, container::DEVICE_CPU> {
    void operator()(
        T_out* arr_out,
        const T_in* arr_in,
        const size_t size) 
    {
        T_in * arr = nullptr;
        cudaMalloc((void **)&arr, sizeof(T_in) * size);
        cudaMemcpy(arr, arr_in, sizeof(T_in) * size, cudaMemcpyHostToDevice);
        const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        cast_memory<<<block, THREADS_PER_BLOCK>>>(arr_out, arr, size);
        cudaFree(arr);
    }
};

template <typename T_out, typename T_in>
struct cast_memory_op<T_out, T_in, container::DEVICE_CPU, container::DEVICE_GPU> {
    void operator()(
        T_out* arr_out,
        const T_in* arr_in,
        const size_t size) 
    {
        auto * arr = (T_in*) malloc(sizeof(T_in) * size);
        cudaMemcpy(arr, arr_in, sizeof(T_in) * size, cudaMemcpyDeviceToHost);
        for (int ii = 0; ii < size; ii++) {
            arr_out[ii] = static_cast<T_out>(arr[ii]);
        }
        free(arr);
    }
};

template <typename T>
void delete_memory_op<T, container::DEVICE_GPU>::operator() (
    T* arr)
{
    cudaFree(arr);
}

template struct resize_memory_op<int, container::DEVICE_GPU>;
template struct resize_memory_op<int64_t, container::DEVICE_GPU>;
template struct resize_memory_op<float, container::DEVICE_GPU>;
template struct resize_memory_op<double, container::DEVICE_GPU>;
template struct resize_memory_op<std::complex<float>, container::DEVICE_GPU>;
template struct resize_memory_op<std::complex<double>, container::DEVICE_GPU>;

template struct set_memory_op<int, container::DEVICE_GPU>;
template struct set_memory_op<int64_t , container::DEVICE_GPU>;
template struct set_memory_op<float, container::DEVICE_GPU>;
template struct set_memory_op<double, container::DEVICE_GPU>;
template struct set_memory_op<std::complex<float>, container::DEVICE_GPU>;
template struct set_memory_op<std::complex<double>, container::DEVICE_GPU>;

template struct synchronize_memory_op<int, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct synchronize_memory_op<int, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct synchronize_memory_op<int, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct synchronize_memory_op<int64_t, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct synchronize_memory_op<int64_t, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct synchronize_memory_op<int64_t, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct synchronize_memory_op<float, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct synchronize_memory_op<float, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct synchronize_memory_op<float, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct synchronize_memory_op<double, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct synchronize_memory_op<double, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct synchronize_memory_op<double, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<float>, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<float>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<float>, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<double>, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<double>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<double>, container::DEVICE_GPU, container::DEVICE_GPU>;

template struct cast_memory_op<float, float, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory_op<double, double, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory_op<float, double, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory_op<double, float, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory_op<float, float, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory_op<double, double, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory_op<float, double, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory_op<double, float, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory_op<float, float, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory_op<double, double, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory_op<float, double, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory_op<double, float, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, container::DEVICE_CPU, container::DEVICE_GPU>;

template struct delete_memory_op<int, container::DEVICE_GPU>;
template struct delete_memory_op<int64_t, container::DEVICE_GPU>;
template struct delete_memory_op<float, container::DEVICE_GPU>;
template struct delete_memory_op<double, container::DEVICE_GPU>;
template struct delete_memory_op<std::complex<float>, container::DEVICE_GPU>;
template struct delete_memory_op<std::complex<double>, container::DEVICE_GPU>;

} // end of namespace container
} // end of namespace op