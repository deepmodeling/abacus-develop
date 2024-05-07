#include "../memory_op.h"

#include <base/macros/macros.h>
#include <cuda_runtime.h>
#include <thrust/complex.h>

#define THREADS_PER_BLOCK 256

namespace device
{
namespace memory
{

template <typename FPTYPE_out, typename FPTYPE_in>
__global__ void cast_memory(FPTYPE_out* out, const FPTYPE_in* in, const int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= size)
    {
        return;
    }
    out[idx] = static_cast<FPTYPE_out>(in[idx]);
}

template <typename FPTYPE_out, typename FPTYPE_in>
__global__ void cast_memory(std::complex<FPTYPE_out>* out, const std::complex<FPTYPE_in>* in, const int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= size)
    {
        return;
    }
    auto* _out = reinterpret_cast<thrust::complex<FPTYPE_out>*>(out);
    const auto* _in = reinterpret_cast<const thrust::complex<FPTYPE_in>*>(in);
    _out[idx] = static_cast<thrust::complex<FPTYPE_out>>(_in[idx]);
}

template <typename FPTYPE>
void resize_memory_op<FPTYPE, device::GPU>::operator()(const device::GPU* dev,
                                                              FPTYPE*& arr,
                                                              const size_t size,
                                                              const char* record_in)
{
    if (arr != nullptr)
    {
        delete_memory_op<FPTYPE, device::GPU>()(dev, arr);
    }
    cudaErrcheck(cudaMalloc((void**)&arr, sizeof(FPTYPE) * size));
}

template <typename FPTYPE>
void set_memory_op<FPTYPE, device::GPU>::operator()(const device::GPU* dev,
                                                           FPTYPE* arr,
                                                           const int var,
                                                           const size_t size)
{
    cudaErrcheck(cudaMemset(arr, var, sizeof(FPTYPE) * size));
}

template <typename FPTYPE>
void synchronize_memory_op<FPTYPE, device::CPU, device::GPU>::operator()(
    const device::CPU* dev_out,
    const device::GPU* dev_in,
    FPTYPE* arr_out,
    const FPTYPE* arr_in,
    const size_t size)
{
    cudaErrcheck(cudaMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, cudaMemcpyDeviceToHost));
}

template <typename FPTYPE>
void synchronize_memory_op<FPTYPE, device::GPU, device::CPU>::operator()(
    const device::GPU* dev_out,
    const device::CPU* dev_in,
    FPTYPE* arr_out,
    const FPTYPE* arr_in,
    const size_t size)
{
    cudaErrcheck(cudaMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, cudaMemcpyHostToDevice));
}

template <typename FPTYPE>
void synchronize_memory_op<FPTYPE, device::GPU, device::GPU>::operator()(
    const device::GPU* dev_out,
    const device::GPU* dev_in,
    FPTYPE* arr_out,
    const FPTYPE* arr_in,
    const size_t size)
{
    cudaErrcheck(cudaMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, cudaMemcpyDeviceToDevice));
}

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, device::GPU, device::GPU>
{
    void operator()(const device::GPU* dev_out,
                    const device::GPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size)
    {
        if (size == 0)
        {
            return;
        }
        const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        cast_memory<<<block, THREADS_PER_BLOCK>>>(arr_out, arr_in, size);

        cudaErrcheck(cudaGetLastError());
        cudaErrcheck(cudaDeviceSynchronize());
    }
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, device::GPU, device::CPU>
{
    void operator()(const device::GPU* dev_out,
                    const device::CPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size)
    {

        if (size == 0)
        {
            return;
        }
        FPTYPE_in* arr = nullptr;
        cudaErrcheck(cudaMalloc((void**)&arr, sizeof(FPTYPE_in) * size));
        cudaErrcheck(cudaMemcpy(arr, arr_in, sizeof(FPTYPE_in) * size, cudaMemcpyHostToDevice));
        const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        cast_memory<<<block, THREADS_PER_BLOCK>>>(arr_out, arr, size);
        cudaErrcheck(cudaGetLastError());
        cudaErrcheck(cudaDeviceSynchronize());
        cudaErrcheck(cudaFree(arr));
    }
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, device::CPU, device::GPU>
{
    void operator()(const device::CPU* dev_out,
                    const device::GPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size)
    {
        auto* arr = (FPTYPE_in*)malloc(sizeof(FPTYPE_in) * size);
        cudaErrcheck(cudaMemcpy(arr, arr_in, sizeof(FPTYPE_in) * size, cudaMemcpyDeviceToHost));
        for (int ii = 0; ii < size; ii++)
        {
            arr_out[ii] = static_cast<FPTYPE_out>(arr[ii]);
        }
        free(arr);
    }
};

template <typename FPTYPE>
void delete_memory_op<FPTYPE, device::GPU>::operator()(const device::GPU* dev, FPTYPE* arr)
{
    cudaErrcheck(cudaFree(arr));
}

template struct resize_memory_op<int, device::GPU>;
template struct resize_memory_op<float, device::GPU>;
template struct resize_memory_op<double, device::GPU>;
template struct resize_memory_op<std::complex<float>, device::GPU>;
template struct resize_memory_op<std::complex<double>, device::GPU>;

template struct set_memory_op<int, device::GPU>;
template struct set_memory_op<float, device::GPU>;
template struct set_memory_op<double, device::GPU>;
template struct set_memory_op<std::complex<float>, device::GPU>;
template struct set_memory_op<std::complex<double>, device::GPU>;

template struct synchronize_memory_op<int, device::CPU, device::GPU>;
template struct synchronize_memory_op<int, device::GPU, device::CPU>;
template struct synchronize_memory_op<int, device::GPU, device::GPU>;
template struct synchronize_memory_op<float, device::CPU, device::GPU>;
template struct synchronize_memory_op<float, device::GPU, device::CPU>;
template struct synchronize_memory_op<float, device::GPU, device::GPU>;
template struct synchronize_memory_op<double, device::CPU, device::GPU>;
template struct synchronize_memory_op<double, device::GPU, device::CPU>;
template struct synchronize_memory_op<double, device::GPU, device::GPU>;
template struct synchronize_memory_op<std::complex<float>, device::CPU, device::GPU>;
template struct synchronize_memory_op<std::complex<float>, device::GPU, device::CPU>;
template struct synchronize_memory_op<std::complex<float>, device::GPU, device::GPU>;
template struct synchronize_memory_op<std::complex<double>, device::CPU, device::GPU>;
template struct synchronize_memory_op<std::complex<double>, device::GPU, device::CPU>;
template struct synchronize_memory_op<std::complex<double>, device::GPU, device::GPU>;

template struct cast_memory_op<float, float, device::GPU, device::GPU>;
template struct cast_memory_op<double, double, device::GPU, device::GPU>;
template struct cast_memory_op<float, double, device::GPU, device::GPU>;
template struct cast_memory_op<double, float, device::GPU, device::GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, device::GPU, device::GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, device::GPU, device::GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, device::GPU, device::GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, device::GPU, device::GPU>;
template struct cast_memory_op<float, float, device::GPU, device::CPU>;
template struct cast_memory_op<double, double, device::GPU, device::CPU>;
template struct cast_memory_op<float, double, device::GPU, device::CPU>;
template struct cast_memory_op<double, float, device::GPU, device::CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, device::GPU, device::CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, device::GPU, device::CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, device::GPU, device::CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, device::GPU, device::CPU>;
template struct cast_memory_op<float, float, device::CPU, device::GPU>;
template struct cast_memory_op<double, double, device::CPU, device::GPU>;
template struct cast_memory_op<float, double, device::CPU, device::GPU>;
template struct cast_memory_op<double, float, device::CPU, device::GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, device::CPU, device::GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, device::CPU, device::GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, device::CPU, device::GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, device::CPU, device::GPU>;

template struct delete_memory_op<int, device::GPU>;
template struct delete_memory_op<float, device::GPU>;
template struct delete_memory_op<double, device::GPU>;
template struct delete_memory_op<std::complex<float>, device::GPU>;
template struct delete_memory_op<std::complex<double>, device::GPU>;

} // namespace memory
} // end of namespace device