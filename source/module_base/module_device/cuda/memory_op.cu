#include "module_base/module_device/memory_op.h"
#include "module_base/memory.h"

#include <base/macros/macros.h>
#include <cuda_runtime.h>
#include <thrust/complex.h>

#include <complex>
#include <type_traits>

#define THREADS_PER_BLOCK 256

namespace base_device
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

template <typename FPTYPE_out, typename FPTYPE_in>
__global__ void cast_memory(std::complex<FPTYPE_out>* out, const FPTYPE_in* in, const int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= size)
    {
        return;
    }
    auto* _out = reinterpret_cast<thrust::complex<FPTYPE_out>*>(out);
    _out[idx] = static_cast<thrust::complex<FPTYPE_out>>(in[idx]);
}

template <typename FPTYPE>
void resize_memory_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* dev,
                                                                   FPTYPE*& arr,
                                                                   const size_t size,
                                                                   const char* record_in)
{
    if (arr != nullptr)
    {
        delete_memory_op<FPTYPE, base_device::DEVICE_GPU>()(dev, arr);
    }
    cudaErrcheck(cudaMalloc((void**)&arr, sizeof(FPTYPE) * size));
    std::string record_string;
    if (record_in != nullptr)
    {
        record_string = record_in;
    }
    else
    {
        record_string = "no_record";
    }

    if (record_string != "no_record")
    {
        ModuleBase::Memory::record_gpu(record_string, sizeof(FPTYPE) * size);
    }
}

template <typename FPTYPE>
void set_memory_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* dev,
                                                                FPTYPE* arr,
                                                                const int var,
                                                                const size_t size)
{
    cudaErrcheck(cudaMemset(arr, var, sizeof(FPTYPE) * size));
}

template <typename FPTYPE>
void synchronize_memory_op<FPTYPE, base_device::DEVICE_CPU, base_device::DEVICE_GPU>::operator()(
    const base_device::DEVICE_CPU* dev_out,
    const base_device::DEVICE_GPU* dev_in,
    FPTYPE* arr_out,
    const FPTYPE* arr_in,
    const size_t size)
{
    cudaErrcheck(cudaMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, cudaMemcpyDeviceToHost));
}

template <typename FPTYPE>
void synchronize_memory_op<FPTYPE, base_device::DEVICE_GPU, base_device::DEVICE_CPU>::operator()(
    const base_device::DEVICE_GPU* dev_out,
    const base_device::DEVICE_CPU* dev_in,
    FPTYPE* arr_out,
    const FPTYPE* arr_in,
    const size_t size)
{
    cudaErrcheck(cudaMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, cudaMemcpyHostToDevice));
}

template <typename FPTYPE>
void synchronize_memory_op<FPTYPE, base_device::DEVICE_GPU, base_device::DEVICE_GPU>::operator()(
    const base_device::DEVICE_GPU* dev_out,
    const base_device::DEVICE_GPU* dev_in,
    FPTYPE* arr_out,
    const FPTYPE* arr_in,
    const size_t size)
{
    cudaErrcheck(cudaMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, cudaMemcpyDeviceToDevice));
}

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, base_device::DEVICE_GPU, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* dev_out,
                    const base_device::DEVICE_GPU* dev_in,
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

        cudaCheckOnDebug();
    }
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, base_device::DEVICE_GPU, base_device::DEVICE_CPU> {
    void operator()(const base_device::DEVICE_GPU* dev_out,
                    const base_device::DEVICE_CPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size) {

        if (size == 0) {return;}
        // No need to cast the memory if the data types are the same.
        if (std::is_same<FPTYPE_out, FPTYPE_in>::value)
        {
            synchronize_memory_op<FPTYPE_out, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(dev_out,
                                                                                                  dev_in,
                                                                                                  arr_out,
                                                                                                  reinterpret_cast<const FPTYPE_out*>(arr_in),
                                                                                                  size);
            return;
        }
        FPTYPE_in * arr = nullptr;
        cudaErrcheck(cudaMalloc((void **)&arr, sizeof(FPTYPE_in) * size));
        cudaErrcheck(cudaMemcpy(arr, arr_in, sizeof(FPTYPE_in) * size, cudaMemcpyHostToDevice));
        const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        cast_memory<<<block, THREADS_PER_BLOCK>>>(arr_out, arr, size);
        cudaCheckOnDebug();
        cudaErrcheck(cudaFree(arr));
    }
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, base_device::DEVICE_CPU, base_device::DEVICE_GPU> {
    void operator()(const base_device::DEVICE_CPU* dev_out,
                    const base_device::DEVICE_GPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size) {
        if (size == 0) {return;}
        // No need to cast the memory if the data types are the same.
        if (std::is_same<FPTYPE_out, FPTYPE_in>::value)
        {
            synchronize_memory_op<FPTYPE_out, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(dev_out,
                                                                                                  dev_in,
                                                                                                  arr_out,
                                                                                                  reinterpret_cast<const FPTYPE_out*>(arr_in),
                                                                                                  size);
            return;
        }
        auto * arr = (FPTYPE_in*) malloc(sizeof(FPTYPE_in) * size);
        cudaErrcheck(cudaMemcpy(arr, arr_in, sizeof(FPTYPE_in) * size, cudaMemcpyDeviceToHost));
        for (int ii = 0; ii < size; ii++) {
            arr_out[ii] = static_cast<FPTYPE_out>(arr[ii]);
        }
        free(arr);
    }
};

template <typename FPTYPE>
void delete_memory_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* dev, FPTYPE* arr)
{
    cudaErrcheck(cudaFree(arr));
}

template struct resize_memory_op<int, base_device::DEVICE_GPU>;
template struct resize_memory_op<float, base_device::DEVICE_GPU>;
template struct resize_memory_op<double, base_device::DEVICE_GPU>;
template struct resize_memory_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>;

template struct set_memory_op<int, base_device::DEVICE_GPU>;
template struct set_memory_op<float, base_device::DEVICE_GPU>;
template struct set_memory_op<double, base_device::DEVICE_GPU>;
template struct set_memory_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct set_memory_op<std::complex<double>, base_device::DEVICE_GPU>;

template struct synchronize_memory_op<int, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;
template struct synchronize_memory_op<int, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
template struct synchronize_memory_op<int, base_device::DEVICE_GPU, base_device::DEVICE_GPU>;
template struct synchronize_memory_op<float, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;
template struct synchronize_memory_op<float, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
template struct synchronize_memory_op<float, base_device::DEVICE_GPU, base_device::DEVICE_GPU>;
template struct synchronize_memory_op<double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;
template struct synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
template struct synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<float>, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<float>, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<float>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<double>, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>;

template struct cast_memory_op<float, float, base_device::DEVICE_GPU, base_device::DEVICE_GPU>;
template struct cast_memory_op<double, double, base_device::DEVICE_GPU, base_device::DEVICE_GPU>;
template struct cast_memory_op<float, double, base_device::DEVICE_GPU, base_device::DEVICE_GPU>;
template struct cast_memory_op<double, float, base_device::DEVICE_GPU, base_device::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>,
                               std::complex<float>,
                               base_device::DEVICE_GPU,
                               base_device::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>,
                               std::complex<double>,
                               base_device::DEVICE_GPU,
                               base_device::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>,
                               std::complex<double>,
                               base_device::DEVICE_GPU,
                               base_device::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>,
                               std::complex<float>,
                               base_device::DEVICE_GPU,
                               base_device::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, float, base_device::DEVICE_GPU, base_device::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, double, base_device::DEVICE_GPU, base_device::DEVICE_GPU>;
template struct cast_memory_op<float, float, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
template struct cast_memory_op<double, double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
template struct cast_memory_op<float, double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
template struct cast_memory_op<double, float, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>,
                               std::complex<float>,
                               base_device::DEVICE_GPU,
                               base_device::DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>,
                               std::complex<double>,
                               base_device::DEVICE_GPU,
                               base_device::DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>,
                               std::complex<double>,
                               base_device::DEVICE_GPU,
                               base_device::DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>,
                               std::complex<float>,
                               base_device::DEVICE_GPU,
                               base_device::DEVICE_CPU>;
template struct cast_memory_op<float, float, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;
template struct cast_memory_op<double, double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;
template struct cast_memory_op<float, double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;
template struct cast_memory_op<double, float, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>,
                               std::complex<float>,
                               base_device::DEVICE_CPU,
                               base_device::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>,
                               std::complex<double>,
                               base_device::DEVICE_CPU,
                               base_device::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>,
                               std::complex<double>,
                               base_device::DEVICE_CPU,
                               base_device::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>,
                               std::complex<float>,
                               base_device::DEVICE_CPU,
                               base_device::DEVICE_GPU>;

template struct delete_memory_op<int, base_device::DEVICE_GPU>;
template struct delete_memory_op<float, base_device::DEVICE_GPU>;
template struct delete_memory_op<double, base_device::DEVICE_GPU>;
template struct delete_memory_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct delete_memory_op<float*, base_device::DEVICE_GPU>;
template struct delete_memory_op<double*, base_device::DEVICE_GPU>;
template struct delete_memory_op<std::complex<float>*, base_device::DEVICE_GPU>;
template struct delete_memory_op<std::complex<double>*, base_device::DEVICE_GPU>;

} // namespace memory
} // end of namespace base_device