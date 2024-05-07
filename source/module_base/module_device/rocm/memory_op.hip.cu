#include "../memory_op.h"

#include <base/macros/macros.h>
#include <hip/hip_runtime.h>
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
    hipErrcheck(hipMalloc((void**)&arr, sizeof(FPTYPE) * size));
}

template <typename FPTYPE>
void set_memory_op<FPTYPE, device::GPU>::operator()(const device::GPU* dev,
                                                        FPTYPE* arr,
                                                        const int var,
                                                        const size_t size)
{
    hipErrcheck(hipMemset(arr, var, sizeof(FPTYPE) * size));
}

template <typename FPTYPE>
void synchronize_memory_op<FPTYPE, device::CPU, device::GPU>::operator()(const device::CPU* dev_out,
                                                                                 const device::GPU* dev_in,
                                                                                 FPTYPE* arr_out,
                                                                                 const FPTYPE* arr_in,
                                                                                 const size_t size)
{
    hipErrcheck(hipMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, hipMemcpyDeviceToHost));
}

template <typename FPTYPE>
void synchronize_memory_op<FPTYPE, device::GPU, device::CPU>::operator()(const device::GPU* dev_out,
                                                                                 const device::CPU* dev_in,
                                                                                 FPTYPE* arr_out,
                                                                                 const FPTYPE* arr_in,
                                                                                 const size_t size)
{
    hipErrcheck(hipMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, hipMemcpyHostToDevice));
}

template <typename FPTYPE>
void synchronize_memory_op<FPTYPE, device::GPU, device::GPU>::operator()(const device::GPU* dev_out,
                                                                                 const device::GPU* dev_in,
                                                                                 FPTYPE* arr_out,
                                                                                 const FPTYPE* arr_in,
                                                                                 const size_t size)
{
    hipErrcheck(hipMemcpy(arr_out, arr_in, sizeof(FPTYPE) * size, hipMemcpyDeviceToDevice));
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
        const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        hipLaunchKernelGGL(cast_memory, dim3(block), dim3(THREADS_PER_BLOCK), 0, 0, arr_out, arr_in, size);
        hipErrcheck(hipGetLastError());
        hipErrcheck(hipDeviceSynchronize());
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
        FPTYPE_in* arr = nullptr;
        hipErrcheck(hipMalloc((void**)&arr, sizeof(FPTYPE_in) * size));
        hipErrcheck(hipMemcpy(arr, arr_in, sizeof(FPTYPE_in) * size, hipMemcpyHostToDevice));
        const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        hipLaunchKernelGGL(cast_memory, dim3(block), dim3(THREADS_PER_BLOCK), 0, 0, arr_out, arr, size);
        hipErrcheck(hipGetLastError());
        hipErrcheck(hipDeviceSynchronize());
        hipErrcheck(hipFree(arr));
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
        hipErrcheck(hipMemcpy(arr, arr_in, sizeof(FPTYPE_in) * size, hipMemcpyDeviceToHost));
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
    hipErrcheck(hipFree(arr));
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