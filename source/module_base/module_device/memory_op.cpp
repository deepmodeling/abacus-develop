#include "memory_op.h"

#include "device.h"
#include "module_base/memory.h"
#include "module_base/tool_threading.h"

#include <complex>
#include <iostream>
#include <string.h>

namespace device
{
namespace memory
{

template <typename FPTYPE>
struct resize_memory_op<FPTYPE, device::CPU>
{
    void operator()(const device::CPU* dev, FPTYPE*& arr, const size_t size, const char* record_in)
    {
        if (arr != nullptr)
        {
            free(arr);
        }
        arr = (FPTYPE*)malloc(sizeof(FPTYPE) * size);
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
            ModuleBase::Memory::record(record_string, sizeof(FPTYPE) * size);
        }
    }
};

template <typename FPTYPE>
struct set_memory_op<FPTYPE, device::CPU>
{
    void operator()(const device::CPU* dev, FPTYPE* arr, const int var, const size_t size)
    {
        ModuleBase::OMP_PARALLEL([&](int num_thread, int thread_id) {
            int beg = 0, len = 0;
            ModuleBase::BLOCK_TASK_DIST_1D(num_thread, thread_id, size, (size_t)4096 / sizeof(FPTYPE), beg, len);
            memset(arr + beg, var, sizeof(FPTYPE) * len);
        });
    }
};

template <typename FPTYPE>
struct synchronize_memory_op<FPTYPE, device::CPU, device::CPU>
{
    void operator()(const device::CPU* dev_out,
                    const device::CPU* dev_in,
                    FPTYPE* arr_out,
                    const FPTYPE* arr_in,
                    const size_t size)
    {
        ModuleBase::OMP_PARALLEL([&](int num_thread, int thread_id) {
            int beg = 0, len = 0;
            ModuleBase::BLOCK_TASK_DIST_1D(num_thread, thread_id, size, (size_t)4096 / sizeof(FPTYPE), beg, len);
            memcpy(arr_out + beg, arr_in + beg, sizeof(FPTYPE) * len);
        });
    }
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, device::CPU, device::CPU>
{
    void operator()(const device::CPU* dev_out,
                    const device::CPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE_out))
#endif
        for (int ii = 0; ii < size; ii++)
        {
            arr_out[ii] = static_cast<FPTYPE_out>(arr_in[ii]);
        }
    }
};

template <typename FPTYPE>
struct delete_memory_op<FPTYPE, device::CPU>
{
    void operator()(const device::CPU* dev, FPTYPE* arr)
    {
        free(arr);
    }
};

template struct resize_memory_op<int, device::CPU>;
template struct resize_memory_op<float, device::CPU>;
template struct resize_memory_op<double, device::CPU>;
template struct resize_memory_op<std::complex<float>, device::CPU>;
template struct resize_memory_op<std::complex<double>, device::CPU>;

template struct set_memory_op<int, device::CPU>;
template struct set_memory_op<float, device::CPU>;
template struct set_memory_op<double, device::CPU>;
template struct set_memory_op<std::complex<float>, device::CPU>;
template struct set_memory_op<std::complex<double>, device::CPU>;

template struct synchronize_memory_op<int, device::CPU, device::CPU>;
template struct synchronize_memory_op<float, device::CPU, device::CPU>;
template struct synchronize_memory_op<double, device::CPU, device::CPU>;
template struct synchronize_memory_op<std::complex<float>, device::CPU, device::CPU>;
template struct synchronize_memory_op<std::complex<double>, device::CPU, device::CPU>;

template struct cast_memory_op<float, float, device::CPU, device::CPU>;
template struct cast_memory_op<double, double, device::CPU, device::CPU>;
template struct cast_memory_op<float, double, device::CPU, device::CPU>;
template struct cast_memory_op<double, float, device::CPU, device::CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, device::CPU, device::CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, device::CPU, device::CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, device::CPU, device::CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, device::CPU, device::CPU>;

template struct delete_memory_op<int, device::CPU>;
template struct delete_memory_op<float, device::CPU>;
template struct delete_memory_op<double, device::CPU>;
template struct delete_memory_op<std::complex<float>, device::CPU>;
template struct delete_memory_op<std::complex<double>, device::CPU>;

#if !(defined(__CUDA) || defined(__ROCM))
// template <typename FPTYPE>
// struct resize_memory_op<FPTYPE, device::GPU>
// {
//     void operator()(const device::GPU* dev, FPTYPE*& arr, const size_t size, const char* record_in = nullptr)
//     {
//     }
// };

// template <typename FPTYPE>
// struct set_memory_op<FPTYPE, device::GPU>
// {
//     void operator()(const device::GPU* dev, FPTYPE* arr, const int var, const size_t size)
//     {
//     }
// };

// template <typename FPTYPE>
// struct synchronize_memory_op<FPTYPE, device::GPU, device::GPU>
// {
//     void operator()(const device::GPU* dev_out,
//                     const device::GPU* dev_in,
//                     FPTYPE* arr_out,
//                     const FPTYPE* arr_in,
//                     const size_t size)
//     {
//     }
// };

// template <typename FPTYPE>
// struct synchronize_memory_op<FPTYPE, device::GPU, device::CPU>
// {
//     void operator()(const device::GPU* dev_out,
//                     const device::CPU* dev_in,
//                     FPTYPE* arr_out,
//                     const FPTYPE* arr_in,
//                     const size_t size)
//     {
//     }
// };

// template <typename FPTYPE>
// struct synchronize_memory_op<FPTYPE, device::CPU, device::GPU>
// {
//     void operator()(const device::CPU* dev_out,
//                     const device::GPU* dev_in,
//                     FPTYPE* arr_out,
//                     const FPTYPE* arr_in,
//                     const size_t size)
//     {
//     }
// };

// template <typename FPTYPE_out, typename FPTYPE_in>
// struct cast_memory_op<FPTYPE_out, FPTYPE_in, device::GPU, device::GPU>
// {
//     void operator()(const device::GPU* dev_out,
//                     const device::GPU* dev_in,
//                     FPTYPE_out* arr_out,
//                     const FPTYPE_in* arr_in,
//                     const size_t size)
//     {
//     }
// };

// template <typename FPTYPE_out, typename FPTYPE_in>
// struct cast_memory_op<FPTYPE_out, FPTYPE_in, device::GPU, device::CPU>
// {
//     void operator()(const device::GPU* dev_out,
//                     const device::CPU* dev_in,
//                     FPTYPE_out* arr_out,
//                     const FPTYPE_in* arr_in,
//                     const size_t size)
//     {
//     }
// };

// template <typename FPTYPE_out, typename FPTYPE_in>
// struct cast_memory_op<FPTYPE_out, FPTYPE_in, device::CPU, device::GPU>
// {
//     void operator()(const device::CPU* dev_out,
//                     const device::GPU* dev_in,
//                     FPTYPE_out* arr_out,
//                     const FPTYPE_in* arr_in,
//                     const size_t size)
//     {
//     }
// };

// template <typename FPTYPE>
// struct delete_memory_op<FPTYPE, device::GPU>
// {
//     void operator()(const device::GPU* dev, FPTYPE* arr)
//     {
//     }
// };

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
#endif

} // namespace memory
} // namespace device