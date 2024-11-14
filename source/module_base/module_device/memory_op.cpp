#include "memory_op.h"

#include "module_base/memory.h"
#include "module_base/tool_threading.h"
#ifdef __DSP
#include "module_base/kernels/dsp/dsp_connector.h"
#include "module_base/global_variable.h"
#endif

#include <complex>
#include <cstring>

namespace base_device
{
namespace memory
{

template <typename FPTYPE>
struct resize_memory_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* dev, FPTYPE*& arr, const size_t size, const char* record_in)
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
struct set_memory_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* dev, FPTYPE* arr, const int var, const size_t size)
    {
        ModuleBase::OMP_PARALLEL([&](int num_thread, int thread_id) {
            int beg = 0, len = 0;
            ModuleBase::BLOCK_TASK_DIST_1D(num_thread, thread_id, size, (size_t)4096 / sizeof(FPTYPE), beg, len);
            memset(arr + beg, var, sizeof(FPTYPE) * len);
        });
    }
};

template <typename FPTYPE>
struct synchronize_memory_op<FPTYPE, base_device::DEVICE_CPU, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* dev_out,
                    const base_device::DEVICE_CPU* dev_in,
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
struct cast_memory_op<FPTYPE_out, FPTYPE_in, base_device::DEVICE_CPU, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* dev_out,
                    const base_device::DEVICE_CPU* dev_in,
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
struct delete_memory_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* dev, FPTYPE* arr)
    {
        free(arr);
    }
};

template struct resize_memory_op<int, base_device::DEVICE_CPU>;
template struct resize_memory_op<float, base_device::DEVICE_CPU>;
template struct resize_memory_op<double, base_device::DEVICE_CPU>;
template struct resize_memory_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct resize_memory_op<std::complex<double>, base_device::DEVICE_CPU>;

template struct set_memory_op<int, base_device::DEVICE_CPU>;
template struct set_memory_op<float, base_device::DEVICE_CPU>;
template struct set_memory_op<double, base_device::DEVICE_CPU>;
template struct set_memory_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct set_memory_op<std::complex<double>, base_device::DEVICE_CPU>;

template struct synchronize_memory_op<int, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
template struct synchronize_memory_op<float, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
template struct synchronize_memory_op<double, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<float>, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<double>, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;

template struct cast_memory_op<float, float, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
template struct cast_memory_op<double, double, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
template struct cast_memory_op<float, double, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
template struct cast_memory_op<double, float, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>,
                               std::complex<float>,
                               base_device::DEVICE_CPU,
                               base_device::DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>,
                               std::complex<double>,
                               base_device::DEVICE_CPU,
                               base_device::DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>,
                               std::complex<double>,
                               base_device::DEVICE_CPU,
                               base_device::DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>,
                               std::complex<float>,
                               base_device::DEVICE_CPU,
                               base_device::DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>, float, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>, double, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;

template struct delete_memory_op<int, base_device::DEVICE_CPU>;
template struct delete_memory_op<float, base_device::DEVICE_CPU>;
template struct delete_memory_op<double, base_device::DEVICE_CPU>;
template struct delete_memory_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct delete_memory_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct delete_memory_op<float*, base_device::DEVICE_CPU>;
template struct delete_memory_op<double*, base_device::DEVICE_CPU>;
template struct delete_memory_op<std::complex<float>*, base_device::DEVICE_CPU>;
template struct delete_memory_op<std::complex<double>*, base_device::DEVICE_CPU>;

#if !(defined(__CUDA) || defined(__ROCM))

template <typename FPTYPE>
struct resize_memory_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* dev,
                    FPTYPE*& arr,
                    const size_t size,
                    const char* record_in = nullptr)
    {
    }
};

template <typename FPTYPE>
struct set_memory_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* dev, FPTYPE* arr, const int var, const size_t size)
    {
    }
};

template <typename FPTYPE>
struct synchronize_memory_op<FPTYPE, base_device::DEVICE_GPU, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* dev_out,
                    const base_device::DEVICE_GPU* dev_in,
                    FPTYPE* arr_out,
                    const FPTYPE* arr_in,
                    const size_t size)
    {
    }
};

template <typename FPTYPE>
struct synchronize_memory_op<FPTYPE, base_device::DEVICE_GPU, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_GPU* dev_out,
                    const base_device::DEVICE_CPU* dev_in,
                    FPTYPE* arr_out,
                    const FPTYPE* arr_in,
                    const size_t size)
    {
    }
};

template <typename FPTYPE>
struct synchronize_memory_op<FPTYPE, base_device::DEVICE_CPU, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_CPU* dev_out,
                    const base_device::DEVICE_GPU* dev_in,
                    FPTYPE* arr_out,
                    const FPTYPE* arr_in,
                    const size_t size)
    {
    }
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, base_device::DEVICE_GPU, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* dev_out,
                    const base_device::DEVICE_GPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size)
    {
    }
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, base_device::DEVICE_GPU, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_GPU* dev_out,
                    const base_device::DEVICE_CPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size)
    {
    }
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, base_device::DEVICE_CPU, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_CPU* dev_out,
                    const base_device::DEVICE_GPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size)
    {
    }
};

template <typename FPTYPE>
struct delete_memory_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* dev, FPTYPE* arr)
    {
    }
};

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
#endif

#ifdef __DSP

template <typename FPTYPE>
struct resize_memory_op_mt<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* dev, FPTYPE*& arr, const size_t size, const char* record_in)
    {
        if (arr != nullptr)
        {
            free_ht(arr);
        }
        arr = (FPTYPE*)malloc_ht(sizeof(FPTYPE) * size, GlobalV::MY_RANK);
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
struct delete_memory_op_mt<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* dev, FPTYPE* arr)
    {
        free_ht(arr);
    }
};


template struct resize_memory_op_mt<int, base_device::DEVICE_CPU>;
template struct resize_memory_op_mt<float, base_device::DEVICE_CPU>;
template struct resize_memory_op_mt<double, base_device::DEVICE_CPU>;
template struct resize_memory_op_mt<std::complex<float>, base_device::DEVICE_CPU>;
template struct resize_memory_op_mt<std::complex<double>, base_device::DEVICE_CPU>;

template struct delete_memory_op_mt<int, base_device::DEVICE_CPU>;
template struct delete_memory_op_mt<float, base_device::DEVICE_CPU>;
template struct delete_memory_op_mt<double, base_device::DEVICE_CPU>;
template struct delete_memory_op_mt<std::complex<float>, base_device::DEVICE_CPU>;
template struct delete_memory_op_mt<std::complex<double>, base_device::DEVICE_CPU>;
#endif

} // namespace memory
} // namespace base_device