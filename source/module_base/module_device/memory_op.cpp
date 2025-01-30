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
    void operator()(FPTYPE*& arr, const size_t size, const char* record_in)
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
    void operator()(FPTYPE* arr, const int var, const size_t size)
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
    void operator()(FPTYPE* arr_out,
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
    void operator()(FPTYPE_out* arr_out,
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
    void operator()(FPTYPE* arr)
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
    void operator()(FPTYPE*& arr,
                    const size_t size,
                    const char* record_in = nullptr)
    {
    }
};

template <typename FPTYPE>
struct set_memory_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(FPTYPE* arr, const int var, const size_t size)
    {
    }
};

template <typename FPTYPE>
struct synchronize_memory_op<FPTYPE, base_device::DEVICE_GPU, base_device::DEVICE_GPU>
{
    void operator()(FPTYPE* arr_out,
                    const FPTYPE* arr_in,
                    const size_t size)
    {
    }
};

template <typename FPTYPE>
struct synchronize_memory_op<FPTYPE, base_device::DEVICE_GPU, base_device::DEVICE_CPU>
{
    void operator()(FPTYPE* arr_out,
                    const FPTYPE* arr_in,
                    const size_t size)
    {
    }
};

template <typename FPTYPE>
struct synchronize_memory_op<FPTYPE, base_device::DEVICE_CPU, base_device::DEVICE_GPU>
{
    void operator()(FPTYPE* arr_out,
                    const FPTYPE* arr_in,
                    const size_t size)
    {
    }
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, base_device::DEVICE_GPU, base_device::DEVICE_GPU>
{
    void operator()(FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size)
    {
    }
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, base_device::DEVICE_GPU, base_device::DEVICE_CPU>
{
    void operator()(FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size)
    {
    }
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, base_device::DEVICE_CPU, base_device::DEVICE_GPU>
{
    void operator()(FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size)
    {
    }
};

template <typename FPTYPE>
struct delete_memory_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(FPTYPE* arr)
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
    void operator()(FPTYPE*& arr, const size_t size, const char* record_in)
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
    void operator()(FPTYPE* arr)
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

template <typename FPTYPE>
void resize_memory(FPTYPE* arr, const size_t size, base_device::AbacusDevice_t device_type)
{
    if (device_type == base_device::AbacusDevice_t::CpuDevice){
        resize_memory_op<FPTYPE, base_device::DEVICE_CPU>()(arr, size);
    }
    else if (device_type == base_device::AbacusDevice_t::GpuDevice){
        resize_memory_op<FPTYPE, base_device::DEVICE_GPU>()(arr, size);
    }
}

template <typename FPTYPE>
void set_memory(FPTYPE* arr, const int var, const size_t size, base_device::AbacusDevice_t device_type){
    if (device_type == base_device::AbacusDevice_t::CpuDevice){
        set_memory_op<FPTYPE, base_device::DEVICE_CPU>()(arr, var, size);
    }
    else if (device_type == base_device::AbacusDevice_t::GpuDevice){
        set_memory_op<FPTYPE, base_device::DEVICE_GPU>()(arr, var, size);
    }
}

template <typename FPTYPE>
void synchronize_memory(FPTYPE* arr_out, const FPTYPE* arr_in, const size_t size, base_device::AbacusDevice_t device_type_out, base_device::AbacusDevice_t device_type_in){
    if (device_type_out == base_device::AbacusDevice_t::CpuDevice || device_type_in == base_device::AbacusDevice_t::CpuDevice){
        synchronize_memory_op<FPTYPE, DEVICE_CPU, DEVICE_CPU>()(arr_out, arr_in, size);
    }
    else if (device_type_out == base_device::AbacusDevice_t::CpuDevice || device_type_in == base_device::AbacusDevice_t::GpuDevice){
        synchronize_memory_op<FPTYPE, DEVICE_CPU, DEVICE_GPU>()(arr_out, arr_in, size);
    }
    else if (device_type_out == base_device::AbacusDevice_t::GpuDevice || device_type_in == base_device::AbacusDevice_t::CpuDevice){
        synchronize_memory_op<FPTYPE, DEVICE_GPU, DEVICE_CPU>()(arr_out, arr_in, size);
    }
    else if (device_type_out == base_device::AbacusDevice_t::GpuDevice || device_type_in == base_device::AbacusDevice_t::GpuDevice){
        synchronize_memory_op<FPTYPE, DEVICE_GPU, DEVICE_GPU>()(arr_out, arr_in, size);
    }
}

template <typename FPTYPE_out, typename FPTYPE_in>
void cast_memory(FPTYPE_out* arr_out, const FPTYPE_in* arr_in, const size_t size, base_device::AbacusDevice_t device_type_out, base_device::AbacusDevice_t device_type_in)
{
    if (device_type_out == base_device::AbacusDevice_t::CpuDevice || device_type_in == base_device::AbacusDevice_t::CpuDevice){
        cast_memory_op<FPTYPE_out, FPTYPE_in, DEVICE_CPU, DEVICE_CPU>()(arr_out, arr_in, size);
    }
    else if (device_type_out == base_device::AbacusDevice_t::CpuDevice || device_type_in == base_device::AbacusDevice_t::GpuDevice){
        cast_memory_op<FPTYPE_out, FPTYPE_in, DEVICE_CPU, DEVICE_GPU>()(arr_out, arr_in, size);
    }
    else if (device_type_out == base_device::AbacusDevice_t::GpuDevice || device_type_in == base_device::AbacusDevice_t::CpuDevice){
        cast_memory_op<FPTYPE_out, FPTYPE_in, DEVICE_GPU, DEVICE_CPU>()(arr_out, arr_in, size);
    }
    else if (device_type_out == base_device::AbacusDevice_t::GpuDevice || device_type_in == base_device::AbacusDevice_t::GpuDevice){
        cast_memory_op<FPTYPE_out, FPTYPE_in, DEVICE_GPU, DEVICE_GPU>()(arr_out, arr_in, size);
    }
}

template <typename FPTYPE>
void delete_memory(FPTYPE* arr, base_device::AbacusDevice_t device_type)
{
    if (device_type == base_device::AbacusDevice_t::CpuDevice){
        delete_memory_op<FPTYPE, DEVICE_CPU>()(arr);
    }
    else if (device_type == base_device::AbacusDevice_t::GpuDevice){
        delete_memory_op<FPTYPE, DEVICE_GPU>()(arr);
    }
}

} // namespace memory
} // namespace base_device