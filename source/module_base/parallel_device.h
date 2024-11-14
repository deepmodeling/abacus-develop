#ifdef __MPI
#include "mpi.h"
#include "module_base/module_device/device.h"
#include <complex>
#include <string>
#include <vector>
namespace Parallel_Common
{
void bcast_complex(std::complex<double>* object, const int& n, const MPI_Comm& comm)
{
    MPI_Bcast(object, n * 2, MPI_DOUBLE, 0, comm);
}
void bcast_complex(std::complex<float>* object, const int& n, const MPI_Comm& comm)
{
    MPI_Bcast(object, n * 2, MPI_FLOAT, 0, comm);
}
void bcast_real(double* object, const int& n, const MPI_Comm& comm)
{
    MPI_Bcast(object, n, MPI_DOUBLE, 0, comm);
}
void bcast_real(float* object, const int& n, const MPI_Comm& comm)
{
    MPI_Bcast(object, n, MPI_FLOAT, 0, comm);
}

template <typename T, typename Device>
/**
 * @brief bcast complex in Device
 * 
 * @param ctx Device ctx
 * @param object complex arrays in Device
 * @param n the size of complex arrays
 * @param comm MPI_Comm
 * @param tmp_space tmp space in CPU
 */
void bcast_complex(const Device* ctx, T* object, const int& n, const MPI_Comm& comm, T* tmp_space = nullptr)
{
    const base_device::DEVICE_CPU* cpu_ctx = {};
    T* object_cpu = nullptr;
    bool alloc = false;
    if (base_device::get_device_type<Device>(ctx) == base_device::GpuDevice)
    {
        if(tmp_space == nullptr)
        {
            base_device::memory::resize_memory_op<T, base_device::DEVICE_CPU>()(cpu_ctx, object_cpu, n);
            alloc = true;
        }
        else
        {
            object_cpu = tmp_space;
        }
        base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>()(cpu_ctx, ctx, object_cpu, object, n);
    }
    else
    {
        object_cpu = object;
    }

    bcast_complex(object_cpu, n, comm);

    if (base_device::get_device_type<Device>(ctx) == base_device::GpuDevice)
    {
        base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(ctx, cpu_ctx, object, object_cpu, n);
        if(alloc)
        {
            base_device::memory::delete_memory_op<T, base_device::DEVICE_CPU>()(cpu_ctx, object_cpu);
        }
    }
    return;
}

template <typename T, typename Device>
void bcast_real(const Device* ctx, T* object, const int& n, const MPI_Comm& comm, T* tmp_space = nullptr)
{
    const base_device::DEVICE_CPU* cpu_ctx = {};
    T* object_cpu = nullptr;
    bool alloc = false;
    if (base_device::get_device_type<Device>(ctx) == base_device::GpuDevice)
    {
        if(tmp_space == nullptr)
        {
            base_device::memory::resize_memory_op<T, base_device::DEVICE_CPU>()(cpu_ctx, object_cpu, n);
            alloc = true;
        }
        else
        {
            object_cpu = tmp_space;
        }
        base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>()(cpu_ctx, ctx, object_cpu, object, n);
    }
    else
    {
        object_cpu = object;
    }

    bcast_real(object_cpu, n, comm);

    if (base_device::get_device_type<Device>(ctx) == base_device::GpuDevice)
    {
        base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(ctx, cpu_ctx, object, object_cpu, n);
        if(alloc)
        {
            base_device::memory::delete_memory_op<T, base_device::DEVICE_CPU>()(cpu_ctx, object_cpu);
        }
    }
    return;
}
}
    

#endif