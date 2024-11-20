#ifndef __PARALLEL_DEVICE_H__
#define __PARALLEL_DEVICE_H__
#ifdef __MPI
#include "mpi.h"
#include "module_base/module_device/device.h"
#include "module_base/module_device/memory_op.h"
#include <complex>
namespace Parallel_Common
{
void bcast_data(std::complex<double>* object, const int& n, const MPI_Comm& comm);
void bcast_data(std::complex<float>* object, const int& n, const MPI_Comm& comm);
void bcast_data(double* object, const int& n, const MPI_Comm& comm);
void bcast_data(float* object, const int& n, const MPI_Comm& comm);
void reduce_data(std::complex<double>* object, const int& n, const MPI_Comm& comm);
void reduce_data(std::complex<float>* object, const int& n, const MPI_Comm& comm);
void reduce_data(double* object, const int& n, const MPI_Comm& comm);
void reduce_data(float* object, const int& n, const MPI_Comm& comm);

/**
 * @brief bcast data in Device
 * 
 * @tparam T: float, double, std::complex<float>, std::complex<double>
 * @tparam Device 
 * @param ctx Device ctx
 * @param object complex arrays in Device
 * @param n the size of complex arrays
 * @param comm MPI_Comm
 * @param tmp_space tmp space in CPU
 */
template <typename T, typename Device>
void bcast_dev(const Device* ctx, T* object, const int& n, const MPI_Comm& comm, T* tmp_space = nullptr)
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

    bcast_data(object_cpu, n, comm);

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
void reduce_dev(const Device* ctx, T* object, const int& n, const MPI_Comm& comm, T* tmp_space = nullptr)
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

    reduce_data(object_cpu, n, comm);

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
#endif