#ifndef MODULE_DEVICE_H_
#define MODULE_DEVICE_H_

#include "types.h"
#include <fstream>

#ifdef __MPI
#include "mpi.h"
#endif

namespace base_device
{

// struct CPU;
// struct GPU;

// enum AbacusDevice_t
// {
//     UnKnown,
//     CpuDevice,
//     GpuDevice
// };

template <typename Device>
base_device::AbacusDevice_t get_device_type(const Device* dev);

template <typename T>
std::string get_current_precision(const T* var);

namespace information
{

/**
 * @brief Get the device info object
 * for module_esolver
 */
std::string get_device_info(std::string device_flag);

/**
 * @brief Get the device kpar object
 * for module_io GlobalV::KPAR
 */
int get_device_kpar(const int& kpar, const int& bndpar);

/**
 * @brief Get the device flag object
 * for module_io PARAM.inp.device
 */
std::string get_device_flag(const std::string& device,
                            const std::string& basis_type);

#if __MPI
/**
 * @brief Get the rank of current node
 *        Note that GPU can only be binded with CPU in the same node
 * 
 * @return int 
 */
int get_node_rank();
int get_node_rank_with_mpi_shared(const MPI_Comm mpi_comm = MPI_COMM_WORLD);
int stringCmp(const void* a, const void* b);

#ifdef __CUDA

int set_device_by_rank(const MPI_Comm mpi_comm = MPI_COMM_WORLD);
#endif

#endif

template <typename Device>
void print_device_info(const Device* dev, std::ofstream& ofs_device)
{
    return;
}

template <typename Device>
void record_device_memory(const Device* dev, std::ofstream& ofs_device, std::string str, size_t size)
{
    return;
}

} // end of namespace information
} // end of namespace base_device

/**
 * @brief for compatibility with __CUDA_ARCH__ 600 and earlier
 *
 */
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
static __inline__ __device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do
    {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
        // Note: uses integer comparison to avoid hang in case of NaN (since NaN !=
        // NaN) } while (assumed != old);
    } while (assumed != old);
    return __longlong_as_double(old);
}
#endif

#endif // MODULE_DEVICE_H_