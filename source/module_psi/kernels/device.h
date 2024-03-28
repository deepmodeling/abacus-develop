// TODO: This is a temperary location for these functions.
// And will be moved to a global module(module base) later.
#ifndef MODULE_PSI_DEVICE_H_
#define MODULE_PSI_DEVICE_H_

#include "types.h"

#include <complex>
#include <iostream>

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
static __inline__ __device__ double atomicAdd(double *address, double val) {
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));
    // Note: uses integer comparison to avoid hang in case of NaN (since NaN !=
    // NaN) } while (assumed != old);
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif

namespace psi {
namespace device {

template<typename Device> AbacusDevice_t get_device_type (const Device* dev);

template<typename T> std::string get_current_precision(const T* var);

template<typename Device> void print_device_info (const Device* dev, std::ofstream& ofs_device) {return;}

template<typename Device> void record_device_memory (const Device* dev, std::ofstream& ofs_device, std::string str, size_t size) {return;}

std::string get_device_info(std::string device_flag);

int get_device_kpar(const int& kpar);
std::string get_device_flag(const std::string& device, const std::string& ks_solver, const std::string& basis_type, const bool& gamma_only);

#if __MPI
int get_node_rank();
int stringCmp(const void *a, const void* b);
#endif

#if ((defined __CUDA) || (defined __ROCM))
int get_device_num();
void set_device(const int rank);
#endif

} // end of namespace device
} // end of namespace psi

#endif  // MODULE_PSI_DEVICE_H_