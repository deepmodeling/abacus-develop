// TODO: This is a temperary location for these functions.
// And will be moved to a global module(module base) later.
#ifndef MODULE_PSI_DEVICE_H_
#define MODULE_PSI_DEVICE_H_

#include <iostream>
#include "module_psi/include/types.h"

namespace psi{
namespace device {

template<typename Device> AbacusDevice_t get_device_type (const Device* dev);

template<typename Device> void print_device_info (const Device* dev, std::ofstream& ofs_device) {return;}

template<typename Device> void record_device_memory (const Device* dev, std::ofstream& ofs_device, std::string str, size_t size) {return;}

int get_device_kpar(const int& kpar);
std::string get_device_flag(const std::string& device, const std::string& ks_solver, const std::string& basis_type);

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