
#include <complex>
#include <iostream>
#include "module_psi/psi.h"
#include "module_psi/include/types.h"
#include "module_psi/include/device.h"

namespace psi{

namespace device{

// functions used in custom ops
template<> AbacusDevice_t get_device_type <DEVICE_CPU> (const DEVICE_CPU* dev) {
    return CpuDevice;
}

#if ((defined __CUDA) || (defined __ROCM))
template<> AbacusDevice_t get_device_type <DEVICE_GPU> (const DEVICE_GPU* dev) {
    return GpuDevice;
}
#endif

} // end of namespace device
} // end of namespace psi