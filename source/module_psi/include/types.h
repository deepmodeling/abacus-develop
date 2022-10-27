// TODO: This is a temperary location for these functions.
// And will be moved to a global module(module base) later.
#ifndef MODULE_PSI_TYPES_H_
#define MODULE_PSI_TYPES_H_

namespace psi {

struct DEVICE_CPU;
struct DEVICE_GPU;
struct DEVICE_SYCL;

enum AbacusDevice_t {UnKnown, CpuDevice, GpuDevice, SyclDevice};

}  // end of namespace psi

#endif  // MODULE_PSI_TYPES_H_