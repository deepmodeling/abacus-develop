// TODO: This is a temperary location for these functions.
// And will be moved to a global module(module base) later.
#ifndef MODULE_PSI_DEVICE_H_
#define MODULE_PSI_DEVICE_H_

#include <iostream>
#include "module_psi/include/types.h"

namespace psi{
namespace device {

template<typename Device> AbacusDevice_t get_device_type (const Device* dev);

} // end of namespace device
} // end of namespace psi

#endif  // MODULE_PSI_DEVICE_H_