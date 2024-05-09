#ifndef MODULE_TYPES_H_
#define MODULE_TYPES_H_

namespace psi
{

struct DEVICE_CPU;
struct DEVICE_GPU;

enum AbacusDevice_t
{
    UnKnown,
    CpuDevice,
    GpuDevice,
    SyclDevice
};

} // namespace psi

#endif // MODULE_TYPES_H_