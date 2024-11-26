#ifndef MODULE_TYPES_H_
#define MODULE_TYPES_H_

namespace base_device
{

struct DEVICE_CPU;
struct DEVICE_GPU;
struct DEVICE_DSP;

enum AbacusDevice_t
{
    UnKnown,
    CpuDevice,
    GpuDevice,
    DspDevice
};

} // namespace base_device

#endif // MODULE_TYPES_H_