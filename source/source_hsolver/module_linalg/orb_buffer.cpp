#include "source_hsolver/module_linalg/orb_buffer.h"

#include <complex>

namespace hsolver
{
template <typename T, typename Device>
void zero_pad(T* data, const OrbLayout& layout)
{
    const int stride = layout.ld / layout.npol;
    if (stride > layout.npw && layout.bands > 0)
    {
        T* tail = data + layout.npw;
        const size_t width = stride - layout.npw;
        const size_t columns = static_cast<size_t>(layout.bands) * layout.npol;
        base_device::memory::set_memory_2d_op<T, Device>()(tail, stride, 0, width, columns);
    }
}

template <typename T, typename Device>
const T* clean_orbs(const T* input, const OrbLayout& layout, ct::Tensor* scratch)
{
    if (layout.npol == 1)
    {
        return input;
    }
    T* data = scratch->template data<T>();
    const size_t count = static_cast<size_t>(layout.ld) * layout.bands;
    if (count > 0)
    {
        base_device::memory::synchronize_memory_op<T, Device, Device>()(data, input, count);
    }
    zero_pad<T, Device>(data, layout);
    return data;
}

template void zero_pad<std::complex<float>, base_device::DEVICE_CPU>(std::complex<float>*, const OrbLayout&);
template void zero_pad<std::complex<double>, base_device::DEVICE_CPU>(std::complex<double>*, const OrbLayout&);
template const std::complex<float>* clean_orbs<std::complex<float>, base_device::DEVICE_CPU>(const std::complex<float>*,
                                                                                             const OrbLayout&,
                                                                                             ct::Tensor*);
template const std::complex<double>* clean_orbs<std::complex<double>, base_device::DEVICE_CPU>(const std::complex<double>*,
                                                                                               const OrbLayout&,
                                                                                               ct::Tensor*);
#if defined(__CUDA) || defined(__ROCM)
template void zero_pad<std::complex<float>, base_device::DEVICE_GPU>(std::complex<float>*, const OrbLayout&);
template void zero_pad<std::complex<double>, base_device::DEVICE_GPU>(std::complex<double>*, const OrbLayout&);
template const std::complex<float>* clean_orbs<std::complex<float>, base_device::DEVICE_GPU>(const std::complex<float>*,
                                                                                             const OrbLayout&,
                                                                                             ct::Tensor*);
template const std::complex<double>* clean_orbs<std::complex<double>, base_device::DEVICE_GPU>(const std::complex<double>*,
                                                                                               const OrbLayout&,
                                                                                               ct::Tensor*);
#endif
} // namespace hsolver
