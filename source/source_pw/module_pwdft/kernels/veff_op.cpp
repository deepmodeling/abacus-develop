#include "source_pw/module_pwdft/kernels/veff_op.h"

namespace hamilt {

template <typename FPTYPE>
struct veff_pw_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* dev, const int& size, std::complex<FPTYPE>* out, const FPTYPE* in)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for (int ir = 0; ir < size; ++ir)
        {
            out[ir] *= in[ir];
        }
    }

    void operator()(const base_device::DEVICE_CPU* dev,
                    const int& size,
                    std::complex<FPTYPE>* out,
                    std::complex<FPTYPE>* out1,
                    const std::complex<FPTYPE>* in)
    {

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int ir = 0; ir < size; ir++) 
        {
            const int base = ir * 4;
            auto sup = out[ir] * (in[base]) + out1[ir] * (in[base + 1]);
            auto sdown = out1[ir] * (in[base + 2]) + out[ir] * (in[base + 3]);
            out[ir] = sup;
            out1[ir] = sdown;
        }
    }
};

template struct veff_pw_op<float, base_device::DEVICE_CPU>;
template struct veff_pw_op<double, base_device::DEVICE_CPU>;

}  // namespace hamilt

