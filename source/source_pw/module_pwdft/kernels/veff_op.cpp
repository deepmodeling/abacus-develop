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
                    std::complex<FPTYPE>* in)
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

template<typename FPTYPE>
struct rearrange<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* dev,
                    const int& size, 
                    const FPTYPE* in, 
                    std::complex<FPTYPE>* out) const
    {
        for (int ir=0; ir < size; ir++)
        {
            const int base = 4 *ir;
            FPTYPE part_1 = in[ir];
            FPTYPE part_2 = in[ir + size];
            FPTYPE part_3 = in[ir + 2*size];
            FPTYPE part_4 = in[ir + 3*size];
            out[base ] = std::complex<FPTYPE>(part_1 + part_4, 0.0);
            out[base + 1] = std::complex<FPTYPE>(part_2 , -part_3);
            out[base + 2] = std::complex<FPTYPE>(part_1 - part_4, 0.0);
            out[base + 3] = std::complex<FPTYPE>(part_2, part_3);
        }
    }
};

template struct rearrange<float, base_device::DEVICE_CPU>;
template struct rearrange<double, base_device::DEVICE_CPU>;
template struct veff_pw_op<float, base_device::DEVICE_CPU>;
template struct veff_pw_op<double, base_device::DEVICE_CPU>;

}  // namespace hamilt

