#include "module_hamilt_pw/hamilt_pwdft/kernels/veff_op.h"

namespace hamilt {

template <typename FPTYPE>
struct veff_pw_op<FPTYPE, psi::DEVICE_CPU> {
    void operator() (
        const psi::DEVICE_CPU* dev,
        const int& size,
        std::complex<FPTYPE>* out,
        const FPTYPE* in)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for (int ir = 0; ir < size; ++ir)
        {
            out[ir] *= in[ir];
        }
    }

    void operator() (
        const psi::DEVICE_CPU* dev,
        const int& size,
        std::complex<FPTYPE>* out,
        std::complex<FPTYPE>* out1,
        const FPTYPE** in)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int ir = 0; ir < size; ir++) {
            auto sup = out[ir] * (in[0][ir] + in[3][ir])
                + out1[ir]
                        * (in[1][ir]
                        - std::complex<FPTYPE>(0.0, 1.0) * in[2][ir]);
            auto sdown = out1[ir] * (in[0][ir] - in[3][ir])
                    + out[ir]
                        * (in[1][ir]
                            + std::complex<FPTYPE>(0.0, 1.0) * in[2][ir]);
            out[ir] = sup;
            out1[ir] = sdown;
        }
    }
};

template struct veff_pw_op<float, psi::DEVICE_CPU>;
template struct veff_pw_op<double, psi::DEVICE_CPU>;

}  // namespace hamilt

