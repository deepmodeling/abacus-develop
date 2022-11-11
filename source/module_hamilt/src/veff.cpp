#include "module_hamilt/include/veff.h"

namespace hamilt{

template <typename FPTYPE>
struct veff_pw_op<FPTYPE, psi::DEVICE_CPU> {
    void operator() (
        const psi::DEVICE_CPU* dev,
        const int& size,
        std::complex<FPTYPE>* out,
        const FPTYPE* in)
    {
        for (int ir = 0; ir < size; ++ir)
        {
            out[ir] *= in[ir];
        }
    }
};

template struct veff_pw_op<double, psi::DEVICE_CPU>;

}  // namespace hamilt

