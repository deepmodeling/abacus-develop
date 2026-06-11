#include "source_pw/module_pwdft/kernels/exx_stress_op.h"

namespace hamilt
{
template <typename FPTYPE>
struct exx_stress_accumulate_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>
{
    using T = std::complex<FPTYPE>;
    using Real = FPTYPE;
    void operator()(const T* density_recip,
                    const Real* pot,
                    const Real* pot_stress,
                    const Real* gcar,
                    Real dkx,
                    Real dky,
                    Real dkz,
                    Real tpiba,
                    Real scalar,
                    int npw,
                    Real* sigma)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            const Real density_recip2 = std::real(density_recip[ig] * std::conj(density_recip[ig]));
            const Real kqg[3] = {(dkx + gcar[ig * 3]) * tpiba,
                                 (dky + gcar[ig * 3 + 1]) * tpiba,
                                 (dkz + gcar[ig * 3 + 2]) * tpiba};
            int idx = 0;
            for (int alpha = 0; alpha < 3; ++alpha)
            {
                for (int beta = alpha; beta < 3; ++beta)
                {
                    const Real delta_ab = alpha == beta ? Real(1.0) : Real(0.0);
                    sigma[idx++] += scalar * density_recip2 * pot[ig]
                                    * (kqg[alpha] * kqg[beta] * pot_stress[ig] - delta_ab);
                }
            }
        }
    }
};

template struct exx_stress_accumulate_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct exx_stress_accumulate_op<std::complex<double>, base_device::DEVICE_CPU>;
} // namespace hamilt
