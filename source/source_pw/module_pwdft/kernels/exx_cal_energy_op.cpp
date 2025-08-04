#include "source_pw/module_pwdft/kernels/exx_cal_energy_op.h"
#include "source_psi/psi.h"

namespace hamilt
{

// #ifdef _OPENMP
// #pragma omp parallel for reduction(+:Eexx_ik_real)
// #endif
// for (int ig = 0; ig < rhopw_dev->npw; ig++)
// {
//     int nks = wfcpw->nks;
//     int npw = rhopw_dev->npw;
//     int nk = nks / nk_fac;
//     Real Fac = pot[ik * nks * npw + iq * npw + ig];

// Eexx_ik_real += Fac * (density_recip[ig] * std::conj(density_recip[ig])).real()
//                 * wg_iqb_real / nqs * wg_ikb_real / kv->wk[ik];
// }

template <typename FPTYPE>
struct exx_cal_energy_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>
{
    using T = std::complex<FPTYPE>;
    FPTYPE operator()(const T *den, const FPTYPE *pot, FPTYPE scalar, int npw)
    {
        FPTYPE energy = 0.0;
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:energy)
        #endif
        for (int ig = 0; ig < npw; ++ig)
        {
            // Calculate the energy contribution from each reciprocal lattice vector
            energy += (den[ig] * std::conj(den[ig])).real() * pot[ig];
        }
        // Scale the energy by the scalar factor
        return scalar * energy;
    }
};

template struct exx_cal_energy_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct exx_cal_energy_op<std::complex<double>, base_device::DEVICE_CPU>;
} // namespace hamilt
