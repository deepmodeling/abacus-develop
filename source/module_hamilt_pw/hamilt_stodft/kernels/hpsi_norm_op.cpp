#include "hpsi_norm_op.h"

#include "module_base/module_device/device.h"
namespace hamilt
{
template <typename FPTYPE>
struct hpsi_norm_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* dev,
                    const int& nbands,
                    const int& npwk_max,
                    const int& npwk,
                    const FPTYPE& Ebar,
                    const FPTYPE& DeltaE,
                    std::complex<FPTYPE>* hpsi_norm,
                    const std::complex<FPTYPE>* psi_in)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int ib = 0; ib < nbands; ++ib)
        {
            const int ig0 = ib * npwk_max;
            for (int ig = 0; ig < npwk; ++ig)
            {
                hpsi_norm[ig + ig0] = (hpsi_norm[ig + ig0] - Ebar * psi_in[ig + ig0]) / DeltaE;
            }
        }
    }
};

// template struct hpsi_norm_op<float, base_device::DEVICE_CPU>;
template struct hpsi_norm_op<double, base_device::DEVICE_CPU>;

} // namespace hamilt