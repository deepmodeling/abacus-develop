#ifndef HPSI_NORM_OP_H
#define HPSI_NORM_OP_H
#include <complex>
#include "module_base/module_device/device.h"
namespace hamilt
{
template <typename FPTYPE, typename Device>
struct hpsi_norm_op
{
    /// @brief normalize hPsi with emin and emax
    ///
    /// Input Parameters
    /// \param dev : the type of computing device
    /// \param nbands : nbands
    /// \param npwk_max : max number of planewaves of all k points
    /// \param npwk : number of planewaves of current k point
    /// \param Ebar : (emin + emax) / 2
    /// \param DeltaE : (emax - emin) / 2
    /// \param hpsi_norm : hPsi
    /// \param psi_in : input psi
    /// Output Parameters
    /// \param tmhpsi : output array
    void operator()(const Device* dev,
                    const int& nbands,
                    const int& npwk_max,
                    const int& npwk,
                    const FPTYPE& Ebar,
                    const FPTYPE& DeltaE,
                    std::complex<FPTYPE>* hpsi_norm,
                    const std::complex<FPTYPE>* psi_in);
};
#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
// Partially specialize functor for base_device::GpuDevice.
template <typename FPTYPE>
struct hpsi_norm_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* dev,
                    const int& nbands,
                    const int& npwk_max,
                    const int& npwk,
                    const FPTYPE& Ebar,
                    const FPTYPE& DeltaE,
                    std::complex<FPTYPE>* hpsi_norm,
                    const std::complex<FPTYPE>* psi_in);
};
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
} // namespace hamilt
#endif