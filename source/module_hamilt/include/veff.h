#ifndef MODULE_HAMILT_VEFF_H
#define MODULE_HAMILT_VEFF_H

#include "module_psi/psi.h"
#include <complex>

namespace hamilt {
template <typename FPTYPE, typename Device>
struct veff_pw_op {
    /// @brief Compute the effective potential of hPsi in real space,
    /// out[ir] *= in[ir];
    ///
    /// Input Parameters
    /// \param dev : the type of computing device
    /// \param size : array size
    /// \param in : input array, elecstate::Potential::v_effective
    ///
    /// Output Parameters
    /// \param out : output array
    void operator() (
        const Device* dev,
        const int& size,
        std::complex<FPTYPE>* out,
        const FPTYPE* in);

    /// @brief Compute the effective potential of hPsi in real space with NSPIN > 2,
    ///
    /// out[ir] = out[ir] * (in[0][ir] + in[3][ir])
    ///       + out1[ir]
    ///               * (in[1][ir]
    ///               - std::complex<FPTYPE>(0.0, 1.0) * in[2][ir]);
    ///
    /// out1[ir] = out1[ir] * (in[0][ir] - in[3][ir])
    ///         + out[ir]
    ///             * (in[1][ir]
    ///                 + std::complex<FPTYPE>(0.0, 1.0) * in[2][ir]);
    ///
    /// Input Parameters
    /// \param dev : the type of computing device
    /// \param size : array size
    /// \param in : input array, elecstate::Potential::v_effective
    ///
    /// Output Parameters
    /// \param out : output array 1
    /// \param out1 : output array 2
    void operator() (
        const Device* dev,
        const int& size,
        std::complex<FPTYPE>* out,
        std::complex<FPTYPE>* out1,
        const FPTYPE** in);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
// Partially specialize functor for psi::GpuDevice.
template <typename FPTYPE>
struct veff_pw_op<FPTYPE, psi::DEVICE_GPU> {
    void operator() (
        const psi::DEVICE_GPU* dev,
        const int& size,
        std::complex<FPTYPE>* out,
        const FPTYPE* in);
    
    void operator() (
        const psi::DEVICE_GPU* dev,
        const int& size,
        std::complex<FPTYPE>* out,
        std::complex<FPTYPE>* out1,
        const FPTYPE** in);
};
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
} // namespace hamilt
#endif //MODULE_HAMILT_VEFF_H