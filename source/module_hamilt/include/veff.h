#ifndef MODULE_HAMILT_VEFF_H
#define MODULE_HAMILT_VEFF_H

#include "module_psi/psi.h"
#include <complex>

namespace hamilt {
template <typename FPTYPE, typename Device>
struct veff_pw_op {
    void operator() (
        const Device* dev,
        const int& size,
        std::complex<FPTYPE>* out,
        const FPTYPE* in);
    
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