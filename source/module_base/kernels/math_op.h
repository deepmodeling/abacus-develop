#ifndef MODULE_BASE_MATH_MULTI_DEVICE_H
#define MODULE_BASE_MATH_MULTI_DEVICE_H

#include "module_psi/psi.h"
#include <complex>

namespace ModuleBase {

template <typename FPTYPE, typename Device>
struct cal_ylm_real_op {
    /// @brief YLM_REAL::Real spherical harmonics ylm(G) up to l=lmax
    /// Use Numerical recursive algorithm as given in Numerical Recipes
    ///
    /// Input Parameters
    /// @param ctx - which device this function runs on
    /// @param ng - number of problem size
    /// @param lmax - determined by lmax2
    /// @param SQRT2 - ModuleBase::SQRT2
    /// @param PI - ModuleBase::PI
    /// @param PI_HALF - ModuleBase::PI_HALF
    /// @param FOUR_PI - ModuleBase::FOUR_PI,
    /// @param SQRT_INVERSE_FOUR_PI - ModuleBase::SQRT_INVERSE_FOUR_PI,
    /// @param g - input array with size npw * 3, GlobalC::wf.get_1qvec_cartesian
    /// @param p - intermediate array
    ///
    /// Output Parameters
    /// @param ylm - output array
    void operator() (
        const Device *ctx,
        const int &ng,
        const int &lmax,
        const FPTYPE &SQRT2,
        const FPTYPE &PI,
        const FPTYPE &PI_HALF,
        const FPTYPE &FOUR_PI,
        const FPTYPE &SQRT_INVERSE_FOUR_PI,
        const FPTYPE *g,
        FPTYPE * p,
        FPTYPE * ylm);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
template <typename FPTYPE>
struct cal_ylm_real_op<FPTYPE, psi::DEVICE_GPU> {
    void operator() (
        const psi::DEVICE_GPU *ctx,
        const int &ng,
        const int &lmax,
        const FPTYPE &SQRT2,
        const FPTYPE &PI,
        const FPTYPE &PI_HALF,
        const FPTYPE &FOUR_PI,
        const FPTYPE &SQRT_INVERSE_FOUR_PI,
        const FPTYPE *g,
        FPTYPE * p,
        FPTYPE * ylm);
};

#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
}  // namespace ModuleBase
#endif //MODULE_BASE_MATH_MULTI_DEVICE_H