#ifndef MODULE_BASE_MATH_MULTI_DEVICE_H
#define MODULE_BASE_MATH_MULTI_DEVICE_H

#include "module_psi/psi.h"
#include <complex>

namespace ModuleBase {

template <typename FPTYPE, typename Device>
struct cal_ylm_real_op {
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
}  // namespace src_pw
#endif //MODULE_BASE_MATH_MULTI_DEVICE_H