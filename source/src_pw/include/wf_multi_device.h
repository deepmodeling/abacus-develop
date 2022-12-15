#ifndef SRC_PW_WF_MULTI_DEVICE_H
#define SRC_PW_WF_MULTI_DEVICE_H

#include "module_psi/psi.h"
#include <complex>

namespace src_pw {

template <typename FPTYPE, typename Device>
struct cal_sk_op {
    void operator() (
        const Device* ctx,
        const int &ik,
        const int &ntype,
        const int &nx,
        const int &ny,
        const int &nz,
        const int &npw,
        const int &npwx,
        const int &fftny,
        const int &eigts1_nc,
        const int &eigts2_nc,
        const int &eigts3_nc,
        const int * atom_na,
        const int * igl2isz,
        const int * is2fftixy,
        const FPTYPE &TWO_PI,
        const FPTYPE *kvec_c,
        const FPTYPE *atom_tau,
        std::complex<FPTYPE> *eigts1,
        std::complex<FPTYPE> *eigts2,
        std::complex<FPTYPE> *eigts3,
        std::complex<FPTYPE> *sk);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
template <typename FPTYPE>
struct cal_sk_op<FPTYPE, psi::DEVICE_GPU> {
    void operator() (
        const psi::DEVICE_GPU *ctx,
        const int &ik,
        const int &ntype,
        const int &nx,
        const int &ny,
        const int &nz,
        const int &npw,
        const int &npwx,
        const int &fftny,
        const int &eigts1_nc,
        const int &eigts2_nc,
        const int &eigts3_nc,
        const int * atom_na,
        const int * igl2isz,
        const int * is2fftixy,
        const FPTYPE &TWO_PI,
        const FPTYPE *kvec_c,
        const FPTYPE *atom_tau,
        std::complex<FPTYPE> *eigts1,
        std::complex<FPTYPE> *eigts2,
        std::complex<FPTYPE> *eigts3,
        std::complex<FPTYPE> *sk);
};
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
}  // namespace src_pw
#endif //SRC_PW_WF_MULTI_DEVICE_H