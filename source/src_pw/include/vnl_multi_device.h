#ifndef SRC_PW_VNL_MULTI_DEVICE_H
#define SRC_PW_VNL_MULTI_DEVICE_H

#include "module_psi/psi.h"
#include <complex>

namespace src_pw {

template <typename FPTYPE, typename Device>
struct cal_vnl_op {
    void operator() (
        const Device* ctx,
        const int &ntype,
        const int &npw,
        const int &npwx,
        const int &nhm,
        const int &NQX,
        const int &tab_2,
        const int &tab_3,
        const int * atom_na,
        const int * atom_nb,
        const int * atom_nh,
        const FPTYPE &DQ,
        const FPTYPE &tpiba,
        const std::complex<FPTYPE> &NEG_IMAG_UNIT,
        const FPTYPE *gk,
        const FPTYPE *ylm,
        const FPTYPE *indv,
        const FPTYPE *nhtol,
        const FPTYPE *nhtolm,
        const FPTYPE *tab,
        FPTYPE *vkb1,
        const std::complex<FPTYPE> *sk,
        std::complex<FPTYPE> *vkb_in);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
template <typename FPTYPE>
struct cal_vnl_op<FPTYPE, psi::DEVICE_GPU> {
    void operator() (
        const psi::DEVICE_GPU *ctx,
        const int &ntype,
        const int &npw,
        const int &npwx,
        const int &nhm,
        const int &NQX,
        const int &tab_2,
        const int &tab_3,
        const int * atom_na,
        const int * atom_nb,
        const int * atom_nh,
        const FPTYPE &DQ,
        const FPTYPE &tpiba,
        const std::complex<FPTYPE> &NEG_IMAG_UNIT,
        const FPTYPE *gk,
        const FPTYPE *ylm,
        const FPTYPE *indv,
        const FPTYPE *nhtol,
        const FPTYPE *nhtolm,
        const FPTYPE *tab,
        FPTYPE *vkb1,
        const std::complex<FPTYPE> *sk,
        std::complex<FPTYPE> *vkb_in);
};
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
}  // namespace src_pw
#endif //SRC_PW_VNL_MULTI_DEVICE_H