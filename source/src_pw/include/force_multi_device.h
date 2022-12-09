#ifndef SRC_PW_FORCE_MULTI_DEVICE_H
#define SRC_PW_FORCE_MULTI_DEVICE_H

#include "module_psi/psi.h"
#include <complex>

namespace src_pw {

template <typename FPTYPE, typename Device>
struct cal_vkb1_nl_op {
    void operator() (
        const Device* ctx,
        const int& nkb,
        const int& npwx,
        const int &npwk_max,
        const int& vkb_nc,
        const int& nbasis,
        const int& ik,
        const int& ipol,
        const std::complex<FPTYPE>& NEG_IMAG_UNIT,
        const std::complex<FPTYPE>* vkb,
        const FPTYPE* gcar,
        std::complex<FPTYPE>* vkb1);
};

template <typename FPTYPE, typename Device>
struct cal_force_nl_op {
    void operator()(
        const psi::DEVICE_CPU *ctx,
        const bool &multi_proj,
        const int &nbands_occ,
        const int &wg_nc,
        const int &ntype,
        const int &spin,
        const int &deeq_2,
        const int &deeq_3,
        const int &deeq_4,
        const int &forcenl_nc,
        const int &nbands,
        const int &ik,
        const int &nkb,
        const int *atom_nh,
        const int *atom_na,
        const FPTYPE &tpiba,
        const FPTYPE *d_wg,
        const FPTYPE *deeq,
        const std::complex<FPTYPE> *becp,
        const std::complex<FPTYPE> *dbecp,
        FPTYPE *force);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
template <typename FPTYPE>
struct cal_vkb1_nl_op<FPTYPE, psi::DEVICE_GPU> {
    void operator() (
        const psi::DEVICE_GPU *ctx,
        const int &nkb,
        const int &npwx,
        const int &npwk_max,
        const int &vkb_nc,
        const int &nbasis,
        const int &ik,
        const int &ipol,
        const std::complex<FPTYPE> &NEG_IMAG_UNIT,
        const std::complex<FPTYPE> *vkb,
        const FPTYPE *gcar,
        std::complex<FPTYPE> *vkb1);
};

template <typename FPTYPE>
struct cal_force_nl_op<FPTYPE, psi::DEVICE_GPU> {
    void operator() (
        const psi::DEVICE_GPU *ctx,
        const bool &multi_proj,
        const int &nbands_occ,
        const int &wg_nc,
        const int &ntype,
        const int &spin,
        const int &deeq_2,
        const int &deeq_3,
        const int &deeq_4,
        const int &forcenl_nc,
        const int &nbands,
        const int &ik,
        const int &nkb,
        const int *atom_nh,
        const int *atom_na,
        const FPTYPE &tpiba,
        const FPTYPE *d_wg,
        const FPTYPE *deeq,
        const std::complex<FPTYPE> *becp,
        const std::complex<FPTYPE> *dbecp,
        FPTYPE *force);
};

#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
}  // namespace src_pw
#endif //SRC_PW_FORCE_MULTI_DEVICE_H