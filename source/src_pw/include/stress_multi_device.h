#ifndef SRC_PW_STRESS_MULTI_DEVICE_H
#define SRC_PW_STRESS_MULTI_DEVICE_H

#include "module_psi/psi.h"
#include <complex>

namespace src_pw {

    template <typename FPTYPE, typename Device>
    struct cal_dbecp_noevc_nl_op {
        void operator() (
                const Device *ctx,
                const int &ipol,
                const int &jpol,
                const int &nkb,
                const int &npw,
                const int &npwx,
                const int &ik,
                const FPTYPE &tpiba,
                const FPTYPE *gcar,
                const FPTYPE *kvec_c,
                std::complex<FPTYPE> *vkbi,
                std::complex<FPTYPE> *vkbj,
                std::complex<FPTYPE> *vkb,
                std::complex<FPTYPE> *vkb1,
                std::complex<FPTYPE> *vkb2,
                std::complex<FPTYPE> *dbecp_noevc);
    };

    template <typename FPTYPE, typename Device>
    struct cal_stress_nl_op {
        void operator()(
                const Device *ctx,
                const bool &multi_proj,
                const int &ipol,
                const int &jpol,
                const int &nkb,
                const int &nbands_occ,
                const int &ntype,
                const int &spin,
                const int &wg_nc,
                const int &ik,
                const int &deeq_2,
                const int &deeq_3,
                const int &deeq_4,
                const int *atom_nh,
                const int *atom_na,
                const FPTYPE *d_wg,
                const FPTYPE *deeq,
                const std::complex<FPTYPE> *becp,
                const std::complex<FPTYPE> *dbecp,
                FPTYPE *stress);
    };

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
template <typename FPTYPE>
struct cal_dbecp_noevc_nl_op<FPTYPE, psi::DEVICE_GPU> {
    void operator() (
        const psi::DEVICE_GPU *ctx,
        const int &ipol,
        const int &jpol,
        const int &nkb,
        const int &npw,
        const int &npwx,
        const int &ik,
        const FPTYPE &tpiba,
        const FPTYPE *gcar,
        const FPTYPE *kvec_c,
        std::complex<FPTYPE> *vkbi,
        std::complex<FPTYPE> *vkbj,
        std::complex<FPTYPE> *vkb,
        std::complex<FPTYPE> *vkb1,
        std::complex<FPTYPE> *vkb2,
        std::complex<FPTYPE> *dbecp_noevc);
};

template <typename FPTYPE>
struct cal_stress_nl_op<FPTYPE, psi::DEVICE_GPU> {
    void operator() (
        const psi::DEVICE_GPU *ctx,
        const bool &multi_proj,
        const int &ipol,
        const int &jpol,
        const int &nkb,
        const int &nbands_occ,
        const int &ntype,
        const int &spin,
        const int &wg_nc,
        const int &ik,
        const int &deeq_2,
        const int &deeq_3,
        const int &deeq_4,
        const int *atom_nh,
        const int *atom_na,
        const FPTYPE *d_wg,
        const FPTYPE *deeq,
        const std::complex<FPTYPE> *becp,
        const std::complex<FPTYPE> *dbecp,
        FPTYPE *stress);
};

#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
}  // namespace src_pw
#endif //SRC_PW_STRESS_MULTI_DEVICE_H