#ifndef SRC_PW_STRESS_MULTI_DEVICE_H
#define SRC_PW_STRESS_MULTI_DEVICE_H

#include "module_psi/psi.h"
#include <complex>

namespace hamilt {

    template <typename FPTYPE, typename Device>
    struct cal_dbecp_noevc_nl_op {
        /// @brief The prestep to calculate the final stresses
        ///
        /// Input Parameters
        /// @param ctx - which device this function runs on
        /// @param ipol - loop of 0, 1, 2
        /// @param jpol - loop of 0, 1, 2
        /// @param nkb - number of k points
        /// @param npw - number of planewaves
        /// @param npwx - max number of planewaves
        /// @param ik - the current k point
        /// @param tpiba - GlobalC::ucell.tpiba
        /// @param gcar - GlobalC::wfcpw->gcar
        /// @param kvec_c - GlobalC::wfcpw->kvec_c
        /// @param vkbi - _vkb0[ipol]
        /// @param vkbj - _vkb0[jpol]
        /// @param vkb - result of getvnl
        /// @param vkb1 - intermeidate matrix with size nkb * GlobalC::wf.npwx
        /// @param vkb2 - intermeidate matrix with size nkb * GlobalC::wf.npwx
        ///
        /// Output Parameters
        /// @param dbecp_noevc - result matrix with size nkb * GlobalC::wf.npwx
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
        /// @brief Calculate the final stresses for multi-device
        ///
        /// Input Parameters
        /// @param ctx - which device this function runs on
        /// @param nondiagonal - control flag
        /// @param ipol - loop of 0, 1, 2
        /// @param jpol - loop of 0, 1, 2
        /// @param nkb - number of k point
        /// @param nbands_occ - number of occupied bands
        /// @param ntype - total atomic type
        /// @param spin - current spin
        /// @param wg_nc - the second dimension of matrix wg
        /// @param ik - current k point
        /// @param deeq_2 - the second dimension of deeq
        /// @param deeq_3 - the third dimension of deeq
        /// @param deeq_4 - the forth dimension of deeq
        /// @param atom_nh - GlobalC::ucell.atoms[ii].ncpp.nh
        /// @param atom_na - GlobalC::ucell.atoms[ii].na
        /// @param d_wg - input parameter wg
        /// @param d_ekb - input parameter ekb
        /// @param qq_nt - GlobalC::ppcell.qq_nt
        /// @param deeq - GlobalC::ppcell.deeq
        /// @param becp - intermediate matrix with GlobalV::NBANDS * nkb
        /// @param dbecp - intermediate matrix with 3 * GlobalV::NBANDS * nkb
        ///
        /// Output Parameters
        /// @param stress - output stresses
        void operator()(const Device* ctx,
                        const bool& nondiagonal,
                        const int& ipol,
                        const int& jpol,
                        const int& nkb,
                        const int& nbands_occ,
                        const int& ntype,
                        const int& spin,
                        const int& wg_nc,
                        const int& ik,
                        const int& deeq_2,
                        const int& deeq_3,
                        const int& deeq_4,
                        const int* atom_nh,
                        const int* atom_na,
                        const FPTYPE* d_wg,
                        const FPTYPE* d_ekb,
                        const FPTYPE* qq_nt,
                        const FPTYPE* deeq,
                        const std::complex<FPTYPE>* becp,
                        const std::complex<FPTYPE>* dbecp,
                        FPTYPE* stress);
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
    void operator()(const psi::DEVICE_GPU* ctx,
                    const bool& nondiagonal,
                    const int& ipol,
                    const int& jpol,
                    const int& nkb,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& spin,
                    const int& wg_nc,
                    const int& ik,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE* d_wg,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const FPTYPE* deeq,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* stress);
};

#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
}  // namespace hamilt
#endif //SRC_PW_STRESS_MULTI_DEVICE_H