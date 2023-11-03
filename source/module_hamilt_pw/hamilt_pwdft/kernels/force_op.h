#ifndef SRC_PW_FORCE_MULTI_DEVICE_H
#define SRC_PW_FORCE_MULTI_DEVICE_H

#include "module_psi/psi.h"
#include <complex>

namespace hamilt {

template <typename FPTYPE, typename Device>
struct cal_vkb1_nl_op {
    /// @brief The prestep to calculate the final forces
    ///
    /// Input Parameters
    /// @param ctx - which device this function runs on
    /// @param nkb - number of k point
    /// @param npwx - number of planewaves
    /// @param npwk_max - number of planewaves
    /// @param vkb_nc - the second dimension of vkb matrix
    /// @param nbasis - number of planewaves of current k point
    /// @param ik - current k point
    /// @param ipol - 0,1,2
    /// @param NEG_IMAG_UNIT - ModuleBase::NEG_IMAG_UNIT
    /// @param vkb - result of getvnl
    /// @param gcar - GlobalC::wfcpw->gcar
    ///
    /// Output Parameters
    /// @param vkb1 - output vkb matrix
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
    /// @brief Calculate the final forces for multi-device
    ///
    /// Input Parameters
    /// @param ctx - which device this function runs on
    /// @param nondiagonal - control flag
    /// @param nbands_occ - number of occupied bands
    /// @param wg_nc - the second dimension of matrix wg
    /// @param ntype - total atomic type
    /// @param spin - current spin
    /// @param deeq_2 - the second dimension of deeq
    /// @param deeq_3 - the third dimension of deeq
    /// @param deeq_4 - the forth dimension of deeq
    /// @param forcenl_nc - the second dimension of matrix forcenl
    /// @param nbands - NBANDS
    /// @param ik - current k point
    /// @param nkb - number of k point
    /// @param atom_nh - GlobalC::ucell.atoms[ii].ncpp.nh
    /// @param atom_na - GlobalC::ucell.atoms[ii].na
    /// @param tpiba - GlobalC::ucell.tpiba
    /// @param d_wg - input parameter wg
    /// @param d_ekb - input parameter ekb
    /// @param qq_nt - GlobalC::ppcell.qq_nt
    /// @param deeq - GlobalC::ppcell.deeq
    /// @param becp - intermediate matrix with GlobalV::NBANDS * nkb
    /// @param dbecp - intermediate matrix with 3 * GlobalV::NBANDS * nkb
    ///
    /// Output Parameters
    /// @param force - output forces
    void operator()(const psi::DEVICE_CPU* ctx,
                    const bool& nondiagonal,
                    const int& nbands_occ,
                    const int& wg_nc,
                    const int& ntype,
                    const int& spin,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& ik,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const FPTYPE* deeq,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* force);
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
    void operator()(const psi::DEVICE_GPU* ctx,
                    const bool& nondiagonal,
                    const int& nbands_occ,
                    const int& wg_nc,
                    const int& ntype,
                    const int& spin,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& ik,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const FPTYPE* deeq,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* force);
};

#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
}  // namespace hamilt
#endif //SRC_PW_FORCE_MULTI_DEVICE_H