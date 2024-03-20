#ifndef SRC_PW_WF_MULTI_DEVICE_H
#define SRC_PW_WF_MULTI_DEVICE_H

#include "module_psi/psi.h"
#include <complex>

namespace hamilt {

template <typename FPTYPE, typename Device>
struct cal_sk_op {
    /// @brief The prestep of the calculation of getvnl for multi device
    ///
    /// Input Parameters
    /// @param ctx - which device this function runs on
    /// @param ik - the current k point
    /// @param ntype - the number of atomic type
    /// @param nx - the first dimension of FFT box
    /// @param ny - the second dimension of FFT box
    /// @param nz - the third dimension of FFT box
    /// @param rho_nx - the first dimension of dense FFT box
    /// @param rho_ny - the second dimension of dense FFT box
    /// @param rho_nz - the third dimension of dense FFT box
    /// @param npw - the number of planewaves of current k point
    /// @param npwx - the number of planewaves of all k point
    /// @param fftny - wfc_basis->fftny
    /// @param eigts1_nc - GlobalC::sf.eigts1.nc
    /// @param eigts2_nc - GlobalC::sf.eigts2.nc
    /// @param eigts3_nc - GlobalC::sf.eigts3.nc
    /// @param atom_na - GlobalC::ucell.atoms[ii].na
    /// @param igl2isz - wfc_basis->igl2isz_k
    /// @param is2fftixy - wfc_basis->is2fftixy
    /// @param TWO_PI - ModuleBase::TWO_PI
    /// @param kvec_c - GlobalC::kv.kvec_c
    /// @param atom_tau - GlobalC::ucell.atoms[it].tau
    /// @param eigts1 - GlobalC::sf.eigts1
    /// @param eigts2 - GlobalC::sf.eigts2
    /// @param eigts3 - GlobalC::sf.eigts3
    ///
    /// Output Parameters
    /// @param sk - output results matrix with size GlobalC::ucell.nat * npw
    void operator()(const Device* ctx,
                    const int& ik,
                    const int& ntype,
                    const int& nx,
                    const int& ny,
                    const int& nz,
                    const int& rho_nx,
                    const int& rho_ny,
                    const int& rho_nz,
                    const int& npw,
                    const int& npwx,
                    const int& fftny,
                    const int& eigts1_nc,
                    const int& eigts2_nc,
                    const int& eigts3_nc,
                    const int* atom_na,
                    const int* igl2isz,
                    const int* is2fftixy,
                    const FPTYPE& TWO_PI,
                    const FPTYPE* kvec_c,
                    const FPTYPE* atom_tau,
                    std::complex<FPTYPE>* eigts1,
                    std::complex<FPTYPE>* eigts2,
                    std::complex<FPTYPE>* eigts3,
                    std::complex<FPTYPE>* sk);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
template <typename FPTYPE>
struct cal_sk_op<FPTYPE, psi::DEVICE_GPU> {
    void operator()(const psi::DEVICE_GPU* ctx,
                    const int& ik,
                    const int& ntype,
                    const int& nx,
                    const int& ny,
                    const int& nz,
                    const int& rho_nx,
                    const int& rho_ny,
                    const int& rho_nz,
                    const int& npw,
                    const int& npwx,
                    const int& fftny,
                    const int& eigts1_nc,
                    const int& eigts2_nc,
                    const int& eigts3_nc,
                    const int* atom_na,
                    const int* igl2isz,
                    const int* is2fftixy,
                    const FPTYPE& TWO_PI,
                    const FPTYPE* kvec_c,
                    const FPTYPE* atom_tau,
                    std::complex<FPTYPE>* eigts1,
                    std::complex<FPTYPE>* eigts2,
                    std::complex<FPTYPE>* eigts3,
                    std::complex<FPTYPE>* sk);
};
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
}  // namespace hamilt
#endif //SRC_PW_WF_MULTI_DEVICE_H