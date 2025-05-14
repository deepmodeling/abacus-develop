#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_PW_HAMILT_PWDFT_KERNELS_VNL_OP_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_PW_HAMILT_PWDFT_KERNELS_VNL_OP_H

#include "module_psi/psi.h"

#include <complex>

namespace hamilt
{

template <typename FPTYPE, typename Device>
struct cal_vnl_op
{
    /// @brief Calculate the getvnl for multi-device
    ///
    /// Input Parameters
    /// @param ctx - which device this function runs on
    /// @param ntype - number of atomic type
    /// @param npw - number of planewaves of current k point
    /// @param npwx - number of planewaves of all k points
    /// @param tab_2 - the second dimension of the input table
    /// @param tab_3 - the third dimension of the input table
    /// @param atom_nh - ucell.atoms[ii].ncpp.nh
    /// @param atom_nb - ucell.atoms[it].ncpp.nbeta
    /// @param atom_na - ucell.atoms[ii].na
    /// @param DQ - PARAM.globalv.dq
    /// @param tpiba - ucell.tpiba
    /// @param NEG_IMAG_UNIT - ModuleBase::NEG_IMAG_UNIT
    /// @param gk - GlobalC::wf.get_1qvec_cartesian
    /// @param ylm - the result of ModuleBase::YlmReal::Ylm_Real
    /// @param indv - pseudopot_cell_vnl::indv
    /// @param nhtol - pseudopot_cell_vnl::nhtol
    /// @param nhtolm - pseudopot_cell_vnl::nhtolm
    /// @param tab - pseudopot_cell_vnl::tab
    /// @param vkb1 - intermedia matrix with size nkb * GlobalC::wf.npwx
    /// @param sk - results of cal_sk
    ///
    /// Output Parameters
    /// @param vkb_in - output results with size nkb * GlobalC::wf.npwx
    void operator()(const Device* ctx,
                    const int& ntype,
                    const int& npw,
                    const int& npwx,
                    const int& nhm,
                    const int& tab_2,
                    const int& tab_3,
                    const int* atom_na,
                    const int* atom_nb,
                    const int* atom_nh,
                    const FPTYPE& DQ,
                    const FPTYPE& tpiba,
                    const std::complex<FPTYPE>& NEG_IMAG_UNIT,
                    const FPTYPE* gk,
                    const FPTYPE* ylm,
                    const FPTYPE* indv,
                    const FPTYPE* nhtol,
                    const FPTYPE* nhtolm,
                    const FPTYPE* tab,
                    FPTYPE* vkb1,
                    const std::complex<FPTYPE>* sk,
                    std::complex<FPTYPE>* vkb_in);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
template <typename FPTYPE>
struct cal_vnl_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const int& ntype,
                    const int& npw,
                    const int& npwx,
                    const int& nhm,
                    const int& tab_2,
                    const int& tab_3,
                    const int* atom_na,
                    const int* atom_nb,
                    const int* atom_nh,
                    const FPTYPE& DQ,
                    const FPTYPE& tpiba,
                    const std::complex<FPTYPE>& NEG_IMAG_UNIT,
                    const FPTYPE* gk,
                    const FPTYPE* ylm,
                    const FPTYPE* indv,
                    const FPTYPE* nhtol,
                    const FPTYPE* nhtolm,
                    const FPTYPE* tab,
                    FPTYPE* vkb1,
                    const std::complex<FPTYPE>* sk,
                    std::complex<FPTYPE>* vkb_in);
};
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
} // namespace hamilt
#endif // W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_PW_HAMILT_PWDFT_KERNELS_VNL_OP_H