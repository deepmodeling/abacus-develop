#include "onsite_proj_tools.h"

#include "source_base/math_polyint.h"
#include "source_base/math_ylmreal.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_pw/module_pwdft/kernels/force_op.h"
#include "source_io/module_parameter/parameter.h"
#include "nonlocal_maths.hpp"

#include <numeric>

namespace hamilt
{

template <typename FPTYPE, typename Device>
void Onsite_Proj_tools<FPTYPE, Device>::cal_force_dftu(int ik,
                                                       int npm,
                                                       FPTYPE* force,
                                                       const int* orbital_corr,
                                                       const std::complex<FPTYPE>* pot_onsite,
                                                       const int size_pot_onsite,
                                                       const FPTYPE* h_wg)
{
    int* orbital_corr_tmp = nullptr;
    std::complex<FPTYPE>* pot_onsite_tmp = nullptr;
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        resmem_int_op()(orbital_corr_tmp, this->ucell_->ntype);
        syncmem_int_h2d_op()(orbital_corr_tmp, orbital_corr, this->ucell_->ntype);
        resmem_complex_op()(pot_onsite_tmp, size_pot_onsite);
        syncmem_complex_h2d_op()(pot_onsite_tmp, pot_onsite, size_pot_onsite);
        syncmem_var_h2d_op()(d_wg, h_wg, this->nbands * (ik+1));
    }
    else
#endif
    {
        orbital_corr_tmp = const_cast<int*>(orbital_corr);
        pot_onsite_tmp = const_cast<std::complex<FPTYPE>*>(pot_onsite);
        d_wg = const_cast<FPTYPE*>(h_wg);
    }
    const int force_nc = 3;
    const int npol = this->ucell_->get_npol();
    cal_force_nl_op<FPTYPE, Device>()(this->ctx,
                                      npm,
                                      this->nbands,
                                      this->ntype,
                                      force_nc,
                                      this->nbands,
                                      ik,
                                      nkb,
                                      npol,
                                      atom_nh,
                                      atom_na,
                                      this->ucell_->tpiba,
                                      d_wg,
                                      pot_onsite_tmp,
                                      orbital_corr_tmp,
                                      becp,
                                      dbecp,
                                      force);
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        delmem_complex_op()(pot_onsite_tmp);
        delmem_int_op()(orbital_corr_tmp);
    }
#endif
}

template <typename FPTYPE, typename Device>
void Onsite_Proj_tools<FPTYPE, Device>::cal_force_dspin(int ik,
                                                        int npm,
                                                        FPTYPE* force,
                                                        const ModuleBase::Vector3<double>* lambda,
                                                        const FPTYPE* h_wg)
{
    std::vector<FPTYPE> lambda_array(this->ucell_->nat * 3);
    for (int iat = 0; iat < this->ucell_->nat; iat++)
    {
        lambda_array[iat * 3] = lambda[iat].x;
        lambda_array[iat * 3 + 1] = lambda[iat].y;
        lambda_array[iat * 3 + 2] = lambda[iat].z;
    }
    FPTYPE* lambda_tmp = nullptr;
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        resmem_var_op()(lambda_tmp, this->ucell_->nat * 3);
        syncmem_var_h2d_op()(lambda_tmp, lambda_array.data(), this->ucell_->nat * 3);
        syncmem_var_h2d_op()(d_wg, h_wg, this->nbands * (ik+1));
    }
    else
#endif
    {
        lambda_tmp = lambda_array.data();
        d_wg = const_cast<FPTYPE*>(h_wg);
    }
    const int force_nc = 3;
    const int npol = this->ucell_->get_npol();
    cal_force_nl_op<FPTYPE, Device>()(this->ctx,
                                       npm,
                                       this->nbands,
                                       this->ntype,
                                       force_nc,
                                       this->nbands,
                                       ik,
                                       nkb,
                                       npol,
                                       atom_nh,
                                       atom_na,
                                       this->ucell_->tpiba,
                                       d_wg,
                                       lambda_tmp,
                                       this->kv_ ? this->kv_->isk.data() : nullptr,
                                       becp,
                                       dbecp,
                                       force);

#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        delmem_var_op()(lambda_tmp);
    }
#endif
}

// template instantiation
template void Onsite_Proj_tools<double, base_device::DEVICE_CPU>::cal_force_dftu(
    int, int, double*, const int*, const std::complex<double>*, const int, const double*);
template void Onsite_Proj_tools<double, base_device::DEVICE_CPU>::cal_force_dspin(
    int, int, double*, const ModuleBase::Vector3<double>*, const double*);
#if ((defined __CUDA) || (defined __ROCM))
template void Onsite_Proj_tools<double, base_device::DEVICE_GPU>::cal_force_dftu(
    int, int, double*, const int*, const std::complex<double>*, const int, const double*);
template void Onsite_Proj_tools<double, base_device::DEVICE_GPU>::cal_force_dspin(
    int, int, double*, const ModuleBase::Vector3<double>*, const double*);
#endif

} // namespace hamilt
