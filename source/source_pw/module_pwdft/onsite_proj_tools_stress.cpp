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
double Onsite_Proj_tools<FPTYPE, Device>::cal_stress_dftu(int ik,
                                                          int npm,
                                                          const int* orb_corr,
                                                          const std::complex<FPTYPE>* pot_onsite,
                                                          const int size_pot_onsite,
                                                          const FPTYPE* h_wg)
{
    double stress_out = 0.0;
    const int npol = this->ucell_->get_npol();
    
    int* orb_corr_tmp = nullptr;
    std::complex<FPTYPE>* pot_onsite_tmp = nullptr;
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
	// orb_corr_tmp
        resmem_int_op()(orb_corr_tmp, this->ucell_->ntype);
        syncmem_int_h2d_op()(orb_corr_tmp, orb_corr, this->ucell_->ntype);

	// pot_onsite_tmp
        resmem_complex_op()(pot_onsite_tmp, size_pot_onsite);
        syncmem_complex_h2d_op()(pot_onsite_tmp, pot_onsite, size_pot_onsite);

	// transfer data from from host to device
        syncmem_var_h2d_op()(d_wg, h_wg, this->nbands * (ik+1));
        
        // Allocate device memory for stress
        FPTYPE* stress_device = nullptr;
        resmem_var_op()(stress_device, 1);
        setmem_var_op()(stress_device, 0, 1);
        
        cal_stress_nl_op()(this->ctx,
                           nkb,
                           npm,
                           this->ntype,
                           this->nbands,
                           ik,
                           npol,
                           atom_nh,
                           atom_na,
                           d_wg,
                           pot_onsite_tmp,
                           orb_corr_tmp,
                           becp,
                           dbecp,
                           stress_device);
        
        // Transfer stress from device to host
        syncmem_var_d2h_op()(&stress_out, stress_device, 1);
        delmem_var_op()(stress_device);
        delmem_complex_op()(pot_onsite_tmp);
        delmem_int_op()(orb_corr_tmp);
    }
    else
#endif
    {
        orb_corr_tmp = const_cast<int*>(orb_corr);
        pot_onsite_tmp = const_cast<std::complex<FPTYPE>*>(pot_onsite);
        d_wg = const_cast<FPTYPE*>(h_wg);
        
        cal_stress_nl_op()(this->ctx,
                           nkb,
                           npm,
                           this->ntype,
                           this->nbands,
                           ik,
                           npol,
                           atom_nh,
                           atom_na,
                           d_wg,
                           pot_onsite_tmp,
                           orb_corr_tmp,
                           becp,
                           dbecp,
                           &stress_out);
//	std::cout << "DFT+U (CPU) stress_out = " << stress_out << std::endl;
    }
    
    return stress_out;
}

template <typename FPTYPE, typename Device>
double Onsite_Proj_tools<FPTYPE, Device>::cal_stress_dspin(int ik,
                                                           int npm,
					          	   const ModuleBase::Vector3<double>* lambda,
                                                           const FPTYPE* h_wg)
{
    double stress_out = 0.0;
    const int npol = this->ucell_->get_npol();
    
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
        
        // Allocate device memory for stress
        FPTYPE* stress_device = nullptr;
        resmem_var_op()(stress_device, 1);
        setmem_var_op()(stress_device, 0, 1);
        
        const int force_nc = 3;
        cal_stress_nl_op()(this->ctx,
                           nkb,
                           npm,
                           this->ntype,
                           this->nbands,
                           ik,
                           npol,
                           atom_nh,
                           atom_na,
                           d_wg,
                           lambda_tmp,
                           this->kv_ ? this->kv_->isk.data() : nullptr,
                           becp,
                           dbecp,
                           stress_device);
        
        // Transfer stress from device to host
        syncmem_var_d2h_op()(&stress_out, stress_device, 1);
        delmem_var_op()(stress_device);
        delmem_var_op()(lambda_tmp);
    }
    else
#endif
    {
        lambda_tmp = lambda_array.data();
        d_wg = const_cast<FPTYPE*>(h_wg);
        
        const int force_nc = 3;
        cal_stress_nl_op()(this->ctx,
                           nkb,
                           npm,
                           this->ntype,
                           this->nbands,
                           ik,
                           npol,
                           atom_nh,
                           atom_na,
                           d_wg,
                           lambda_tmp,
                           this->kv_ ? this->kv_->isk.data() : nullptr,
                           becp,
                           dbecp,
                           &stress_out);
    }
    
    return stress_out;
}

// template instantiation
template double Onsite_Proj_tools<double, base_device::DEVICE_CPU>::cal_stress_dftu(
    int, int, const int*, const std::complex<double>*, const int, const double*);
template double Onsite_Proj_tools<double, base_device::DEVICE_CPU>::cal_stress_dspin(
    int, int, const ModuleBase::Vector3<double>*, const double*);
#if ((defined __CUDA) || (defined __ROCM))
template double Onsite_Proj_tools<double, base_device::DEVICE_GPU>::cal_stress_dftu(
    int, int, const int*, const std::complex<double>*, const int, const double*);
template double Onsite_Proj_tools<double, base_device::DEVICE_GPU>::cal_stress_dspin(
    int, int, const ModuleBase::Vector3<double>*, const double*);
#endif

} // namespace hamilt
