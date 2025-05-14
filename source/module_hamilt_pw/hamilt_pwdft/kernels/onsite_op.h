#ifndef MODULE_HAMILT_OPERATOR_KERNELS_ONSITE_H
#define MODULE_HAMILT_OPERATOR_KERNELS_ONSITE_H

#include "module_psi/psi.h"
#include <complex>

namespace hamilt {
template <typename FPTYPE, typename Device> 
struct onsite_ps_op {
  void operator() (
      const Device* dev,
      const int& npm,
      const int npol,
      const int* ip_iat,
      const int& tnp,
      const std::complex<FPTYPE>* lambda_coeff,
      std::complex<FPTYPE>* ps,
      const std::complex<FPTYPE>* becp);

  void operator() (
      const Device* dev,
      const int& npm,
      const int npol,
      const int* orb_l_iat,
      const int* ip_iat,
      const int* ip_m,
      const int* vu_begin_iat,
      const int& tnp,
      const std::complex<FPTYPE>* vu,
      std::complex<FPTYPE>* ps,
      const std::complex<FPTYPE>* becp);
};
                      
#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
// Partially specialize functor for base_device::GpuDevice.
template <typename FPTYPE>
struct onsite_ps_op<FPTYPE, base_device::DEVICE_GPU> {
  void operator() (
      const base_device::DEVICE_GPU* dev,
      const int& npm,
      const int npol,
      const int* ip_iat,
      const int& tnp,
      const std::complex<FPTYPE>* lambda_coeff,
      std::complex<FPTYPE>* ps,
      const std::complex<FPTYPE>* becp);

  void operator() (
      const base_device::DEVICE_GPU* dev,
      const int& npm,
      const int npol,
      const int* orb_l_iat,
      const int* ip_iat,
      const int* ip_m,
      const int* vu_begin_iat,
      const int& tnp,
      const std::complex<FPTYPE>* vu,
      std::complex<FPTYPE>* ps,
      const std::complex<FPTYPE>* becp);
};
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
} // namespace hamilt
#endif //MODULE_HAMILT_OPERATOR_KERNELS_ONSITE_H