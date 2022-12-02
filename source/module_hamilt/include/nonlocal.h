#ifndef MODULE_HAMILT_NONLOCAL_H
#define MODULE_HAMILT_NONLOCAL_H

#include "module_psi/psi.h"
#include <complex>

namespace hamilt {
template <typename FPTYPE, typename Device> 
struct nonlocal_pw_op {
  void operator() (
      const Device* dev,
      const int& l1,
      const int& l2,
      const int& l3,
      int& sum,
      int& iat,
      const int& spin,
      const int& nkb,
      const int& deeq_x,
      const int& deeq_y,
      const int& deeq_z,
      const FPTYPE* deeq,
      std::complex<FPTYPE>* ps,
      const std::complex<FPTYPE>* becp);
  
  void operator() (
      const Device* dev,
      const int& l1,
      const int& l2,
      const int& l3,
      int& sum,
      int& iat,
      const int& nkb,
      const int& deeq_x,
      const int& deeq_y,
      const int& deeq_z,
      const std::complex<FPTYPE>* deeq_nc,
      std::complex<FPTYPE>* ps,
      const std::complex<FPTYPE>* becp);
};
                      
#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
// Partially specialize functor for psi::GpuDevice.
template <typename FPTYPE> 
struct nonlocal_pw_op<FPTYPE, psi::DEVICE_GPU> {
  void operator() (
      const psi::DEVICE_GPU* dev,
      const int& l1,
      const int& l2,
      const int& l3,
      int& sum,
      int& iat,
      const int& spin,
      const int& nkb,
      const int& deeq_x,
      const int& deeq_y,
      const int& deeq_z,
      const FPTYPE* deeq,
      std::complex<FPTYPE>* ps,
      const std::complex<FPTYPE>* becp);
  
  void operator() (
      const psi::DEVICE_GPU* dev,
      const int& l1,
      const int& l2,
      const int& l3,
      int& sum,
      int& iat,
      const int& nkb,
      const int& deeq_x,
      const int& deeq_y,
      const int& deeq_z,
      const std::complex<FPTYPE>* deeq_nc,
      std::complex<FPTYPE>* ps,
      const std::complex<FPTYPE>* becp);
};
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
} // namespace hamilt
#endif //MODULE_HAMILT_NONLOCAL_H