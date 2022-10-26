#ifndef MODULE_HAMILT_EKINETIC_H
#define MODULE_HAMILT_EKINETIC_H

#include "module_psi/psi.h"
#include <complex>

namespace hamilt {
template <typename FPTYPE, typename Device> 
struct ekinetic_pw_op {
  void operator() (
      const Device* dev,
      const int& nband,
      const int& npw,
      const int& max_npw,
      const FPTYPE& tpiba2,
      const FPTYPE* gk2_ik,
      std::complex<FPTYPE>* tmhpsi,
      const std::complex<FPTYPE>* tmpsi_in);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
// Partially specialize functor for psi::GpuDevice.
template <typename FPTYPE> 
struct ekinetic_pw_op<FPTYPE, psi::DEVICE_GPU> {
  void operator() (
      const psi::DEVICE_GPU* dev,
      const int& nband,
      const int& npw,
      const int& max_npw,
      const FPTYPE& tpiba2,
      const FPTYPE* gk2_ik,
      std::complex<FPTYPE>* tmhpsi,
      const std::complex<FPTYPE>* tmpsi_in);
};
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
} // namespace hamilt
#endif //MODULE_HAMILT_EKINETIC_H