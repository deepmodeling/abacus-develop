// TODO: This is a temperary location for these functions.
// And will be moved to a global module(module base) later.
#ifndef MODULE_ELECSTATE_ELECSTATE_MULTI_DEVICE_H
#define MODULE_ELECSTATE_ELECSTATE_MULTI_DEVICE_H
#include <complex>
#include "module_psi/psi.h"

namespace elecstate{

template <typename FPTYPE, typename Device> 
struct elecstate_pw_op {
  void operator() (
      const Device* ctx,
      const int& spin,
      const int& nrxx,
      const FPTYPE& weight,
      FPTYPE** rho,
      const std::complex<FPTYPE>* wfcr);
  
  void operator() (
      const Device* ctx,
      const bool& DOMAG,
      const bool& DOMAG_Z,
      const int& nrxx,
      const FPTYPE& weight,
      FPTYPE** rho,
      const std::complex<FPTYPE>* wfcr,
      const std::complex<FPTYPE>* wfcr_another_spin);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
template <typename FPTYPE> 
struct elecstate_pw_op<FPTYPE, psi::DEVICE_GPU> {
  void operator()(
    const psi::DEVICE_GPU* ctx,
    const int& spin,
    const int& nrxx,
    const FPTYPE& w1,
    FPTYPE** rho,
    const std::complex<FPTYPE>* wfcr);

  void operator()(
    const psi::DEVICE_GPU* ctx,
    const bool& DOMAG,
    const bool& DOMAG_Z,
    const int& nrxx,
    const FPTYPE& w1,
    FPTYPE** rho,
    const std::complex<FPTYPE>* wfcr,
    const std::complex<FPTYPE>* wfcr_another_spin);
};
#endif
} // namespace elecstate

#endif // MODULE_ELECSTATE_ELECSTATE_MULTI_DEVICE_H