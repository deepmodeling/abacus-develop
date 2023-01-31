#ifndef MODULE_HAMILT_EKINETIC_H
#define MODULE_HAMILT_EKINETIC_H

#include "module_psi/psi.h"
#include <complex>

namespace hamilt {
template <typename FPTYPE, typename Device> 
struct ekinetic_pw_op {
  /// @brief Compute the ekinetic potential of hPsi
  ///
  /// Input Parameters
  /// \param dev : the type of computing device
  /// \param nband : nbands
  /// \param npw : number of planewaves of current k point
  /// \param max_npw : max number of planewaves of all k points
  /// \param tpiba2 : GlobalC::ucell.tpiba2
  /// \param spin : current spin
  /// \param gk2_ik : GlobalC::wfcpw->gk2
  /// \param tmpsi_in : intermediate array
  ///
  /// Output Parameters
  /// \param tmhpsi : output array
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