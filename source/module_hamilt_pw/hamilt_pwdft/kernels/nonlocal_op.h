#ifndef MODULE_HAMILT_NONLOCAL_H
#define MODULE_HAMILT_NONLOCAL_H

#include "module_psi/psi.h"
#include <complex>

namespace hamilt {
template <typename FPTYPE, typename Device> 
struct nonlocal_pw_op {
  /// @brief Compute the nonlocal potential of hPsi
  ///
  /// Input Parameters
  /// \param dev : the type of computing device
  /// \param l1 : ucell->atoms[it].na
  /// \param l2 : nbands
  /// \param l3 : ucell->atoms[it].ncpp.nh
  /// \param sum : intermediate value
  /// \param iat : intermediate value
  /// \param spin : current spin
  /// \param nkb : ppcell->nkb, number of kpoints
  /// \param deeq_x : second dimension of deeq
  /// \param deeq_y : third dimension of deeq
  /// \param deeq_z : forth dimension of deeq
  /// \param deeq : ppcell->deeq
  /// \param becp : intermediate array
  ///
  /// Output Parameters
  /// \param ps : output array
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

  /// @brief Compute the nonlocal potential of hPsi, with NSPIN > 2
  ///
  /// Input Parameters
  /// \param dev : the type of computing device
  /// \param l1 : ucell->atoms[it].na
  /// \param l2 : nbands
  /// \param l3 : ucell->atoms[it].ncpp.nh
  /// \param sum : intermediate value
  /// \param iat : intermediate value
  /// \param nkb : ppcell->nkb, number of kpoints
  /// \param deeq_x : second dimension of deeq
  /// \param deeq_y : third dimension of deeq
  /// \param deeq_z : forth dimension of deeq
  /// \param deeq_nc : ppcell->deeq_nc
  /// \param becp : intermediate array
  ///
  /// Output Parameters
  /// \param ps : output array
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