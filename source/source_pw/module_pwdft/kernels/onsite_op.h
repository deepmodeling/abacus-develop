#ifndef MODULE_HAMILT_OPERATOR_KERNELS_ONSITE_H
#define MODULE_HAMILT_OPERATOR_KERNELS_ONSITE_H

#include "source_base/module_device/types.h"
#include <complex>

namespace hamilt {
#if defined(__CUDA)
// Applies the scalar or Pauli constraint field to column-major projector
// coefficients with shape (nkb * npol, bands). Lambda has shape (nkb, 3)
// in projector-major order. The output is fully overwritten.
void spin_constraint_projectors_gpu(int bands, int npol, int nkb, int sign,
                                    const double* lambda, const std::complex<double>* becp,
                                    std::complex<double>* ps);

// Accumulate one k point's weighted atomic magnetic moments into mi_xyz.
// projections has shape (nkb * npol, bands), offsets has atoms + 1 cumulative
// projector counts, and weights has nbands entries for this k point.
void accumulate_projector_moments_gpu(int bands, int npol, int atoms, int nkb, int sign,
                                      const int* offsets, const double* weights,
                                      const std::complex<double>* projections, double* mi_xyz);
#endif

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
      const int* pot_onsite_begin_iat,
      const int& tnp,
      const std::complex<FPTYPE>* pot_onsite,
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
      const int* pot_onsite_begin_iat,
      const int& tnp,
      const std::complex<FPTYPE>* pot_onsite,
      std::complex<FPTYPE>* ps,
      const std::complex<FPTYPE>* becp);
};
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
} // namespace hamilt
#endif //MODULE_HAMILT_OPERATOR_KERNELS_ONSITE_H