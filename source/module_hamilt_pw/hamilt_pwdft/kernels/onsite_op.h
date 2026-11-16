#ifndef MODULE_HAMILT_OPERATOR_KERNELS_ONSITE_H
#define MODULE_HAMILT_OPERATOR_KERNELS_ONSITE_H

#include "module_psi/psi.h"
#include <complex>

namespace hamilt {
#if defined(__CUDA)
// Applies the scalar or Pauli constraint field to column-major projector
// coefficients with shape (nkb * npol, bands). Lambda has shape (nkb, 3)
// in projector-major order. The output is fully overwritten.
void spin_constraint_projectors_gpu(int bands, int npol, int nkb, int sign,
                                    const double* lambda, const std::complex<double>* becp,
                                    std::complex<double>* ps);

// Forms unweighted Cartesian moments with shape (bands, atoms, 3) from
// column-major projections with shape (nkb * npol, bands). Offsets contains
// atoms + 1 cumulative projector counts. Each atom's projector order is retained.
void projector_band_moments_gpu(int bands, int npol, int atoms, int nkb, const int* offsets,
                                const std::complex<double>* projections, double* moments);

// Weights moments with shape (points, bands, atoms, 3) by occupations with
// shape (points, bands), retaining point-then-band accumulation order.
// Spin labels select +z for zero and -z otherwise when npol is one.
void sum_projector_moments_gpu(int points, int bands, int atoms, int npol,
                              const double* moments, const double* weights,
                              const int* spins, double* output);
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
