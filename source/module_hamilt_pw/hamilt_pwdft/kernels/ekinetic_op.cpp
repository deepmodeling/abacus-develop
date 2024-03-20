#include "module_hamilt_pw/hamilt_pwdft/kernels/ekinetic_op.h"

namespace hamilt {

template <typename FPTYPE> 
struct ekinetic_pw_op<FPTYPE, psi::DEVICE_CPU> {
  void operator() (
      const psi::DEVICE_CPU* /*dev*/,
      const int& nband,
      const int& npw,
      const int& max_npw,
      const FPTYPE& tpiba2,
      const FPTYPE* gk2_ik,
      std::complex<FPTYPE>* tmhpsi,
      const std::complex<FPTYPE>* tmpsi_in)
  {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 4096/sizeof(FPTYPE))
#endif
    for (int ib = 0; ib < nband; ++ib) {
      for (int ig = 0; ig < npw; ++ig) {
        tmhpsi[ib * max_npw + ig] += gk2_ik[ig] * tpiba2 * tmpsi_in[ib * max_npw + ig];
      }
    }
  }
};

template struct ekinetic_pw_op<float, psi::DEVICE_CPU>;
template struct ekinetic_pw_op<double, psi::DEVICE_CPU>;

}  // namespace hamilt

