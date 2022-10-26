#include "module_hamilt/include/ekinetic.h"

#include <iomanip>
#include <iostream>

using namespace hamilt; 

template <typename FPTYPE> 
struct hamilt::ekinetic_pw_op<FPTYPE, psi::DEVICE_CPU> {
  void operator() (
      const psi::DEVICE_CPU* dev,
      const int& nband,
      const int& npw,
      const int& max_npw,
      const FPTYPE& tpiba2,
      const FPTYPE* gk2_ik,
      std::complex<FPTYPE>* tmhpsi,
      const std::complex<FPTYPE>* tmpsi_in)
  {
    for (int ib = 0; ib < nband; ++ib) {
      for (int ig = 0; ig < npw; ++ig) {
        tmhpsi[ig] += gk2_ik[ig] * tpiba2 * tmpsi_in[ig];
      }
      tmhpsi += max_npw;
      tmpsi_in += max_npw;
    }
  }
};

namespace hamilt{
template struct ekinetic_pw_op<double, psi::DEVICE_CPU>;
}

