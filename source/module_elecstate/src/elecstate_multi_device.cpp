#include "module_elecstate/include/elecstate_multi_device.h"
#include <iostream>

namespace elecstate{

template <typename FPTYPE> 
struct elecstate_pw_op<FPTYPE, psi::DEVICE_CPU> {
  void operator()(
    const psi::DEVICE_CPU * ctx,
    const int& spin,
    const int& nrxx,
    const FPTYPE& w1,
    FPTYPE** rho,
    const std::complex<FPTYPE>* wfcr) 
    {
      // for (int ir = 0; ir < nrxx; ir++)
      // {
      //   rho[spin][ir] += weight * norm(wfcr[ir]);
      // }
      for (int ir = 0; ir < nrxx; ir++)
      {
        rho[spin][ir] += w1 * norm(wfcr[ir]);
      }
    }

   void operator()(
    const psi::DEVICE_CPU * ctx,
    const bool& DOMAG,
    const bool& DOMAG_Z,
    const int& nrxx,
    const FPTYPE& w1,
    FPTYPE** rho,
    const std::complex<FPTYPE>* wfcr,
    const std::complex<FPTYPE>* wfcr_another_spin) 
    {
      for (int ir = 0; ir < nrxx; ir++) {
          rho[0][ir] += w1 * (norm(wfcr[ir]) + norm(wfcr_another_spin[ir]));
      }
      // In this case, calculate the three components of the magnetization
      if (DOMAG)
      {
          for (int ir = 0; ir < nrxx; ir++) {
              rho[1][ir] += w1 * 2.0
                                          * (wfcr[ir].real() * wfcr_another_spin[ir].real()
                                             + wfcr[ir].imag() * wfcr_another_spin[ir].imag());
              rho[2][ir] += w1 * 2.0
                                          * (wfcr[ir].real() * wfcr_another_spin[ir].imag()
                                             - wfcr_another_spin[ir].real() * wfcr[ir].imag());
              rho[3][ir] += w1 * (norm(wfcr[ir]) - norm(wfcr_another_spin[ir]));
          }
      }
      else if (DOMAG_Z)
      {
          for (int ir = 0; ir < nrxx; ir++)
          {
              rho[1][ir] = 0;
              rho[2][ir] = 0;
              rho[3][ir] += w1 * (norm(wfcr[ir]) - norm(wfcr_another_spin[ir]));
          }
      }
      else {
          for (int is = 1; is < 4; is++)
          {
              for (int ir = 0; ir < nrxx; ir++)
                  rho[is][ir] = 0;
          }
      }
    } 
};


template struct elecstate_pw_op<double, psi::DEVICE_CPU>;
}