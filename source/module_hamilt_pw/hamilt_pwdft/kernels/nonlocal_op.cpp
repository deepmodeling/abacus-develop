#include "module_hamilt_pw/hamilt_pwdft/kernels/nonlocal_op.h"

namespace hamilt {

template <typename FPTYPE> 
struct nonlocal_pw_op<FPTYPE, psi::DEVICE_CPU> {
  void operator() (
      const psi::DEVICE_CPU* /*dev*/,
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
      const std::complex<FPTYPE>* becp)
  {
#ifdef _OPENMP
#pragma omp parallel for collapse(3)
#endif
    for (int ii = 0; ii < l1; ii++) {
      // each atom has nproj, means this is with structure factor;
      // each projector (each atom) must multiply coefficient
      // with all the other projectors.
      for (int jj = 0; jj < l2; ++jj) 
        for (int kk = 0; kk < l3; kk++) 
          for (int xx = 0; xx < l3; xx++) 
            ps[(sum + ii * l3 + kk) * l2 + jj]
              += deeq[((spin * deeq_x + iat + ii) * deeq_y + xx) * deeq_z + kk] 
              *  becp[jj * nkb + sum + ii * l3 + xx];
    }
    sum += l1 * l3;
    iat += l1;
  }

  void operator() (
      const psi::DEVICE_CPU* dev,
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
      const std::complex<FPTYPE>* becp)
  {
#ifdef _OPENMP
#pragma omp parallel for collapse(3)
#endif
    for (int ii = 0; ii < l1; ii++)
      // each atom has nproj, means this is with structure factor;
      // each projector (each atom) must multiply coefficient
      // with all the other projectors.
      for (int jj = 0; jj < l2; jj+=2)
        for (int kk = 0; kk < l3; kk++)
          for (int xx = 0; xx < l3; xx++)
          {
              int psind = (sum + ii * l3 + kk) * l2 + jj;
              int becpind = jj * nkb + sum + ii * l3 + xx;
              auto &becp1 = becp[becpind];
              auto &becp2 = becp[becpind + nkb];
              ps[psind] += deeq_nc[((0 * deeq_x + iat + ii) * deeq_y + kk) * deeq_z + xx] * becp1
                           + deeq_nc[((1 * deeq_x + iat + ii) * deeq_y + kk) * deeq_z + xx] * becp2;
              ps[psind + 1] += deeq_nc[((2 * deeq_x + iat + ii) * deeq_y + kk) * deeq_z + xx] * becp1
                               + deeq_nc[((3 * deeq_x + iat + ii) * deeq_y + kk) * deeq_z + xx] * becp2;
          } // end jj
    iat += l1;
    sum += l1 * l3;
  }
};

template struct nonlocal_pw_op<float, psi::DEVICE_CPU>;
template struct nonlocal_pw_op<double, psi::DEVICE_CPU>;

}  // namespace hamilt