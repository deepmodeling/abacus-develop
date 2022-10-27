#include "module_hamilt/include/nonlocal.h"
#include <complex>
#include <thrust/complex.h>

using namespace hamilt; 

#define THREADS_PER_BLOCK 256

template <typename FPTYPE>
__global__ void nonlocal_pw(
    const int l1,
    const int l2,
    const int l3,
    const int sum,
    const int iat,
    const int spin,
    const int nkb,
    const int deeq_x,
    const int deeq_y,
    const int deeq_z,
    const FPTYPE* deeq,
    thrust::complex<FPTYPE>* ps,
    const thrust::complex<FPTYPE>* becp)
{
  const int ii = blockIdx.x / l2;
  const int jj = blockIdx.x % l2;
  for (int kk = threadIdx.x; kk < l3; kk += blockDim.x) {
    thrust::complex<FPTYPE> res(0.0, 0.0);
    for (int xx = 0; xx < l3; xx++) {
      res
        += deeq[((spin * deeq_x + iat + ii) * deeq_y + xx) * deeq_z + kk] 
        *  becp[jj * nkb + sum + ii * l3 + xx];
    }
    ps[(sum + ii * l3 + kk) * l2 + jj] += res;
  }
}

template <typename FPTYPE> 
void hamilt::nonlocal_pw_op<FPTYPE, psi::DEVICE_GPU>::operator() (
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
    const std::complex<FPTYPE>* becp)
{
  // denghui implement 20221019
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  nonlocal_pw<FPTYPE><<<l1 * l2, THREADS_PER_BLOCK>>>(
    l1, l2, l3, // loop size
    sum, iat, spin, nkb,   // control params
    deeq_x, deeq_y, deeq_z, deeq,  // deeq realArray operator()
    reinterpret_cast<thrust::complex<FPTYPE>*>(ps), // array of data
    reinterpret_cast<const thrust::complex<FPTYPE>*>(becp)); // array of data
  iat += l1;
  sum += l1 * l3;
  // for (int ii = 0; ii < l1; ii++) {
  //   // each atom has nproj, means this is with structure factor;
  //   // each projector (each atom) must multiply coefficient
  //   // with all the other projectors.
  //   for (int jj = 0; jj < l2; ++jj) 
  //     for (int kk = 0; kk < l3; kk++) 
  //       for (int xx = 0; xx < l3; xx++) 
  //         ps[(sum + kk) * l2 + jj]
  //             += deeq[((current_spin * deeq_x + iat) * deeq_y + xx) * deeq_z + kk] 
  //             *  becp[jj * nkb + sum + xx];
  //   sum += l3;
  //   ++iat;
  // }
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
}

namespace hamilt{
template struct nonlocal_pw_op<double, psi::DEVICE_GPU>;
}