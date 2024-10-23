#include "module_hamilt_pw/hamilt_pwdft/kernels/ekinetic_op.h"

#include <complex>
#include <thrust/complex.h>

#include <hip/hip_runtime.h>
#include <base/macros/macros.h>

using namespace hamilt;
#define THREADS_PER_BLOCK 256

template <typename FPTYPE>
__global__ void ekinetic_pw(
    const int npw,
    const int max_npw,
    const bool is_first_node,
    const FPTYPE tpiba2,
    const FPTYPE* gk2,
    thrust::complex<FPTYPE>* hpsi,
    const thrust::complex<FPTYPE>* psi)
{
  const int block_idx = blockIdx.x;
  const int thread_idx = threadIdx.x;
  const int start_idx = block_idx * max_npw;
  if(is_first_node)
  {
      for (int ii = thread_idx; ii < npw; ii += blockDim.x)
      {
          hpsi[start_idx + ii] = gk2[ii] * tpiba2 * psi[start_idx + ii];
      }
      for (int ii = npw + thread_idx; ii < max_npw; ii += blockDim.x)
      {
          hpsi[start_idx + ii] = 0.0;
      }
  }
  else
  {
      for (int ii = thread_idx; ii < npw; ii += blockDim.x)
      {
          hpsi[start_idx + ii] += gk2[ii] * tpiba2 * psi[start_idx + ii];
      }
  }
  
}

template <typename FPTYPE>
void hamilt::ekinetic_pw_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* dev,
                                                                         const int& nband,
                                                                         const int& npw,
                                                                         const int& max_npw,
                                                                         const bool& is_first_node,
                                                                         const FPTYPE& tpiba2,
                                                                         const FPTYPE* gk2_ik,
                                                                         std::complex<FPTYPE>* tmhpsi,
                                                                         const std::complex<FPTYPE>* tmpsi_in)
{
  // denghui implement 20221019
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  hipLaunchKernelGGL(HIP_KERNEL_NAME(ekinetic_pw<FPTYPE>), dim3(nband), dim3(THREADS_PER_BLOCK), 0, 0,
    npw, max_npw, is_first_node, tpiba2, // control params
    gk2_ik, // array of data
    reinterpret_cast<thrust::complex<FPTYPE>*>(tmhpsi), // array of data
    reinterpret_cast<const thrust::complex<FPTYPE>*>(tmpsi_in)); // array of data

   hipCheckOnDebug();
}

namespace hamilt{
template struct ekinetic_pw_op<float, base_device::DEVICE_GPU>;
template struct ekinetic_pw_op<double, base_device::DEVICE_GPU>;
}