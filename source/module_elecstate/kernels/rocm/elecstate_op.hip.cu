#include "module_elecstate/kernels/elecstate_op.h"
#include <thrust/complex.h>

#include <hip/hip_runtime.h>
#include <base/macros/macros.h>

#define THREADS_PER_BLOCK 256

namespace elecstate {

template<typename FPTYPE>
__global__ void elecstate_pw(
    const int spin,
    const int nrxx,
    const FPTYPE w1,
    FPTYPE* rho,
    const thrust::complex<FPTYPE>* wfcr)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if(idx >= nrxx) {return;}
  rho[spin * nrxx + idx] += w1 * norm(wfcr[idx]);
}

template<typename FPTYPE>
__global__ void elecstate_pw(
    const bool DOMAG,
    const bool DOMAG_Z,
    const int nrxx,
    const FPTYPE w1,
    FPTYPE* rho,
    const thrust::complex<FPTYPE>* wfcr,
    const thrust::complex<FPTYPE>* wfcr_another_spin)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if(idx >= nrxx) {return;}
  rho[0 * nrxx + idx] += w1 * (norm(wfcr[idx]) + norm(wfcr_another_spin[idx]));

  if (DOMAG) {
    rho[1 * nrxx + idx] += w1 * 2.0
                  * (wfcr[idx].real() * wfcr_another_spin[idx].real()
                  +  wfcr[idx].imag() * wfcr_another_spin[idx].imag());
    rho[2 * nrxx + idx] += w1 * 2.0
                  * (wfcr[idx].real() * wfcr_another_spin[idx].imag()
                  - wfcr_another_spin[idx].real() * wfcr[idx].imag());
    rho[3 * nrxx + idx] += w1 * (norm(wfcr[idx]) - norm(wfcr_another_spin[idx]));
  }
  else if(DOMAG_Z) {
    rho[1 * nrxx + idx] = 0;
    rho[2 * nrxx + idx] = 0;
    rho[3 * nrxx + idx] += w1 * (norm(wfcr[idx]) - norm(wfcr_another_spin[idx]));
  }
  else {
    rho[0 * nrxx + idx] = 0;
    rho[1 * nrxx + idx] = 0;
    rho[2 * nrxx + idx] = 0;
    rho[3 * nrxx + idx] = 0;
  }
}

template <typename FPTYPE>
void elecstate_pw_op<FPTYPE, psi::DEVICE_GPU>::operator() (
    const psi::DEVICE_GPU* ctx,
    const int& spin,
    const int& nrxx,
    const FPTYPE& w1,
    FPTYPE** rho,
    const std::complex<FPTYPE>* wfcr)
{
  const int block = (nrxx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  hipLaunchKernelGGL(HIP_KERNEL_NAME(elecstate_pw<FPTYPE>), dim3(block), dim3(THREADS_PER_BLOCK), 0, 0, 
    spin, nrxx, w1, rho[0], 
    reinterpret_cast<const thrust::complex<FPTYPE>*>(wfcr)
  );

  hipErrcheck(hipGetLastError());
  hipErrcheck(hipDeviceSynchronize());
}

template <typename FPTYPE>
void elecstate_pw_op<FPTYPE, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU* ctx,
    const bool& DOMAG,
    const bool& DOMAG_Z,
    const int& nrxx,
    const FPTYPE& w1,
    FPTYPE** rho,
    const std::complex<FPTYPE>* wfcr,
    const std::complex<FPTYPE>* wfcr_another_spin)
{
  const int block = (nrxx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  hipLaunchKernelGGL(HIP_KERNEL_NAME(elecstate_pw<FPTYPE>), dim3(block), dim3(THREADS_PER_BLOCK), 0, 0, 
    DOMAG, DOMAG_Z, nrxx, w1, rho[0], 
    reinterpret_cast<const thrust::complex<FPTYPE>*>(wfcr), 
    reinterpret_cast<const thrust::complex<FPTYPE>*>(wfcr_another_spin)
  );

  hipErrcheck(hipGetLastError());
  hipErrcheck(hipDeviceSynchronize());
}

template struct elecstate_pw_op<float, psi::DEVICE_GPU>;
template struct elecstate_pw_op<double, psi::DEVICE_GPU>;
}