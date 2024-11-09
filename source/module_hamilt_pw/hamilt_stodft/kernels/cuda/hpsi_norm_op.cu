#include "module_hamilt_pw/hamilt_stodft/kernels/hpsi_norm_op.h"

#include <base/macros/macros.h>
#include <cuda_runtime.h>
#include <thrust/complex.h>

namespace hamilt
{
#define THREADS_PER_BLOCK 256

template <typename FPTYPE>
__global__ void hpsi_norm(const int npwk_max,
                          const int npwk,
                          const FPTYPE Ebar,
                          const FPTYPE DeltaE,
                          thrust::complex<FPTYPE>* hpsi,
                          const thrust::complex<FPTYPE>* psi_in)
{
    const int block_idx = blockIdx.x;
    const int thread_idx = threadIdx.x;
    const int start_idx = block_idx * npwk_max;
    for (int ii = thread_idx; ii < npwk; ii += blockDim.x)
    {
        hpsi[start_idx + ii] = (hpsi[start_idx + ii] - Ebar * psi_in[start_idx + ii]) / DeltaE;
    }
}

template <typename FPTYPE>
void hamilt::hpsi_norm_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* dev,
                                                                       const int& nbands,
                                                                       const int& npwk_max,
                                                                       const int& npwk,
                                                                       const FPTYPE& Ebar,
                                                                       const FPTYPE& DeltaE,
                                                                       std::complex<FPTYPE>* hpsi,
                                                                       const std::complex<FPTYPE>* psi_in)
{
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    hpsi_norm<FPTYPE><<<nbands, THREADS_PER_BLOCK>>>(
      npwk_max, npwk, Ebar, DeltaE,
      reinterpret_cast<thrust::complex<FPTYPE>*>(hpsi),
      reinterpret_cast<const thrust::complex<FPTYPE>*>(psi_in));
    cudaCheckOnDebug();
}

template struct hpsi_norm_op<float, base_device::DEVICE_GPU>;
template struct hpsi_norm_op<double, base_device::DEVICE_GPU>;

} // namespace hamilt