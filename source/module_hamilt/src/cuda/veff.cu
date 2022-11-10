#include "module_hamilt/include/veff.h"
#include <complex>
#include <thrust/complex.h>
#include "cuda_runtime.h"

namespace hamilt{

#define THREADS_PER_BLOCK 256

template <typename FPTYPE>
__global__ void veff_pw(
    const int size,
    thrust::complex<FPTYPE>* out,
    const FPTYPE* in)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size) {return;}
    out[idx] *= in[idx];
}

template <typename FPTYPE>
void veff_pw_op<FPTYPE, psi::DEVICE_GPU>::operator() (
    const psi::DEVICE_GPU* dev,
    const int& size,
    std::complex<FPTYPE>* out,
    const FPTYPE* in)
{
    const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    veff_pw<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        size, // control params
        reinterpret_cast<thrust::complex<FPTYPE>*>(out), // array of data
        in); // array of data
    // cpu part:
    // for (int ir = 0; ir < size; ++ir)
    // {
    //     out[ir] *= in[ir];
    // }
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
}

template struct veff_pw_op<double, psi::DEVICE_GPU>;

}  // namespace hamilt