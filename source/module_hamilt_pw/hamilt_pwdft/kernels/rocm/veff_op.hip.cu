#include "module_hamilt_pw/hamilt_pwdft/kernels/veff_op.h"

#include <complex>
#include <thrust/complex.h>

#include <hip/hip_runtime.h>
#include <base/macros/macros.h>

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
__global__ void veff_pw(
    const int size,
    thrust::complex<FPTYPE>* out,
    thrust::complex<FPTYPE>* out1,
    const FPTYPE* in)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size) {return;}
    thrust::complex<FPTYPE> sup = 
        out[idx] * (in[0 * size + idx] + in[3 * size + idx])
            + out1[idx] * (in[1 * size + idx] - thrust::complex<FPTYPE>(0.0, 1.0) * in[2 * size + idx]);
    thrust::complex<FPTYPE> sdown = 
        out1[idx] * (in[0 * size + idx] - in[3 * size + idx])
            + out[idx] * (in[1 * size + idx] + thrust::complex<FPTYPE>(0.0, 1.0) * in[2 * size + idx]);
    out[idx] = sup;
    out1[idx] = sdown;
}

template <typename FPTYPE>
void veff_pw_op<FPTYPE, psi::DEVICE_GPU>::operator() (
    const psi::DEVICE_GPU* dev,
    const int& size,
    std::complex<FPTYPE>* out,
    const FPTYPE* in)
{
    const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(veff_pw<FPTYPE>), dim3(block), dim3(THREADS_PER_BLOCK), 0, 0, 
        size, // control params
        reinterpret_cast<thrust::complex<FPTYPE>*>(out), // array of data
        in); // array of data

    hipErrcheck(hipGetLastError());
    hipErrcheck(hipDeviceSynchronize());
}

template <typename FPTYPE>
void veff_pw_op<FPTYPE, psi::DEVICE_GPU>::operator() (
    const psi::DEVICE_GPU* dev,
    const int& size,
    std::complex<FPTYPE>* out,
    std::complex<FPTYPE>* out1,
    const FPTYPE** in)
{
    const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(veff_pw<FPTYPE>), dim3(block), dim3(THREADS_PER_BLOCK), 0, 0, 
        size, // control params
        reinterpret_cast<thrust::complex<FPTYPE>*>(out), // array of data
        reinterpret_cast<thrust::complex<FPTYPE>*>(out1), // array of data
        in[0]); // array of data
    
    hipErrcheck(hipGetLastError());
    hipErrcheck(hipDeviceSynchronize());
}

template struct veff_pw_op<float, psi::DEVICE_GPU>;
template struct veff_pw_op<double, psi::DEVICE_GPU>;

}  // namespace hamilt