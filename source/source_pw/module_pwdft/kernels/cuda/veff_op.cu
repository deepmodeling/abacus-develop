#include "source_pw/module_pwdft/kernels/veff_op.h"

#include <complex>

#include <cuda_runtime.h>
#include <thrust/complex.h>
#include <base/macros/macros.h>

namespace hamilt {

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
    thrust::complex<FPTYPE>* in)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size) {return;}
    const int base = idx * 4;
    thrust::complex<FPTYPE> sup = out[idx] * in[base] + out1[idx] * in[base+1];
    thrust::complex<FPTYPE> sdown = out1[idx] * in[base+2] + out[idx] * in[base+3];
    out[idx] = sup;
    out1[idx] = sdown;
}

template <typename FPTYPE>
__global__ void rearrange_op(
    const int size,
    const FPTYPE* in,
    thrust::complex<FPTYPE>* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size) {return;}
    const int base = idx * 4;
    const FPTYPE part_1 = in[idx];
    const FPTYPE part_2 = in[idx + size];
    const FPTYPE part_3 = in[idx + 2 * size];
    const FPTYPE part_4 = in[idx + 3 * size];
    out[base] = thrust::complex<FPTYPE>(part_1 + part_4, 0.0);
    out[base + 1] = thrust::complex<FPTYPE>(part_2 , -part_3);
    out[base + 2] = thrust::complex<FPTYPE>(part_1 - part_4, 0.0);
    out[base + 3] = thrust::complex<FPTYPE>(part_2, part_3);

}
template <typename FPTYPE>
void rearrange<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* device, 
                    const int& size, 
                    const FPTYPE* in, 
                    std::complex<FPTYPE>* out) const
{
    const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    rearrange_op<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        size, // control params
        in, // array of data
        reinterpret_cast<thrust::complex<FPTYPE>*>(out)); // array of data   
    cudaCheckOnDebug();
}

template <typename FPTYPE>
void veff_pw_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* dev,
                                                             const int& size,
                                                             std::complex<FPTYPE>* out,
                                                             const FPTYPE* in)
{
    const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    veff_pw<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        size, // control params
        reinterpret_cast<thrust::complex<FPTYPE>*>(out), // array of data
        in); // array of data

    cudaCheckOnDebug();
}

template <typename FPTYPE>
void veff_pw_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* dev,
                                                             const int& size,
                                                             std::complex<FPTYPE>* out,
                                                             std::complex<FPTYPE>* out1,
                                                             std::complex<FPTYPE>* in)
{
    const int block = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    veff_pw<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        size, // control params
        reinterpret_cast<thrust::complex<FPTYPE>*>(out), // array of data
        reinterpret_cast<thrust::complex<FPTYPE>*>(out1), // array of data
        reinterpret_cast<thrust::complex<FPTYPE>*>(in)); // array of data

    cudaCheckOnDebug();
}

template struct rearrange<float, base_device::DEVICE_GPU>;
template struct rearrange<double, base_device::DEVICE_GPU>;
template struct veff_pw_op<float, base_device::DEVICE_GPU>;
template struct veff_pw_op<double, base_device::DEVICE_GPU>;

}  // namespace hamilt