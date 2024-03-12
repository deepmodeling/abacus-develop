#include "module_basis/module_pw/kernels/pw_op.h"

#include <thrust/complex.h>
#include <cuda_runtime.h>
#include <base/macros/macros.h>

namespace ModulePW {

#define THREADS_PER_BLOCK 256

template<class FPTYPE>
__global__ void set_3d_fft_box(
    const int npwk,
    const int* box_index,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < npwk)
    {
        int xx = box_index[idx];
        out[xx] = in[idx];
    }
}

template<class FPTYPE>
__global__ void set_recip_to_real_output(
    const int nrxx,
    const bool add,
    const FPTYPE factor,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= nrxx) {return;}
    if(add) {
        out[idx] += factor * in[idx];
    }
    else {
        out[idx] = in[idx];
    }
}

template<class FPTYPE>
__global__ void set_real_to_recip_output(
    const int npwk,
    const int nxyz,
    const bool add,
    const FPTYPE factor,
    const int* box_index,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= npwk) {return;}
    if(add) {
        out[idx] += factor / nxyz * in[box_index[idx]];
    }
    else {
        out[idx] = in[box_index[idx]] / nxyz;
    }
}

template <typename FPTYPE>
void set_3d_fft_box_op<FPTYPE, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU*  /*dev*/,
    const int npwk,
    const int* box_index,
    const std::complex<FPTYPE>* in,
    std::complex<FPTYPE>* out)
{
    const int block = (npwk + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_3d_fft_box<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        npwk,
        box_index,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));
    
    cudaErrcheck(cudaGetLastError());
    cudaErrcheck(cudaDeviceSynchronize());
}

template <typename FPTYPE>
void set_recip_to_real_output_op<FPTYPE, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU*  /*dev*/,
    const int nrxx,
    const bool add,
    const FPTYPE factor,
    const std::complex<FPTYPE>* in,
    std::complex<FPTYPE>* out)
{
    const int block = (nrxx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_recip_to_real_output<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        nrxx,
        add,
        factor,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));
    
    cudaErrcheck(cudaGetLastError());
    cudaErrcheck(cudaDeviceSynchronize());
}

template <typename FPTYPE>
void set_real_to_recip_output_op<FPTYPE, psi::DEVICE_GPU>::operator()(
    const psi::DEVICE_GPU*  /*dev*/,
    const int npwk,
    const int nxyz,
    const bool add,
    const FPTYPE factor,
    const int* box_index,
    const std::complex<FPTYPE>* in,
    std::complex<FPTYPE>* out)
{
    const int block = (npwk + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_real_to_recip_output<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        npwk,
        nxyz,
        add,
        factor,
        box_index,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));
    
    cudaErrcheck(cudaGetLastError());
    cudaErrcheck(cudaDeviceSynchronize());
}

template struct set_3d_fft_box_op<float, psi::DEVICE_GPU>;
template struct set_recip_to_real_output_op<float, psi::DEVICE_GPU>;
template struct set_real_to_recip_output_op<float, psi::DEVICE_GPU>;
template struct set_3d_fft_box_op<double, psi::DEVICE_GPU>;
template struct set_recip_to_real_output_op<double, psi::DEVICE_GPU>;
template struct set_real_to_recip_output_op<double, psi::DEVICE_GPU>;

}  // namespace ModulePW

