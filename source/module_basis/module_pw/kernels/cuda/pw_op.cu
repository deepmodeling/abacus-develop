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
__global__ void set_recip_to_real_output(
    const int nrxx,
    const bool add,
    const FPTYPE factor,
    const thrust::complex<FPTYPE>* in,
    FPTYPE* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= nrxx) {return;}
    if(add) {
        out[idx] += factor * in[idx].real();
    }
    else {
        out[idx] = in[idx].real();
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

template<class FPTYPE>
__global__ void set_real_to_recip_output_batch(
    const int total_elements,
    const int npwk,
    const int nxyz,
    const bool add,
    const FPTYPE factor,
    const int* box_index,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= total_elements) return;
    const int batch_idx = idx / npwk;
    const int element_idx = idx % npwk;
    const thrust::complex<FPTYPE>* batch_in = in + batch_idx * nxyz;
    thrust::complex<FPTYPE>* batch_out = out + batch_idx * npwk;
    
    const int box_idx = box_index[element_idx];  
    
    const FPTYPE scale_factor = factor / nxyz;
    const thrust::complex<FPTYPE> input_val = batch_in[box_idx];
    batch_out[element_idx] = add ? (batch_out[element_idx] + scale_factor * input_val) 
                                 : (scale_factor * input_val);
}

template<class FPTYPE>
__global__ void set_real_to_recip_output(
    const int npwk,
    const int nxyz,
    const bool add,
    const FPTYPE factor,
    const int* box_index,
    const thrust::complex<FPTYPE>* in,
    FPTYPE* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= npwk) {return;}
    if(add) {
        out[idx] += factor / nxyz * in[box_index[idx]].real();
    }
    else {
        out[idx] = in[box_index[idx]].real() / nxyz;
    }
}

template <typename FPTYPE>
void set_3d_fft_box_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const int npwk,
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

    cudaCheckOnDebug();
}

template <typename FPTYPE>
void set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const int nrxx,
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

    cudaCheckOnDebug();
}

template <typename FPTYPE>
void set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const int nrxx,
                                                                              const bool add,
                                                                              const FPTYPE factor,
                                                                              const std::complex<FPTYPE>* in,
                                                                              FPTYPE* out)
{
    const int block = (nrxx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_recip_to_real_output<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        nrxx,
        add,
        factor,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<FPTYPE*>(out));

    cudaCheckOnDebug();
}

template <typename FPTYPE>
void set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const int npwk,
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

    cudaCheckOnDebug();
}
template <typename FPTYPE>
void set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>::operator() (
        const int npwk,
        const int nxyz,
        const bool add,
        const FPTYPE factor,
        const int* box_index,
        const std::complex<FPTYPE>* in,
        std::complex<FPTYPE>* out,
        const int batch)
{
    const int total_elements = batch * npwk;
    const int grid_size = (total_elements + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_real_to_recip_output_batch<FPTYPE><<<grid_size, THREADS_PER_BLOCK>>>(
            total_elements,
            npwk,
            nxyz,
            add,
            factor,
            box_index,
            reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
            reinterpret_cast<thrust::complex<FPTYPE>*>(out));

    cudaCheckOnDebug();
}

template <typename FPTYPE>
void set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const int npwk,
                                                                              const int nxyz,
                                                                              const bool add,
                                                                              const FPTYPE factor,
                                                                              const int* box_index,
                                                                              const std::complex<FPTYPE>* in,
                                                                              FPTYPE* out)
{
    const int block = (npwk + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_real_to_recip_output<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        npwk,
        nxyz,
        add,
        factor,
        box_index,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<FPTYPE*>(out));

    cudaCheckOnDebug();
}

template struct set_3d_fft_box_op<float, base_device::DEVICE_GPU>;
template struct set_recip_to_real_output_op<float, base_device::DEVICE_GPU>;
template struct set_real_to_recip_output_op<float, base_device::DEVICE_GPU>;
template struct set_3d_fft_box_op<double, base_device::DEVICE_GPU>;
template struct set_recip_to_real_output_op<double, base_device::DEVICE_GPU>;
template struct set_real_to_recip_output_op<double, base_device::DEVICE_GPU>;

}  // namespace ModulePW
