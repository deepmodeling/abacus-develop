#include "source_basis/module_pw/kernels/pw_op.h"

#include <thrust/complex.h>
#include <cuda_runtime.h>
#include <base/macros/macros.h>
#include <type_traits>

namespace ModulePW {

#define THREADS_PER_BLOCK 256

template <class FPTYPE>
__device__ void atomic_add_complex(thrust::complex<FPTYPE>* address, const thrust::complex<FPTYPE>& value)
{
    static_assert(std::is_same<FPTYPE, float>::value || std::is_same<FPTYPE, double>::value,
                  "atomic_add_complex supports float and double only");
    FPTYPE* scalar_address = reinterpret_cast<FPTYPE*>(address);
    atomicAdd(scalar_address, value.real());
    atomicAdd(scalar_address + 1, value.imag());
}

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
__global__ void set_3d_fft_box_batch(
    const int npwk,
    const int nxyz,
    const int in_stride,
    const int* box_index,
    const int nbatch,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;
    for (int ib = 0; ib < nbatch; ib++)
    {
        for (int i = idx; i < npwk; i += stride)
        {
            int xx = box_index[i];
            out[xx] = in[i];

        }
        // Increment pointer
        out = out + nxyz;
        in = in + in_stride;
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
__global__ void set_recip_to_real_output_batch(
    const int nrxx,
    const int npw,
    const int nbatch,
    const int add,
    const FPTYPE* factor,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    if (add) {
        for (int ib = 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < nrxx; i += stride)
            {
                out[i] += factor[ib] * in[i];
            }
            out = out + nrxx;
            in = in + npw;

        }
    }
    else
    {
        for (int ib = 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < nrxx; i += stride)
            {
                out[i] =  in[i];
            }
            out = out + nrxx;
            in = in + npw;
        }
    }
}

template<class FPTYPE>
__global__ void set_recip_to_real_output_batch(
    const int nrxx,
    const int npw,
    const int nbatch,
    const int add,
    const FPTYPE factor,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    if (add) {
        for (int ib = 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < nrxx; i += stride)
            {
                out[i] += factor * in[i];
            }
            out = out + nrxx;
            in = in + npw;

        }
    }
    else
    {
        for (int ib = 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < nrxx; i += stride)
            {
                out[i] =  in[i];
            }
            out = out + nrxx;
            in = in + npw;
        }
    }
}

template<class FPTYPE>
__global__ void set_recip_to_real_output_batch(
    const int nrxx,
    const int npw,
    const int nbatch,
    const int add,
    const FPTYPE* factor,
    const thrust::complex<FPTYPE>* in,
    FPTYPE* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    if (add) {
        for (int ib = 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < nrxx; i += stride)
            {
                out[i] += factor[ib] * in[i].real();
            }
            out = out + nrxx;
            in = in + npw;
        }
    }
    else
    {
        for (int ib = 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < nrxx; i += stride)
            {
                out[i] =  in[i].real();
            }
            out = out + nrxx;
            in = in + npw;

        }
    }
}

template<class FPTYPE>
__global__ void set_recip_to_real_output_batch(
    const int nrxx,
    const int npw,
    const int nbatch,
    const int add,
    const FPTYPE factor,
    const thrust::complex<FPTYPE>* in,
    FPTYPE* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    if (add) {
        for (int ib = 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < nrxx; i += stride)
            {
                out[i] += factor * in[i].real();
            }
            out = out + nrxx;
            in = in + npw;
        }
    }
    else
    {
        for (int ib = 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < nrxx; i += stride)
            {
                out[i] =  in[i].real();
            }
            out = out + nrxx;
            in = in + npw;

        }
    }
}


template<class FPTYPE>
__global__ void set_real_to_recip_output(
    const int npwk,
    const FPTYPE one_over_nxyz,
    const bool add,
    const FPTYPE factor,
    const int* box_index,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= npwk) {return;}
    if(add) {
        out[idx] += factor * one_over_nxyz * in[box_index[idx]];
    }
    else {
        out[idx] = in[box_index[idx]] * one_over_nxyz;
    }
}

template<class FPTYPE>
__global__ void set_real_to_recip_output(
    const int npwk,
    const FPTYPE one_over_nxyz,
    const bool add,
    const FPTYPE factor,
    const int* box_index,
    const thrust::complex<FPTYPE>* in,
    FPTYPE* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= npwk) {return;}
    if(add) {
        out[idx] += factor * one_over_nxyz * in[box_index[idx]].real();
    }
    else {
        out[idx] = in[box_index[idx]].real() * one_over_nxyz;
    }
}

template<class FPTYPE>
__global__ void set_real_to_recip_output_batch(
    const int npwk,
    const int nxyz,
    const int out_stride,
    const FPTYPE one_over_nxyz,
    const int nbatch,
    const int add,
    const FPTYPE* factor,
    const int* box_index,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    if (add) {
        for (int ib= 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < npwk; i += stride)
            {
                out[i] += factor[ib] * one_over_nxyz * in[box_index[i]];
            }
            in = in + nxyz;
            out = out + out_stride;

        }
    }
    else
    {
        for (int ib= 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < npwk; i += stride)
            {
                out[i] = in[box_index[i]] * one_over_nxyz;
            }
            in = in + nxyz;
            out = out + out_stride;
        }
    }
}

template<class FPTYPE>
__global__ void set_real_to_recip_output_batch(
    const int npwk,
    const int nxyz,
    const int out_stride,
    const FPTYPE one_over_nxyz,
    const int nbatch,
    const int add,
    const FPTYPE factor,
    const int* box_index,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    if (add) {
        for (int ib= 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < npwk; i += stride)
            {
                out[i] += factor * one_over_nxyz * in[box_index[i]];
            }
            in = in + nxyz;
            out = out + out_stride;

        }
    }
    else
    {
        for (int ib= 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < npwk; i += stride)
            {
                out[i] = in[box_index[i]] * one_over_nxyz;
            }
            in = in + nxyz;
            out = out + out_stride;
        }
    }
}

template<class FPTYPE>
__global__ void set_real_to_recip_output_batch(
    const int npwk,
    const int nxyz,
    const int out_stride,
    const FPTYPE one_over_nxyz,
    const int nbatch,
    const int add,
    const FPTYPE* factor,
    const int* box_index,
    const thrust::complex<FPTYPE>* in,
    FPTYPE* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    if (add) {
        for (int ib= 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < npwk; i += stride)
            {
                out[i] += factor[ib] * one_over_nxyz * in[box_index[i]].real();
            }
            in = in + nxyz;
            out = out + out_stride;

        }
    }
    else
    {
        for (int ib= 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < npwk; i += stride)
            {
                out[i] = in[box_index[i]].real() * one_over_nxyz;
            }
            in = in + nxyz;
            out = out + out_stride;
        }
    }
}

template<class FPTYPE>
__global__ void set_real_to_recip_output_batch(
    const int npwk,
    const int nxyz,
    const int out_stride,
    const FPTYPE one_over_nxyz,
    const int nbatch,
    const int add,
    const FPTYPE factor,
    const int* box_index,
    const thrust::complex<FPTYPE>* in,
    FPTYPE* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    if (add) {
        for (int ib= 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < npwk; i += stride)
            {
                out[i] += factor * one_over_nxyz * in[box_index[i]].real();
            }
            in = in + nxyz;
            out = out + out_stride;

        }
    }
    else
    {
        for (int ib= 0; ib < nbatch; ib++)
        {
            for (int i = idx; i < npwk; i += stride)
            {
                out[i] = in[box_index[i]].real() * one_over_nxyz;
            }
            in = in + nxyz;
            out = out + out_stride;
        }
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

    CHECK_CUDA_SYNC();
}

template <typename FPTYPE>
void set_3d_fft_box_op<FPTYPE, base_device::DEVICE_GPU>::operator_batch(const int npwk,
                                                                    const int nxyz,
                                                                    const int in_stride,
                                                                    const int* box_index,
                                                                    const int nbatch,
                                                                    const std::complex<FPTYPE>* in,
                                                                    std::complex<FPTYPE>* out)
{
    const int block = (npwk + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_3d_fft_box_batch<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        npwk,
        nxyz,
        in_stride,
        box_index,
        nbatch,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));

    CHECK_CUDA_SYNC();
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

    CHECK_CUDA_SYNC();
}

template <typename FPTYPE>
void set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>::operator_batch(const int nrxx,
                                                                              const int npw,
                                                                              const int nbatch,
                                                                              const bool add,
                                                                              const FPTYPE* factor,
                                                                              const std::complex<FPTYPE>* in,
                                                                              std::complex<FPTYPE>* out)
{
    const int block = (nrxx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_recip_to_real_output_batch<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        nrxx,
        npw,
        nbatch,
        add,
        factor,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));

    CHECK_CUDA_SYNC();
}

template <typename FPTYPE>
void set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>::operator_batch(const int nrxx,
                                                                              const int npw,
                                                                              const int nbatch,
                                                                              const bool add,
                                                                              const FPTYPE factor,
                                                                              const std::complex<FPTYPE>* in,
                                                                              std::complex<FPTYPE>* out)
{
    const int block = (nrxx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_recip_to_real_output_batch<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        nrxx,
        npw,
        nbatch,
        add,
        factor,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));

    CHECK_CUDA_SYNC();
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

    CHECK_CUDA_SYNC();
}

template <typename FPTYPE>
void set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>::operator_batch(const int nrxx,
                                                                              const int npw,
                                                                              const int nbatch,
                                                                              const bool add,
                                                                              const FPTYPE* factor,
                                                                              const std::complex<FPTYPE>* in,
                                                                              FPTYPE* out)
{
    const int block = (nrxx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_recip_to_real_output_batch<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        nrxx,
        npw,
        nbatch,
        add,
        factor,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<FPTYPE*>(out)
    );

    CHECK_CUDA_SYNC();
}

template <typename FPTYPE>
void set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>::operator_batch(const int nrxx,
                                                                              const int npw,
                                                                              const int nbatch,
                                                                              const bool add,
                                                                              const FPTYPE factor,
                                                                              const std::complex<FPTYPE>* in,
                                                                              FPTYPE* out)
{
    const int block = (nrxx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_recip_to_real_output_batch<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        nrxx,
        npw,
        nbatch,
        add,
        factor,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<FPTYPE*>(out)
    );

    CHECK_CUDA_SYNC();
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
        static_cast<FPTYPE>(1.0/nxyz),
        add,
        factor,
        box_index,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));

    CHECK_CUDA_SYNC();
}

template <typename FPTYPE>
void set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>::operator_batch(const int npwk,
                                                                              const int nxyz,
                                                                              const int out_stride,
                                                                              const int nbatch,
                                                                              const bool add,
                                                                              const FPTYPE* factor,
                                                                              const int* box_index,
                                                                              const std::complex<FPTYPE>* in,
                                                                              std::complex<FPTYPE>* out)
{
    const int block = (npwk + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_real_to_recip_output_batch<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        npwk,
        nxyz,
        out_stride,
        static_cast<FPTYPE>(1.0/nxyz),
        nbatch,
        add,
        factor,
        box_index,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));

    CHECK_CUDA_SYNC();
}

template <typename FPTYPE>
void set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>::operator_batch(const int npwk,
                                                                              const int nxyz,
                                                                              const int out_stride,
                                                                              const int nbatch,
                                                                              const bool add,
                                                                              const FPTYPE factor,
                                                                              const int* box_index,
                                                                              const std::complex<FPTYPE>* in,
                                                                              std::complex<FPTYPE>* out)
{
    const int block = (npwk + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_real_to_recip_output_batch<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        npwk,
        nxyz,
        out_stride,
        static_cast<FPTYPE>(1.0/nxyz),
        nbatch,
        add,
        factor,
        box_index,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));

    CHECK_CUDA_SYNC();
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
        static_cast<FPTYPE>(1.0/nxyz),
        add,
        factor,
        box_index,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<FPTYPE*>(out));

    CHECK_CUDA_SYNC();
}


template <typename FPTYPE>
void set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>::operator_batch(const int npwk,
                                                                              const int nxyz,
                                                                              const int out_stride,
                                                                              const int nbatch,
                                                                              const bool add,
                                                                              const FPTYPE* factor,
                                                                              const int* box_index,
                                                                              const std::complex<FPTYPE>* in,
                                                                              FPTYPE* out)
{
    const int block = (npwk + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_real_to_recip_output_batch<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        npwk,
        nxyz,
        out_stride,
        static_cast<FPTYPE>(1.0/nxyz),
        nbatch,
        add,
        factor,
        box_index,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<FPTYPE*>(out));

    CHECK_CUDA_SYNC();
}

template <typename FPTYPE>
void set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>::operator_batch(const int npwk,
                                                                              const int nxyz,
                                                                              const int out_stride,
                                                                              const int nbatch,
                                                                              const bool add,
                                                                              const FPTYPE factor,
                                                                              const int* box_index,
                                                                              const std::complex<FPTYPE>* in,
                                                                              FPTYPE* out)
{
    const int block = (npwk + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_real_to_recip_output_batch<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        npwk,
        nxyz,
        out_stride,
        static_cast<FPTYPE>(1.0/nxyz),
        nbatch,
        add,
        factor,
        box_index,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<FPTYPE*>(out));

    CHECK_CUDA_SYNC();
}

template<class FPTYPE>
__global__ void set_3d_fft_box_remapped(
    const int npw_full,
    const int* rep_igl,
    const int* box_index,
    const thrust::complex<FPTYPE>* phase,
    const bool conjugate,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= npw_full) {return;}
    const thrust::complex<FPTYPE> value = conjugate ? thrust::conj(in[rep_igl[idx]]) * thrust::conj(phase[idx])
                                                    : in[rep_igl[idx]] * phase[idx];
    out[box_index[idx]] = value;
}

template<class FPTYPE>
__global__ void set_3d_fft_box_remapped_batch(
    const int npw_full,
    const int nxyz,
    const int in_stride,
    const int* rep_igl,
    const int* box_index,
    const thrust::complex<FPTYPE>* phase,
    const bool conjugate,
    const int nbatch,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int total = npw_full * nbatch;
    if (idx >= total) {return;}
    const int ib = idx / npw_full;
    const int ig = idx - ib * npw_full;
    const thrust::complex<FPTYPE>* in_ib = in + static_cast<std::size_t>(ib) * in_stride;
    thrust::complex<FPTYPE>* out_ib = out + static_cast<std::size_t>(ib) * nxyz;
    const thrust::complex<FPTYPE> value = conjugate ? thrust::conj(in_ib[rep_igl[ig]]) * thrust::conj(phase[ig])
                                                    : in_ib[rep_igl[ig]] * phase[ig];
    out_ib[box_index[ig]] = value;
}

template<class FPTYPE>
__global__ void zero_remapped_recip_output(
    const int npw_full,
    const int* rep_igl,
    thrust::complex<FPTYPE>* out)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < npw_full)
    {
        out[rep_igl[idx]] = thrust::complex<FPTYPE>(0, 0);
    }
}

template<class FPTYPE>
__global__ void set_real_to_recip_remapped_output(
    const int npw_full,
    const int nxyz,
    const FPTYPE factor,
    const int* rep_igl,
    const int* box_index,
    const thrust::complex<FPTYPE>* phase,
    const bool conjugate,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= npw_full) {return;}
    const thrust::complex<FPTYPE> value = conjugate ? thrust::conj(in[box_index[idx]]) * phase[idx]
                                                    : in[box_index[idx]] * thrust::conj(phase[idx]);
    atomic_add_complex(out + rep_igl[idx], factor * value / static_cast<FPTYPE>(nxyz));
}

template <typename FPTYPE>
void set_3d_fft_box_remapped_op<FPTYPE, base_device::DEVICE_GPU>::operator()(
    const int npw_full,
    const int* rep_igl,
    const int* box_index,
    const std::complex<FPTYPE>* phase,
    const bool conjugate,
    const std::complex<FPTYPE>* in,
    std::complex<FPTYPE>* out)
{
    const int block = (npw_full + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_3d_fft_box_remapped<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        npw_full,
        rep_igl,
        box_index,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(phase),
        conjugate,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));

    CHECK_CUDA_SYNC();
}

template <typename FPTYPE>
void set_3d_fft_box_remapped_op<FPTYPE, base_device::DEVICE_GPU>::operator_batch(
    const int npw_full,
    const int nxyz,
    const int in_stride,
    const int* rep_igl,
    const int* box_index,
    const std::complex<FPTYPE>* phase,
    const bool conjugate,
    const int nbatch,
    const std::complex<FPTYPE>* in,
    std::complex<FPTYPE>* out)
{
    const int total = npw_full * nbatch;
    const int block = (total + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_3d_fft_box_remapped_batch<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        npw_full,
        nxyz,
        in_stride,
        rep_igl,
        box_index,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(phase),
        conjugate,
        nbatch,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));

    CHECK_CUDA_SYNC();
}

template <typename FPTYPE>
void set_real_to_recip_remapped_output_op<FPTYPE, base_device::DEVICE_GPU>::operator()(
    const int npw_full,
    const int nxyz,
    const bool add,
    const FPTYPE factor,
    const int* rep_igl,
    const int* box_index,
    const std::complex<FPTYPE>* phase,
    const bool conjugate,
    const std::complex<FPTYPE>* in,
    std::complex<FPTYPE>* out)
{
    const int block = (npw_full + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    if (!add)
    {
        zero_remapped_recip_output<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
            npw_full,
            rep_igl,
            reinterpret_cast<thrust::complex<FPTYPE>*>(out));
        CHECK_CUDA_SYNC();
    }
    const FPTYPE tmpfac = add ? factor : FPTYPE(1.0);
    set_real_to_recip_remapped_output<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        npw_full,
        nxyz,
        tmpfac,
        rep_igl,
        box_index,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(phase),
        conjugate,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));

    CHECK_CUDA_SYNC();
}



template struct set_3d_fft_box_op<float, base_device::DEVICE_GPU>;
template struct set_recip_to_real_output_op<float, base_device::DEVICE_GPU>;
template struct set_real_to_recip_output_op<float, base_device::DEVICE_GPU>;
template struct set_3d_fft_box_remapped_op<float, base_device::DEVICE_GPU>;
template struct set_real_to_recip_remapped_output_op<float, base_device::DEVICE_GPU>;
template struct set_3d_fft_box_op<double, base_device::DEVICE_GPU>;
template struct set_recip_to_real_output_op<double, base_device::DEVICE_GPU>;
template struct set_real_to_recip_output_op<double, base_device::DEVICE_GPU>;
template struct set_3d_fft_box_remapped_op<double, base_device::DEVICE_GPU>;
template struct set_real_to_recip_remapped_output_op<double, base_device::DEVICE_GPU>;

}  // namespace ModulePW
