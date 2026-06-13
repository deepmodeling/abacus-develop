#include "source_pw/module_pwdft/kernels/exx_q_state_op.h"
#include "source_base/module_device/device_check.h"

#include <stdexcept>
#include <string>
#include <thrust/complex.h>

namespace hamilt
{
template <typename FPTYPE>
__global__ void exx_conjugate_real_kernel(const thrust::complex<FPTYPE>* in,
                                          thrust::complex<FPTYPE>* out,
                                          int nrxx)
{
    const int ir = blockIdx.x * blockDim.x + threadIdx.x;
    if (ir < nrxx)
    {
        out[ir] = thrust::conj(in[ir]);
    }
}

template <typename FPTYPE>
struct exx_conjugate_real_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>
{
    using T = std::complex<FPTYPE>;
    void operator()(const T* in, T* out, int nrxx)
    {
        const int threads_per_block = 256;
        const int num_blocks = (nrxx + threads_per_block - 1) / threads_per_block;
        exx_conjugate_real_kernel<FPTYPE><<<num_blocks, threads_per_block>>>(
            reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
            reinterpret_cast<thrust::complex<FPTYPE>*>(out),
            nrxx);

        CHECK_LAST_CUDA_ERROR("exx_conjugate_real_kernel");
        CHECK_CUDA_SYNC();
    }
};

template struct exx_conjugate_real_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct exx_conjugate_real_op<std::complex<double>, base_device::DEVICE_GPU>;

template <typename FPTYPE>
__global__ void exx_gather_recip_kernel(const thrust::complex<FPTYPE>* in,
                                        thrust::complex<FPTYPE>* out,
                                        const int* map,
                                        int nout)
{
    const int ig = blockIdx.x * blockDim.x + threadIdx.x;
    if (ig < nout)
    {
        const int src = map[ig];
        out[ig] = src >= 0 ? in[src] : thrust::complex<FPTYPE>(0, 0);
    }
}

template <typename FPTYPE>
__global__ void exx_scatter_add_recip_kernel(const thrust::complex<FPTYPE>* in,
                                             thrust::complex<FPTYPE>* out,
                                             const int* map,
                                             int nin,
                                             thrust::complex<FPTYPE> factor)
{
    const int ig = blockIdx.x * blockDim.x + threadIdx.x;
    if (ig < nin)
    {
        const int dst = map[ig];
        if (dst >= 0)
        {
            out[dst] += factor * in[ig];
        }
    }
}

template <typename FPTYPE>
struct exx_gather_recip_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>
{
    using T = std::complex<FPTYPE>;
    void operator()(const T* in, T* out, const int* map, int nout)
    {
        const int threads_per_block = 256;
        const int num_blocks = (nout + threads_per_block - 1) / threads_per_block;
        exx_gather_recip_kernel<FPTYPE><<<num_blocks, threads_per_block>>>(
            reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
            reinterpret_cast<thrust::complex<FPTYPE>*>(out),
            map,
            nout);
        CHECK_LAST_CUDA_ERROR("exx_gather_recip_kernel");
        CHECK_CUDA_SYNC();
    }
};

template <typename FPTYPE>
struct exx_scatter_add_recip_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>
{
    using T = std::complex<FPTYPE>;
    void operator()(const T* in, T* out, const int* map, int nin, T factor)
    {
        const int threads_per_block = 256;
        const int num_blocks = (nin + threads_per_block - 1) / threads_per_block;
        exx_scatter_add_recip_kernel<FPTYPE><<<num_blocks, threads_per_block>>>(
            reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
            reinterpret_cast<thrust::complex<FPTYPE>*>(out),
            map,
            nin,
            thrust::complex<FPTYPE>(static_cast<FPTYPE>(factor.real()), static_cast<FPTYPE>(factor.imag())));
        CHECK_LAST_CUDA_ERROR("exx_scatter_add_recip_kernel");
        CHECK_CUDA_SYNC();
    }
};

template struct exx_gather_recip_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct exx_gather_recip_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct exx_scatter_add_recip_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct exx_scatter_add_recip_op<std::complex<double>, base_device::DEVICE_GPU>;
} // namespace hamilt
