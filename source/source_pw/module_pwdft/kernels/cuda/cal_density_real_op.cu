#include "source_pw/module_pwdft/kernels/cal_density_real_op.h"
#include "source_psi/psi.h"

#include <thrust/complex.h>

namespace hamilt
{
template <typename FPTYPE>
__global__ void cal_density_real_kernel(
    const thrust::complex<FPTYPE> *in1,
    const thrust::complex<FPTYPE> *in2,
    thrust::complex<FPTYPE> *out,
    const FPTYPE one_over_omega,
    int nrxx)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = idx; i < nrxx ; i+=stride)
    {
        out[i] = in1[i] * thrust::conj(in2[i]) * static_cast<thrust::complex<FPTYPE>>(one_over_omega);
    }
}

// Batched version: process multiple density calculations in one launch
template <typename FPTYPE>
__global__ void cal_density_real_kernel_batch(
    const thrust::complex<FPTYPE> *psi_nk,          // Constant input (nrxx elements)
    const thrust::complex<FPTYPE> *psi_mq_batch,    // Batch input (batch_size × nrxx)
    thrust::complex<FPTYPE> *density_batch,         // Batch output (batch_size × nrxx)
    const FPTYPE one_over_omega,
    int nrxx,
    int batch_size)
{
    // 1D thread indexing across both spatial and batch dimensions
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    int ntotal = nrxx * batch_size;

    for (int i = idx; i < ntotal; i+=stride)
    {
        // Calculate: density[batch,spatial] = psi_nk[spatial] * conj(psi_mq[batch,spatial]) / omega
        density_batch[i] = psi_nk[i % nrxx]
                                  * thrust::conj(psi_mq_batch[i])
                                  * static_cast<thrust::complex<FPTYPE>>(one_over_omega);
    }
}

template <typename FPTYPE>
struct cal_density_real_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>
{
    using T = std::complex<FPTYPE>;
    void operator()(const T *psi1, const T *psi2, T *out, double omega, int nrxx)
    {
        int threads_per_block = 256;
        int num_blocks = (nrxx + threads_per_block - 1) / threads_per_block;
        double one_over_omega = 1 / omega;

        cal_density_real_kernel<FPTYPE><<<num_blocks, threads_per_block>>>(
            reinterpret_cast<const thrust::complex<FPTYPE> *>(psi1),
            reinterpret_cast<const thrust::complex<FPTYPE> *>(psi2),
            reinterpret_cast<thrust::complex<FPTYPE> *>(out),
            static_cast<FPTYPE>(one_over_omega), nrxx);

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            throw std::runtime_error("CUDA error in cal_density_real_kernel: " + std::string(cudaGetErrorString(err)));
        }
    }

    // Batched operator
    void operator_batch(
        const T *psi_nk,
        const T *psi_mq_batch,
        T *density_batch,
        double omega,
        int nrxx,
        int batch_size)
    {
        const int total_elements = nrxx * batch_size;
        const int threads_per_block = 256;
        const int blocks = (total_elements + threads_per_block - 1) / threads_per_block;

        double one_over_omega = 1 / omega;
        cal_density_real_kernel_batch<FPTYPE><<<blocks, threads_per_block>>>(
            reinterpret_cast<const thrust::complex<FPTYPE>*>(psi_nk),
            reinterpret_cast<const thrust::complex<FPTYPE>*>(psi_mq_batch),
            reinterpret_cast<thrust::complex<FPTYPE>*>(density_batch),
            static_cast<FPTYPE>(one_over_omega),
            nrxx,
            batch_size);

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            throw std::runtime_error("CUDA error in cal_density_real_kernel_batch: " + std::string(cudaGetErrorString(err)));
        }
    }
};

template struct cal_density_real_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct cal_density_real_op<std::complex<double>, base_device::DEVICE_GPU>;
} // namespace hamilt