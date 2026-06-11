#include "source_pw/module_pwdft/kernels/mul_potential_op.h"
#include "source_io/module_parameter/parameter.h"
#include "source_psi/psi.h"

#include <thrust/complex.h>

namespace hamilt {
template <typename FPTYPE>
__global__ void mul_potential_kernel(
    const FPTYPE *pot_shifted,
    thrust::complex<FPTYPE> *density_recip,
    int npw)
{
    int ig = blockIdx.x * blockDim.x + threadIdx.x;
    if (ig < npw)
    {
        density_recip[ig] *= pot_shifted[ig];
    }
}

// Batched kernel: process multiple density arrays with same potential
template <typename FPTYPE>
__global__ void mul_potential_kernel_batch(
    const FPTYPE *pot,                            // Potential (npw elements, shared)
    thrust::complex<FPTYPE> *density_recip_batch, // Batch arrays (batch_size × npw)
    int npw,
    int batch_size)
{
    // 1D thread indexing across both batch and spatial dimensions
    int linear_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total_elements = npw * batch_size;

    if (linear_idx < total_elements)
    {
        int batch_idx = linear_idx / npw;           // Which batch element (0 to batch_size-1)
        int ig = linear_idx % npw;                  // Which plane wave (0 to npw-1)

        // Apply same potential to all batch elements
        density_recip_batch[linear_idx] *= pot[ig];
    }
}

template <typename FPTYPE>
struct mul_potential_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>
{
    using T = std::complex<FPTYPE>;

    // Single operator
    void operator()(const FPTYPE *pot, T *density_recip, int npw, int nks, int ik, int iq)
    {
// #ifdef _OPENMP
// #pragma omp parallel for schedule(static)
// #endif
//         for (int ig = 0; ig < npw; ig++)
//         {
//             int ig_kq = ik * nks * npw + iq * npw + ig;
//             density_recip[ig] *= pot[ig_kq];
//
//         }
        int threads_per_block = 256;
        int num_blocks = (npw + threads_per_block - 1) / threads_per_block;

        mul_potential_kernel<<<num_blocks, threads_per_block>>>(
            pot,
            reinterpret_cast<thrust::complex<FPTYPE>*>(density_recip),
            npw);

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            throw std::runtime_error("CUDA error in mul_potential_kernel: " + std::string(cudaGetErrorString(err)));
        }
    }

    // Batched operator
    void operator_batch(
        const FPTYPE *pot,
        T *density_recip_batch,
        int npw,
        int batch_size)
    {
        const int total_elements = npw * batch_size;
        const int threads_per_block = 256;
        const int num_blocks = (total_elements + threads_per_block - 1) / threads_per_block;

        mul_potential_kernel_batch<FPTYPE><<<num_blocks, threads_per_block>>>(
            pot,
            reinterpret_cast<thrust::complex<FPTYPE>*>(density_recip_batch),
            npw,
            batch_size);

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            throw std::runtime_error("CUDA error in mul_potential_kernel_batch: "
                                   + std::string(cudaGetErrorString(err)));
        }
    }
};
template struct mul_potential_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct mul_potential_op<std::complex<double>, base_device::DEVICE_GPU>;
} // hamilt