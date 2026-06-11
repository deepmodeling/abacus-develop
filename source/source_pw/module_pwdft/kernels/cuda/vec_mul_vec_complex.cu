#include "source_pw/module_pwdft/kernels/vec_mul_vec_complex_op.h"
#include "source_io/module_parameter/parameter.h"
#include "source_psi/psi.h"

#include <thrust/complex.h>

namespace hamilt {
template <typename FPTYPE>
__global__ void vec_mul_vec_complex_kernel(
    const thrust::complex<FPTYPE> *vec1,
    const thrust::complex<FPTYPE> *vec2,
    thrust::complex<FPTYPE> *result,
    int size)
{
    int ig = blockIdx.x * blockDim.x + threadIdx.x;
    if (ig < size)
    {
        result[ig] = vec1[ig] * vec2[ig];
    }
}

// Batched kernel: process multiple vector pairs in parallel
template <typename FPTYPE>
__global__ void vec_mul_vec_complex_kernel_batch(
    const thrust::complex<FPTYPE> *vec1_batch,
    const thrust::complex<FPTYPE> *vec2_batch,
    thrust::complex<FPTYPE> *out_batch,
    int n,
    int batch_size)
{
    // 1D thread indexing across both batch and spatial dimensions
    int linear_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total_elements = n * batch_size;

    if (linear_idx < total_elements)
    {
        // Element-wise multiplication (same index in vec1 and vec2)
        out_batch[linear_idx] = vec1_batch[linear_idx] * vec2_batch[linear_idx];
    }
}

template <typename FPTYPE>
struct vec_mul_vec_complex_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>
{
    using T = std::complex<FPTYPE>;
    void operator()(const T *a, const T *b, T* out, int size)
    {
        int threads_per_block = 256;
        int num_blocks = (size + threads_per_block - 1) / threads_per_block;

        vec_mul_vec_complex_kernel<FPTYPE><<<num_blocks, threads_per_block>>>(
            reinterpret_cast<const thrust::complex<FPTYPE> *>(a),
            reinterpret_cast<const thrust::complex<FPTYPE> *>(b),
            reinterpret_cast<thrust::complex<FPTYPE> *>(out),
            size);

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            throw std::runtime_error("CUDA error in vec_mul_vec_kernel: " + std::string(cudaGetErrorString(err)));
        }
    }

    // Batched operator: process multiple vector pairs in a single kernel launch
    void operator_batch(
        const T *vec1_batch,
        const T *vec2_batch,
        T *out_batch,
        int n,
        int batch_size)
    {
        const int total_elements = n * batch_size;
        const int threads_per_block = 256;
        const int num_blocks = (total_elements + threads_per_block - 1) / threads_per_block;

        vec_mul_vec_complex_kernel_batch<FPTYPE><<<num_blocks, threads_per_block>>>(
            reinterpret_cast<const thrust::complex<FPTYPE>*>(vec1_batch),
            reinterpret_cast<const thrust::complex<FPTYPE>*>(vec2_batch),
            reinterpret_cast<thrust::complex<FPTYPE>*>(out_batch),
            n,
            batch_size);

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            throw std::runtime_error("CUDA error in vec_mul_vec_complex_kernel_batch: "
                                   + std::string(cudaGetErrorString(err)));
        }
    }

};
template struct vec_mul_vec_complex_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct vec_mul_vec_complex_op<std::complex<double>, base_device::DEVICE_GPU>;
} // hamilt