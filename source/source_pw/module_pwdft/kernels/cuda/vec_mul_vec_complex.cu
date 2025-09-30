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

};
template struct vec_mul_vec_complex_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct vec_mul_vec_complex_op<std::complex<double>, base_device::DEVICE_GPU>;
} // hamilt