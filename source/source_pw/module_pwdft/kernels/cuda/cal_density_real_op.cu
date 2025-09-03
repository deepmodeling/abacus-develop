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
    const FPTYPE omega,
    int nrxx)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < nrxx)
    {
        out[idx] = in1[idx] * thrust::conj(in2[idx]) / static_cast<thrust::complex<FPTYPE>>(omega);
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

        cal_density_real_kernel<FPTYPE><<<num_blocks, threads_per_block>>>(
            reinterpret_cast<const thrust::complex<FPTYPE> *>(psi1),
            reinterpret_cast<const thrust::complex<FPTYPE> *>(psi2),
            reinterpret_cast<thrust::complex<FPTYPE> *>(out),
            static_cast<FPTYPE>(omega), nrxx);

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            throw std::runtime_error("CUDA error in cal_density_real_kernel: " + std::string(cudaGetErrorString(err)));
        }
    }
};

template struct cal_density_real_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct cal_density_real_op<std::complex<double>, base_device::DEVICE_GPU>;
} // namespace hamilt