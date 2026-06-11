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
} // namespace hamilt
