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

template <typename FPTYPE>
struct mul_potential_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>
{
    using T = std::complex<FPTYPE>;
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
            pot + ik * nks * npw + iq * npw,
            reinterpret_cast<thrust::complex<FPTYPE>*>(density_recip),
            npw);

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            throw std::runtime_error("CUDA error in mul_potential_kernel: " + std::string(cudaGetErrorString(err)));
        }
    }
};
template struct mul_potential_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct mul_potential_op<std::complex<double>, base_device::DEVICE_GPU>;
} // hamilt