#include "source_pw/module_pwdft/kernels/exx_stress_op.h"
#include "source_base/module_device/device_check.h"

#include <stdexcept>
#include <string>
#include <thrust/complex.h>

namespace hamilt
{
template <typename FPTYPE>
__global__ void exx_stress_accumulate_kernel(const thrust::complex<FPTYPE>* density_recip,
                                             const FPTYPE* pot,
                                             const FPTYPE* pot_stress,
                                             const FPTYPE* gcar,
                                             FPTYPE dkx,
                                             FPTYPE dky,
                                             FPTYPE dkz,
                                             FPTYPE tpiba,
                                             FPTYPE scalar,
                                             int npw,
                                             FPTYPE* sigma)
{
    const int ig = blockIdx.x * blockDim.x + threadIdx.x;
    if (ig >= npw) {return;}
    const FPTYPE density_recip2 = density_recip[ig].real() * density_recip[ig].real()
                                  + density_recip[ig].imag() * density_recip[ig].imag();
    const FPTYPE kqg[3] = {(dkx + gcar[ig * 3]) * tpiba,
                           (dky + gcar[ig * 3 + 1]) * tpiba,
                           (dkz + gcar[ig * 3 + 2]) * tpiba};
    int idx = 0;
    for (int alpha = 0; alpha < 3; ++alpha)
    {
        for (int beta = alpha; beta < 3; ++beta)
        {
            const FPTYPE delta_ab = alpha == beta ? FPTYPE(1.0) : FPTYPE(0.0);
            const FPTYPE contrib = scalar * density_recip2 * pot[ig]
                                   * (kqg[alpha] * kqg[beta] * pot_stress[ig] - delta_ab);
            atomicAdd(sigma + idx, contrib);
            ++idx;
        }
    }
}

template <typename FPTYPE>
struct exx_stress_accumulate_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>
{
    using T = std::complex<FPTYPE>;
    using Real = FPTYPE;
    void operator()(const T* density_recip,
                    const Real* pot,
                    const Real* pot_stress,
                    const Real* gcar,
                    Real dkx,
                    Real dky,
                    Real dkz,
                    Real tpiba,
                    Real scalar,
                    int npw,
                    Real* sigma)
    {
        const int threads_per_block = 256;
        const int num_blocks = (npw + threads_per_block - 1) / threads_per_block;
        exx_stress_accumulate_kernel<FPTYPE><<<num_blocks, threads_per_block>>>(
            reinterpret_cast<const thrust::complex<FPTYPE>*>(density_recip),
            pot,
            pot_stress,
            gcar,
            dkx,
            dky,
            dkz,
            tpiba,
            scalar,
            npw,
            sigma);

        CHECK_LAST_CUDA_ERROR("exx_stress_accumulate_kernel");
        CHECK_CUDA_SYNC();
    }
};

template struct exx_stress_accumulate_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct exx_stress_accumulate_op<std::complex<double>, base_device::DEVICE_GPU>;
} // namespace hamilt
