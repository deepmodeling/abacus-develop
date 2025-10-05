#include "source_pw/module_pwdft/kernels/exx_cal_energy_op.h"
#include "source_psi/psi.h"

#include <thrust/complex.h>
#include "source_base/module_device/device.h"

namespace hamilt
{

// #ifdef _OPENMP
// #pragma omp parallel for reduction(+:Eexx_ik_real)
// #endif
// for (int ig = 0; ig < rhopw_dev->npw; ig++)
// {
//     int nks = wfcpw->nks;
//     int npw = rhopw_dev->npw;
//     int nk = nks / nk_fac;
//     Real Fac = pot[ik * nks * npw + iq * npw + ig];

// Eexx_ik_real += Fac * (density_recip[ig] * std::conj(density_recip[ig])).real()
//                 * wg_iqb_real / nqs * wg_ikb_real / kv->wk[ik];
// }

template <typename FPTYPE>
__global__ void cal_vec_norm_kernel(
    const thrust::complex<FPTYPE> *den,
    const FPTYPE *pot,
    FPTYPE *result,
    int npw)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < npw)
    {
        // atomicAdd(result, (den[idx].real() * den[idx].real() + den[idx].imag() * den[idx].imag()) * pot[idx]);
        FPTYPE tmp =(den[idx] * thrust::conj(den[idx])).real() * pot[idx];
        atomicAdd(result, tmp);
    }
    __syncthreads();
}

template <typename FPTYPE>
struct exx_cal_energy_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>
{
    using T = std::complex<FPTYPE>;
    FPTYPE operator()(const T *den, const FPTYPE *pot, FPTYPE scalar, int npw)
    {
        // T *den_cpu = new T[npw];
        // FPTYPE *pot_cpu = new FPTYPE[npw];
        // cudaMemcpy(den_cpu, den, npw * sizeof(T), cudaMemcpyDeviceToHost);
        // cudaMemcpy(pot_cpu, pot, npw * sizeof(FPTYPE), cudaMemcpyDeviceToHost);
        // FPTYPE result = exx_cal_energy_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>()(den_cpu, pot_cpu, scalar, npw);
        // delete[] den_cpu;
        // delete[] pot_cpu;
        // return result;
        FPTYPE result = 0.0;

        int threads_per_block = 256;
        int num_blocks = (npw + threads_per_block - 1) / threads_per_block;

        FPTYPE *d_result;
        cudaMalloc(&d_result, sizeof(FPTYPE));
        cudaMemset(d_result, 0, sizeof(FPTYPE));

        cal_vec_norm_kernel<FPTYPE><<<num_blocks, threads_per_block>>>(
            reinterpret_cast<const thrust::complex<FPTYPE> *>(den),
            pot,
            d_result,
            npw);

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            throw std::runtime_error("CUDA error in cal_vec_norm_kernel: " + std::string(cudaGetErrorString(err)));
        }

        cudaMemcpy(&result, d_result, sizeof(FPTYPE), cudaMemcpyDeviceToHost);
        cudaFree(d_result);

        return scalar * result;
    }
};

template struct exx_cal_energy_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct exx_cal_energy_op<std::complex<double>, base_device::DEVICE_GPU>;
} // namespace hamilt
