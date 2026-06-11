#include "source_pw/module_pwdft/kernels/exx_cal_energy_op.h"
#include "source_psi/psi.h"
#include "source_base/kernels/math_kernel_op.h"

#include <thrust/complex.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include "source_base/module_device/device.h"
#include "source_base/module_device/kernel_compat.h"
#include "source_base/module_device/device_check.h"

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
__global__ void cal_vec_norm_accumulate_kernel(
    const thrust::complex<FPTYPE> *den,
    const FPTYPE *pot,
    const FPTYPE scalar,
    FPTYPE *result,
    int npw)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < npw)
    {
        FPTYPE tmp = (den[idx] * thrust::conj(den[idx])).real() * pot[idx] * scalar;
        atomicAdd(result, tmp);
    }
}

// Kernel to compute element-wise norm of a complex vector and save to real vector
template <typename FPTYPE>
__global__ void cal_vec_elem_norm_squared(
    const thrust::complex<FPTYPE> *vector_in,
    FPTYPE *vector_out,
    int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = idx; i < n; i += stride)
    {
        vector_out[i] = (vector_in[i] * thrust::conj(vector_in[i])).real();
    }
}


template <typename FPTYPE>
struct exx_density_potential_mul_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>
{
    using T = std::complex<FPTYPE>;
    FPTYPE operator()(const T *vector_in,
        FPTYPE *vector_buffer, const FPTYPE *pot,
        FPTYPE *vec_temp, FPTYPE *weights,
        int npw, int batch_idx)
    {
        const FPTYPE one{1}; // Scalar one for gemv
        const FPTYPE zero{0}; // Scalar one for gemv
        const int inc{1}; // Increment for gemv
        int threads_per_block = 256;
        int num_blocks = (npw * batch_idx + threads_per_block - 1) / threads_per_block;
        cal_vec_elem_norm_squared<FPTYPE><<<num_blocks, threads_per_block>>>(
            reinterpret_cast<const thrust::complex<FPTYPE> *>(vector_in),
            vector_buffer,
            npw * batch_idx
        );
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            throw std::runtime_error("CUDA error in cal_vec_norm_real_kernel: " + std::string(cudaGetErrorString(err)));
        }
        ModuleBase::gemv_op<FPTYPE, base_device::DEVICE_GPU>()(
            'T',
            npw,          // m
            batch_idx,           // n
            &one,            // alpha
            vector_buffer,   // matrix as (npw, batch_idx ) in column major layour
            npw,         // lda
            pot,  // vector (npw, )
            inc,                   // incx
            &zero,            // beta
            vec_temp,          // batch_idx vector (batch_idx)
            inc);                  // incy
        // Obtain the energy
        return ModuleBase::dot_real_op<FPTYPE, base_device::DEVICE_GPU>()(
            batch_idx,
            vec_temp,
            weights,
            false
        );
    }
};

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
        CHECK_CUDA(cudaMalloc(&d_result, sizeof(FPTYPE)));
        CHECK_CUDA(cudaMemset(d_result, 0, sizeof(FPTYPE)));

        cal_vec_norm_kernel<FPTYPE><<<num_blocks, threads_per_block>>>(
            reinterpret_cast<const thrust::complex<FPTYPE> *>(den),
            pot,
            d_result,
            npw);

        CHECK_LAST_CUDA_ERROR("cal_vec_norm_kernel");
        CHECK_CUDA_SYNC();

        CHECK_CUDA(cudaMemcpy(&result, d_result, sizeof(FPTYPE), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaFree(d_result));

        return scalar * result;
    }
};

template struct exx_cal_energy_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct exx_cal_energy_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct exx_density_potential_mul_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct exx_density_potential_mul_op<std::complex<double>, base_device::DEVICE_GPU>;

template <typename FPTYPE>
struct exx_cal_energy_accumulate_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>
{
    using T = std::complex<FPTYPE>;
    void operator()(const T *den, const FPTYPE *pot, FPTYPE scalar, int npw, FPTYPE* result)
    {
        int threads_per_block = 256;
        int num_blocks = (npw + threads_per_block - 1) / threads_per_block;

        cal_vec_norm_accumulate_kernel<FPTYPE><<<num_blocks, threads_per_block>>>(
            reinterpret_cast<const thrust::complex<FPTYPE> *>(den),
            pot,
            scalar,
            result,
            npw);

        CHECK_LAST_CUDA_ERROR("cal_vec_norm_accumulate_kernel");
        CHECK_CUDA_SYNC();
    }
};

template struct exx_cal_energy_accumulate_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct exx_cal_energy_accumulate_op<std::complex<double>, base_device::DEVICE_GPU>;
} // namespace hamilt
