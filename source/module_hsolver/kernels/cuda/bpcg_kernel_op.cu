#include "module_base/kernels/math_kernel_op.h"
#include "module_hsolver/kernels/bpcg_kernel_op.h"

#include <base/macros/macros.h>
#include <thrust/complex.h>
namespace hsolver
{
const int warp_size = 32;
const int thread_per_block = 256;

template <typename Real>
__global__ void line_minimize_with_block(
        thrust::complex<Real>* grad,
        thrust::complex<Real>* hgrad,
        thrust::complex<Real>* psi,
        thrust::complex<Real>* hpsi,
        const int n_basis,
        const int n_basis_max)
{
    int band_idx = blockIdx.x; // band_idx
    int tid = threadIdx.x; // basis_idx
    int item = 0;
    Real epsilo_0 = 0.0, epsilo_1 = 0.0, epsilo_2 = 0.0;
    Real theta = 0.0, cos_theta = 0.0, sin_theta = 0.0;
    __shared__ Real data[thread_per_block * 3];

    data[tid] = 0;

    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += thread_per_block) {
        item = band_idx * n_basis_max + basis_idx;
        data[tid] += (grad[item] * thrust::conj(grad[item])).real();
    }
    __syncthreads();
    // just do some parallel reduction in shared memory
    for (int ii = thread_per_block >> 1; ii > warp_size; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
        }
        __syncthreads();
    }

    // For threads in the same warp, it is better that they process the same work
    // Also, __syncwarp() should be used instead of __syncthreads()
    // Therefore we unroll the loop and ensure that the threads does the same work
    if (tid < warp_size) {
        data[tid] += data[tid + 32]; __syncwarp();
        data[tid] += data[tid + 16]; __syncwarp();
        data[tid] += data[tid + 8]; __syncwarp();
        data[tid] += data[tid + 4]; __syncwarp();
        data[tid] += data[tid + 2]; __syncwarp();
        data[tid] += data[tid + 1]; __syncwarp();
    }

    __syncthreads();

    Real norm = 1.0 / sqrt(data[0]);
    __syncthreads();

    data[tid] = 0;
    data[thread_per_block + tid] = 0;
    data[2 * thread_per_block + tid] = 0;
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += thread_per_block) {
        item = band_idx * n_basis_max + basis_idx;
        grad[item] *= norm;
        hgrad[item] *= norm;
        data[tid] += (hpsi[item] * thrust::conj(psi[item])).real();
        data[thread_per_block + tid] += (grad[item] * thrust::conj(hpsi[item])).real();
        data[2 * thread_per_block + tid] += (grad[item] * thrust::conj(hgrad[item])).real();
    }
    __syncthreads();

    // just do some parallel reduction in shared memory
    for (int ii = thread_per_block >> 1; ii > warp_size; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
            data[thread_per_block + tid] += data[thread_per_block + tid + ii];
            data[2 * thread_per_block + tid] += data[2 * thread_per_block + tid + ii];
        }
        __syncthreads();
    }
    if (tid < warp_size) {
        data[tid] += data[tid + 32]; __syncwarp();
        data[tid] += data[tid + 16]; __syncwarp();
        data[tid] += data[tid + 8]; __syncwarp();
        data[tid] += data[tid + 4]; __syncwarp();
        data[tid] += data[tid + 2]; __syncwarp();
        data[tid] += data[tid + 1]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 32]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 16]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 8]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 4]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 2]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 1]; __syncwarp();
        data[2 * thread_per_block + tid] += data[2 * thread_per_block + tid + 32]; __syncwarp();
        data[2 * thread_per_block + tid] += data[2 * thread_per_block + tid + 16]; __syncwarp();
        data[2 * thread_per_block + tid] += data[2 * thread_per_block + tid + 8]; __syncwarp();
        data[2 * thread_per_block + tid] += data[2 * thread_per_block + tid + 4]; __syncwarp();
        data[2 * thread_per_block + tid] += data[2 * thread_per_block + tid + 2]; __syncwarp();
        data[2 * thread_per_block + tid] += data[2 * thread_per_block + tid + 1]; __syncwarp();
    }
    __syncthreads();
    epsilo_0 = data[0];
    epsilo_1 = data[thread_per_block];
    epsilo_2 = data[2 * thread_per_block];

    theta = 0.5 * abs(atan(2 * epsilo_1/(epsilo_0 - epsilo_2)));
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += thread_per_block) {
        item = band_idx * n_basis_max + basis_idx;
        psi [item] = psi [item] * cos_theta + grad [item] * sin_theta;
        hpsi[item] = hpsi[item] * cos_theta + hgrad[item] * sin_theta;
    }
}

template <typename Real>
__global__ void calc_grad_with_block(
        const Real* prec,
        Real* err,
        Real* beta,
        thrust::complex<Real>* psi,
        thrust::complex<Real>* hpsi,
        thrust::complex<Real>* grad,
        thrust::complex<Real>* grad_old,
        const int n_basis,
        const int n_basis_max)
{
    int band_idx = blockIdx.x; // band_idx
    int tid = threadIdx.x; // basis_idx
    int item = 0;
    Real err_st = 0.0;
    Real beta_st = 0.0;
    Real epsilo = 0.0;
    Real grad_2 = 0.0;
    thrust::complex<Real> grad_1 = {0, 0};
    __shared__ Real data[thread_per_block * 2];

    // Init shared memory
    data[tid] = 0;

    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += thread_per_block) {
        item = band_idx * n_basis_max + basis_idx;
        data[tid] += (psi[item] * thrust::conj(psi[item])).real();
    }
    __syncthreads();
    // just do some parallel reduction in shared memory
    for (int ii = thread_per_block >> 1; ii > warp_size; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
        }
        __syncthreads();
    }

    if (tid < warp_size) {
        data[tid] += data[tid + 32]; __syncwarp();
        data[tid] += data[tid + 16]; __syncwarp();
        data[tid] += data[tid + 8]; __syncwarp();
        data[tid] += data[tid + 4]; __syncwarp();
        data[tid] += data[tid + 2]; __syncwarp();
        data[tid] += data[tid + 1]; __syncwarp();
    }

    __syncthreads();

    Real norm = 1.0 / sqrt(data[0]);
    __syncthreads();

    data[tid] = 0;
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += thread_per_block) {
        item = band_idx * n_basis_max + basis_idx;
        psi[item] *= norm;
        hpsi[item] *= norm;
        data[tid] += (hpsi[item] * thrust::conj(psi[item])).real();
    }
    __syncthreads();

    // just do some parallel reduction in shared memory
    for (int ii = thread_per_block >> 1; ii > warp_size; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
        }
        __syncthreads();
    }

    if (tid < warp_size) {
        data[tid] += data[tid + 32]; __syncwarp();
        data[tid] += data[tid + 16]; __syncwarp();
        data[tid] += data[tid + 8]; __syncwarp();
        data[tid] += data[tid + 4]; __syncwarp();
        data[tid] += data[tid + 2]; __syncwarp();
        data[tid] += data[tid + 1]; __syncwarp();
    }

    __syncthreads();
    epsilo = data[0];
    __syncthreads();

    data[tid] = 0;
    data[thread_per_block + tid] = 0;
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += thread_per_block) {
        item = band_idx * n_basis_max + basis_idx;
        grad_1 = hpsi[item] - epsilo * psi[item];
        grad_2 = thrust::norm(grad_1);
        data[tid] += grad_2;
        data[thread_per_block + tid] += grad_2 / prec[basis_idx];
    }
    __syncthreads();

    // just do some parallel reduction in shared memory
    for (int ii = thread_per_block >> 1; ii > warp_size; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
            data[thread_per_block + tid] += data[thread_per_block + tid + ii];
        }
        __syncthreads();
    }

    if (tid < warp_size) {
        data[tid] += data[tid + 32]; __syncwarp();
        data[tid] += data[tid + 16]; __syncwarp();
        data[tid] += data[tid + 8]; __syncwarp();
        data[tid] += data[tid + 4]; __syncwarp();
        data[tid] += data[tid + 2]; __syncwarp();
        data[tid] += data[tid + 1]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 32]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 16]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 8]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 4]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 2]; __syncwarp();
        data[thread_per_block + tid] += data[thread_per_block + tid + 1]; __syncwarp();
    }

    __syncthreads();
    err_st = data[0];
    beta_st = data[thread_per_block];
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += thread_per_block) {
        item = band_idx * n_basis_max + basis_idx;
        grad_1 = hpsi[item] - epsilo * psi[item];
        grad[item] = -grad_1 / prec[basis_idx] + beta_st / beta[band_idx] * grad_old[item];
    }

    __syncthreads();
    if (tid == 0) {
        beta[band_idx] = beta_st;
        err[band_idx] = sqrt(err_st);
    }
}

template <typename T>
void line_minimize_with_block_op<T, base_device::DEVICE_GPU>::operator()(T* grad_out,
                                                                         T* hgrad_out,
                                                                         T* psi_out,
                                                                         T* hpsi_out,
                                                                         const int& n_basis,
                                                                         const int& n_basis_max,
                                                                         const int& n_band)
{
    auto A = reinterpret_cast<thrust::complex<Real>*>(grad_out);
    auto B = reinterpret_cast<thrust::complex<Real>*>(hgrad_out);
    auto C = reinterpret_cast<thrust::complex<Real>*>(psi_out);
    auto D = reinterpret_cast<thrust::complex<Real>*>(hpsi_out);

    line_minimize_with_block<Real><<<n_band, thread_per_block>>>(
            A, B, C, D,
            n_basis, n_basis_max);

    cudaCheckOnDebug();
}

template <typename T>
void calc_grad_with_block_op<T, base_device::DEVICE_GPU>::operator()(const Real* prec_in,
                                                                     Real* err_out,
                                                                     Real* beta_out,
                                                                     T* psi_out,
                                                                     T* hpsi_out,
                                                                     T* grad_out,
                                                                     T* grad_old_out,
                                                                     const int& n_basis,
                                                                     const int& n_basis_max,
                                                                     const int& n_band)
{
    auto A = reinterpret_cast<thrust::complex<Real>*>(psi_out);
    auto B = reinterpret_cast<thrust::complex<Real>*>(hpsi_out);
    auto C = reinterpret_cast<thrust::complex<Real>*>(grad_out);
    auto D = reinterpret_cast<thrust::complex<Real>*>(grad_old_out);

    calc_grad_with_block<Real><<<n_band, thread_per_block>>>(
            prec_in, err_out, beta_out,
            A, B, C, D,
            n_basis, n_basis_max);

    cudaCheckOnDebug();
}

template struct calc_grad_with_block_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct line_minimize_with_block_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct calc_grad_with_block_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct line_minimize_with_block_op<std::complex<double>, base_device::DEVICE_GPU>;

}