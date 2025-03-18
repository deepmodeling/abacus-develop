#include "module_base/kernels/math_kernel_op.h"
#include "module_hsolver/kernels/bpcg_kernel_op.h"

#include <base/macros/macros.h>
#include <thrust/complex.h>
#define WARP_SIZE 32
#define THREAD_PER_BLOCK 256
namespace hsolver
{
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
    __shared__ Real data[THREAD_PER_BLOCK * 3];

    data[tid] = 0;

    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += THREAD_PER_BLOCK) {
        item = band_idx * n_basis_max + basis_idx;
        data[tid] += (grad[item] * thrust::conj(grad[item])).real();
    }
    __syncthreads();
    // just do some parallel reduction in shared memory
    for (int ii = THREAD_PER_BLOCK >> 1; ii > 0; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
        }
        __syncthreads();
    }

    Real norm = 1.0 / sqrt(data[0]);
    __syncthreads();

    data[tid] = 0;
    data[THREAD_PER_BLOCK + tid] = 0;
    data[2 * THREAD_PER_BLOCK + tid] = 0;
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += THREAD_PER_BLOCK) {
        item = band_idx * n_basis_max + basis_idx;
        grad[item] *= norm;
        hgrad[item] *= norm;
        data[tid] += (hpsi[item] * thrust::conj(psi[item])).real();
        data[THREAD_PER_BLOCK + tid] += (grad[item] * thrust::conj(hpsi[item])).real();
        data[2 * THREAD_PER_BLOCK + tid] += (grad[item] * thrust::conj(hgrad[item])).real();
    }
    __syncthreads();

    // just do some parallel reduction in shared memory
    for (int ii = THREAD_PER_BLOCK >> 1; ii > 0; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
            data[THREAD_PER_BLOCK + tid] += data[THREAD_PER_BLOCK + tid + ii];
            data[2 * THREAD_PER_BLOCK + tid] += data[2 * THREAD_PER_BLOCK + tid + ii];
        }
        __syncthreads();
    }
    epsilo_0 = data[0];
    epsilo_1 = data[THREAD_PER_BLOCK];
    epsilo_2 = data[2 * THREAD_PER_BLOCK];

    theta = 0.5 * abs(atan(2 * epsilo_1/(epsilo_0 - epsilo_2)));
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += THREAD_PER_BLOCK) {
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
    __shared__ Real data[THREAD_PER_BLOCK * 2];

    // Init shared memory
    data[tid] = 0;

    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += THREAD_PER_BLOCK) {
        item = band_idx * n_basis_max + basis_idx;
        data[tid] += (psi[item] * thrust::conj(psi[item])).real();
    }
    __syncthreads();
    // just do some parallel reduction in shared memory
    for (int ii = THREAD_PER_BLOCK >> 1; ii > 0; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
        }
        __syncthreads();
    }

    Real norm = 1.0 / sqrt(data[0]);
    __syncthreads();

    data[tid] = 0;
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += THREAD_PER_BLOCK) {
        item = band_idx * n_basis_max + basis_idx;
        psi[item] *= norm;
        hpsi[item] *= norm;
        data[tid] += (hpsi[item] * thrust::conj(psi[item])).real();
    }
    __syncthreads();

    // just do some parallel reduction in shared memory
    for (int ii = THREAD_PER_BLOCK >> 1; ii > 0; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
        }
        __syncthreads();
    }
    epsilo = data[0];
    __syncthreads();

    data[tid] = 0;
    data[THREAD_PER_BLOCK + tid] = 0;
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += THREAD_PER_BLOCK) {
        item = band_idx * n_basis_max + basis_idx;
        grad_1 = hpsi[item] - epsilo * psi[item];
        grad_2 = thrust::norm(grad_1);
        data[tid] += grad_2;
        data[THREAD_PER_BLOCK + tid] += grad_2 / prec[basis_idx];
    }
    __syncthreads();

    // just do some parallel reduction in shared memory
    for (int ii = THREAD_PER_BLOCK >> 1; ii > 0; ii >>= 1) {
        if (tid < ii) {
            data[tid] += data[tid + ii];
            data[THREAD_PER_BLOCK + tid] += data[THREAD_PER_BLOCK + tid + ii];
        }
        __syncthreads();
    }
    err_st = data[0];
    beta_st = data[THREAD_PER_BLOCK];
    for (int basis_idx = tid; basis_idx < n_basis; basis_idx += THREAD_PER_BLOCK) {
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

    line_minimize_with_block<Real><<<n_band, THREAD_PER_BLOCK>>>(
            A, B, C, D,
            n_basis, n_basis_max);

    hipCheckOnDebug();
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

    calc_grad_with_block<Real><<<n_band, THREAD_PER_BLOCK>>>(
            prec_in, err_out, beta_out,
            A, B, C, D,
            n_basis, n_basis_max);

    hipCheckOnDebug();
}

template struct calc_grad_with_block_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct line_minimize_with_block_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct calc_grad_with_block_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct line_minimize_with_block_op<std::complex<double>, base_device::DEVICE_GPU>;
}