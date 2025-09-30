#include "source_base/kernels/math_kernel_op.h"
#include "source_hsolver/kernels/bpcg_kernel_op.h"

#include <base/macros/macros.h>
#include <thrust/complex.h>
namespace hsolver
{
const int warp_size = 32;
const int thread_per_block = 256;
#define FULL_MASK 0xffffffff
#define WARP_SIZE 32

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

template <typename Real>
__global__ void apply_eigenvalues_kernel(
        const thrust::complex<Real>* vectors,
        thrust::complex<Real>* result,
        const Real* eigenvalues,
        const int nbase,
        const int nbase_x,
        const int notconv)
{
    int m = blockIdx.x;
    int idx = threadIdx.x + blockIdx.y * blockDim.x;

    if (m < notconv && idx < nbase) {
        result[m * nbase_x + idx] = eigenvalues[m] * vectors[m * nbase_x + idx];
    }
}

template <typename Real>
__global__ void precondition_kernel(
        thrust::complex<Real>* psi_iter,
        const Real* precondition,
        const Real* eigenvalues,
        const int dim,
        const int nbase,
        const int notconv)
{
    int m = blockIdx.x;
    int i = threadIdx.x + blockIdx.y * blockDim.x;

    if (m < notconv && i < dim) {
        Real x = abs(precondition[i] - eigenvalues[m]);
        Real pre = 0.5 * (1.0 + x + sqrt(1 + (x - 1.0) * (x - 1.0)));
        psi_iter[(nbase + m) * dim + i] = psi_iter[(nbase + m) * dim + i] / pre;
    }
}

template <typename Real>
__device__ Real warpReduceSum(Real val) {
    for (int offset = WARP_SIZE / 2; offset > 0; offset >>= 1)
        val += __shfl_down_sync(FULL_MASK, val, offset);
    return val;
}

template <typename Real>
__device__ Real blockReduceSum(Real val, volatile Real* shared) {
    int lane = threadIdx.x % WARP_SIZE;
    int wid  = threadIdx.x / WARP_SIZE;

    val = warpReduceSum(val);

    if (lane == 0)
        shared[wid] = val;

    __syncthreads();

    Real sum = 0.0;
    if (wid == 0) {
        sum = (threadIdx.x < blockDim.x / 32) ? shared[lane] : 0.0;
        sum = warpReduceSum(sum);
        if (lane == 0) shared[0] = sum;
    }

    __syncthreads();
    return shared[0];
}


template <typename Real>
__global__ void normalize_kernel(
        thrust::complex<Real>* psi_iter,
        Real* psi_norm,
        const int dim,
        const int nbase,
        const int notconv)
{
    int m = blockIdx.x;
    int tid = threadIdx.x;
    extern __shared__ char s_char[];
    Real* shared = reinterpret_cast<Real*>(s_char);

    Real local_sum = 0.0;

    // Calculate the sum for normalization
    for (int i = tid; i < dim; i += thread_per_block) {
        auto val = psi_iter[(nbase + m) * dim + i];
        local_sum += (val * thrust::conj(val)).real();
    }

    Real l2_sq = blockReduceSum(local_sum, shared);
    Real norm = sqrt(l2_sq);

    // Normalize the vector
    for (int i = tid; i < dim; i += thread_per_block) {
        psi_iter[(nbase + m) * dim + i] /= norm;
    }

    // Store the norm if needed
    if (tid == 0 && psi_norm != nullptr) {
        psi_norm[m] = norm;
    }
}

template <typename T, typename Real>
__global__ void refresh_hcc_scc_vcc_kernel(
        const int n,
        T *hcc,
        T *scc,
        T *vcc,
        const int ldh,
        const Real *eigenvalue,
        const T one)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n)
    {
        hcc[i * ldh + i] = eigenvalue[i];
        scc[i * ldh + i] = one;
        vcc[i * ldh + i] = one;
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

template <typename T>
void apply_eigenvalues_op<T, base_device::DEVICE_GPU>::operator()(const int& nbase,
                                                                const int& nbase_x,
                                                                const int& notconv,
                                                                T* result,
                                                                const T* vectors,
                                                                const Real* eigenvalues)
{
    const int threads_per_block = 256;
    const int blocks_per_grid_y = (nbase + threads_per_block - 1) / threads_per_block;

    dim3 grid(notconv, blocks_per_grid_y);

    auto vec_complex = reinterpret_cast<const thrust::complex<Real>*>(vectors);
    auto res_complex = reinterpret_cast<thrust::complex<Real>*>(result);

    apply_eigenvalues_kernel<Real><<<grid, threads_per_block>>>(
        vec_complex, res_complex, eigenvalues, nbase, nbase_x, notconv);

    cudaCheckOnDebug();
}

template <typename T>
void precondition_op<T, base_device::DEVICE_GPU>::operator()(const int& dim,
                                                           T* psi_iter,
                                                           const int& nbase,
                                                           const int& notconv,
                                                           const Real* precondition,
                                                           const Real* eigenvalues)
{
    const int threads_per_block = 256;
    const int blocks_per_grid_y = (dim + threads_per_block - 1) / threads_per_block;

    dim3 grid(notconv, blocks_per_grid_y);

    auto psi_complex = reinterpret_cast<thrust::complex<Real>*>(psi_iter);

    precondition_kernel<Real><<<grid, threads_per_block>>>(
        psi_complex, precondition, eigenvalues, dim, nbase, notconv);

    cudaCheckOnDebug();
}

template <typename T>
void normalize_op<T, base_device::DEVICE_GPU>::operator()(const int& dim,
                                                        T* psi_iter,
                                                        const int& nbase,
                                                        const int& notconv,
                                                        Real* psi_norm)
{
    auto psi_complex = reinterpret_cast<thrust::complex<Real>*>(psi_iter);
    int sharedMemSize = (thread_per_block / WARP_SIZE) * sizeof(Real);

    normalize_kernel<Real><<<notconv, thread_per_block, sharedMemSize, 0>>>(
        psi_complex, psi_norm, dim, nbase, notconv);

    cudaCheckOnDebug();
}

template <>
void refresh_hcc_scc_vcc_op<double, base_device::DEVICE_GPU>::operator()(const int &n,
                  double *hcc,
                  double *scc,
                  double *vcc,
                  const int &ldh,
                  const double *eigenvalue,
                  const double& one)
{
    int thread = 512;
    int block = (n + thread - 1) / thread;
    refresh_hcc_scc_vcc_kernel<double, double> <<<block, thread >>> (n, hcc, scc, vcc, ldh, eigenvalue, one);

    cudaCheckOnDebug();
}

template <>
void refresh_hcc_scc_vcc_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const int &n,
                  std::complex<float> *hcc,
                  std::complex<float> *scc,
                  std::complex<float> *vcc,
                  const int &ldh,
                  const float *eigenvalue,
                  const std::complex<float>& one)
{
    int thread = 512;
    int block = (n + thread - 1) / thread;
    refresh_hcc_scc_vcc_kernel<thrust::complex<float>, float> <<<block, thread >>> (n, reinterpret_cast<thrust::complex<float>*>(hcc),
                    reinterpret_cast<thrust::complex<float>*>(scc), reinterpret_cast<thrust::complex<float>*>(vcc), ldh, eigenvalue,
                    thrust::complex<float>(one));

    cudaCheckOnDebug();
}

template <>
void refresh_hcc_scc_vcc_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(const int &n,
                  std::complex<double> *hcc,
                  std::complex<double> *scc,
                  std::complex<double> *vcc,
                  const int &ldh,
                  const double *eigenvalue,
                  const std::complex<double>& one)
{
    int thread = 512;
    int block = (n + thread - 1) / thread;
    refresh_hcc_scc_vcc_kernel<thrust::complex<double>, double> <<<block, thread >>> (n, reinterpret_cast<thrust::complex<double>*>(hcc),
                    reinterpret_cast<thrust::complex<double>*>(scc), reinterpret_cast<thrust::complex<double>*>(vcc), ldh, eigenvalue,
                    thrust::complex<double>(one));

    cudaCheckOnDebug();
}

template struct calc_grad_with_block_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct line_minimize_with_block_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct calc_grad_with_block_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct line_minimize_with_block_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct apply_eigenvalues_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct apply_eigenvalues_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct apply_eigenvalues_op<double, base_device::DEVICE_GPU>;
template struct precondition_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct precondition_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct precondition_op<double, base_device::DEVICE_GPU>;
template struct normalize_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct normalize_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct normalize_op<double, base_device::DEVICE_GPU>;
template struct refresh_hcc_scc_vcc_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct refresh_hcc_scc_vcc_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct refresh_hcc_scc_vcc_op<double, base_device::DEVICE_GPU>;
}