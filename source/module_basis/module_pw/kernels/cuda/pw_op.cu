#include "module_basis/module_pw/kernels/pw_op.h"

#include <thrust/complex.h>
#include <cuda_runtime.h>
#include <base/macros/macros.h>

namespace ModulePW {

#define THREADS_PER_BLOCK 256

template <typename Real>
__global__ void pack_real_grid(int size, const Real* input, thrust::complex<Real>* grid)
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < size)
    {
        grid[index] = thrust::complex<Real>(input[index], 0);
    }
}

template <typename Real>
__global__ void extract_real_grid(int size, bool add, Real factor, const thrust::complex<Real>* grid, Real* output)
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < size)
    {
        output[index] = add ? output[index] + factor * grid[index].real() : grid[index].real();
    }
}

template <typename Real>
void pack_real_grid_gpu(int size, const Real* input, std::complex<Real>* grid)
{
    pack_real_grid<<<(size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK>>>(
        size, input, reinterpret_cast<thrust::complex<Real>*>(grid));
    cudaCheckOnDebug();
}

template <typename Real>
void extract_real_grid_gpu(int size, bool add, Real factor, const std::complex<Real>* grid, Real* output)
{
    extract_real_grid<<<(size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK>>>(
        size, add, factor, reinterpret_cast<const thrust::complex<Real>*>(grid), output);
    cudaCheckOnDebug();
}

template void pack_real_grid_gpu<float>(int, const float*, std::complex<float>*);
template void pack_real_grid_gpu<double>(int, const double*, std::complex<double>*);
template void extract_real_grid_gpu<float>(int, bool, float, const std::complex<float>*, float*);
template void extract_real_grid_gpu<double>(int, bool, double, const std::complex<double>*, double*);

template <typename Real>
__global__ void gather_fft_coefficients(int size, bool add, Real scale, const int* index,
                                       const thrust::complex<Real>* grid, thrust::complex<Real>* output)
{
    const int ig = blockIdx.x * blockDim.x + threadIdx.x;
    if (ig < size)
    {
        const auto value = scale * grid[index[ig]];
        output[ig] = add ? output[ig] + value : value;
    }
}

template <typename Real>
void gather_fft_coefficients_gpu(int size, bool add, Real scale, const int* index,
                                const std::complex<Real>* grid, std::complex<Real>* output)
{
    gather_fft_coefficients<<<(size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK>>>(
        size, add, scale, index, reinterpret_cast<const thrust::complex<Real>*>(grid),
        reinterpret_cast<thrust::complex<Real>*>(output));
    cudaCheckOnDebug();
}

template void gather_fft_coefficients_gpu<float>(int, bool, float, const int*, const std::complex<float>*,
                                                std::complex<float>*);
template void gather_fft_coefficients_gpu<double>(int, bool, double, const int*, const std::complex<double>*,
                                                 std::complex<double>*);

template<class FPTYPE>
__global__ void set_3d_fft_box(
    const int npwk,
    const int* box_index,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < npwk)
    {
        int xx = box_index[idx];
        out[xx] = in[idx];
    }
}

template<class FPTYPE>
__global__ void set_recip_to_real_output(
    const int nrxx,
    const bool add,
    const FPTYPE factor,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= nrxx) {return;}
    if(add) {
        out[idx] += factor * in[idx];
    }
    else {
        out[idx] = in[idx];
    }
}

template<class FPTYPE>
__global__ void set_real_to_recip_output(
    const int npwk,
    const int nxyz,
    const bool add,
    const FPTYPE factor,
    const int* box_index,
    const thrust::complex<FPTYPE>* in,
    thrust::complex<FPTYPE>* out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= npwk) {return;}
    if(add) {
        out[idx] += factor / nxyz * in[box_index[idx]];
    }
    else {
        out[idx] = in[box_index[idx]] / nxyz;
    }
}

template <typename FPTYPE>
void set_3d_fft_box_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* /*dev*/,
                                                                    const int npwk,
                                                                    const int* box_index,
                                                                    const std::complex<FPTYPE>* in,
                                                                    std::complex<FPTYPE>* out)
{
    const int block = (npwk + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_3d_fft_box<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        npwk,
        box_index,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));

    cudaCheckOnDebug();
}

template <typename FPTYPE>
void set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* /*dev*/,
                                                                              const int nrxx,
                                                                              const bool add,
                                                                              const FPTYPE factor,
                                                                              const std::complex<FPTYPE>* in,
                                                                              std::complex<FPTYPE>* out)
{
    const int block = (nrxx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_recip_to_real_output<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        nrxx,
        add,
        factor,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));

    cudaCheckOnDebug();
}

template <typename FPTYPE>
void set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* /*dev*/,
                                                                              const int npwk,
                                                                              const int nxyz,
                                                                              const bool add,
                                                                              const FPTYPE factor,
                                                                              const int* box_index,
                                                                              const std::complex<FPTYPE>* in,
                                                                              std::complex<FPTYPE>* out)
{
    const int block = (npwk + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    set_real_to_recip_output<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
        npwk,
        nxyz,
        add,
        factor,
        box_index,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(in),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));

    cudaCheckOnDebug();
}

template struct set_3d_fft_box_op<float, base_device::DEVICE_GPU>;
template struct set_recip_to_real_output_op<float, base_device::DEVICE_GPU>;
template struct set_real_to_recip_output_op<float, base_device::DEVICE_GPU>;
template struct set_3d_fft_box_op<double, base_device::DEVICE_GPU>;
template struct set_recip_to_real_output_op<double, base_device::DEVICE_GPU>;
template struct set_real_to_recip_output_op<double, base_device::DEVICE_GPU>;

}  // namespace ModulePW
