#include "module_elecstate/kernels/hartree_op.h"

#include <base/macros/macros.h>
#include <cuda_runtime.h>
#include <thrust/complex.h>

namespace elecstate
{
namespace
{

constexpr int block_size = 256;

__global__ void pack_density(int size, int channels, const double* density, thrust::complex<double>* grid)
{
    const int ir = blockIdx.x * blockDim.x + threadIdx.x;
    if (ir < size)
    {
        double value = 0;
        for (int spin = 0; spin < channels; ++spin) { value += density[spin * size + ir]; }
        grid[ir] = thrust::complex<double>(value, 0);
    }
}

__global__ void poisson_coefficients(int npw, const int* index, const double* gg, double scale,
                                     double numerator, double tpiba2, const thrust::complex<double>* grid,
                                     thrust::complex<double>* potential, double* energy)
{
    const int ig = blockIdx.x * blockDim.x + threadIdx.x;
    double contribution = 0;
    if (ig < npw)
    {
        potential[ig] = thrust::complex<double>(0, 0);
        if (gg[ig] >= 1e-8)
        {
            const auto density = scale * grid[index[ig]];
            const double factor = numerator / (tpiba2 * gg[ig]);
            potential[ig] = factor * density;
            contribution = (thrust::conj(density) * density).real() * factor;
        }
    }
    __shared__ double partial[block_size];
    partial[threadIdx.x] = contribution;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2)
    {
        if (threadIdx.x < stride) { partial[threadIdx.x] += partial[threadIdx.x + stride]; }
        __syncthreads();
    }
    if (threadIdx.x == 0) { energy[blockIdx.x] = partial[0]; }
}

void check_launch()
{
    const cudaError_t status = cudaGetLastError();
    cudaErrcheck(status);
    cudaCheckOnDebug();
}

} // namespace

int hartree_reduction_size_gpu(int npw)
{
    return (npw + block_size - 1) / block_size;
}

void hartree_pack_density_gpu(int size, int channels, const double* density, std::complex<double>* grid)
{
    pack_density<<<(size + block_size - 1) / block_size, block_size>>>(
        size, channels, density, reinterpret_cast<thrust::complex<double>*>(grid));
    check_launch();
}

void hartree_coefficients_gpu(int npw, const int* index, const double* gg, double scale,
                              double coulomb_numerator, double tpiba2, const std::complex<double>* grid,
                              std::complex<double>* potential, double* energy)
{
    poisson_coefficients<<<hartree_reduction_size_gpu(npw), block_size>>>(
        npw, index, gg, scale, coulomb_numerator, tpiba2, reinterpret_cast<const thrust::complex<double>*>(grid),
        reinterpret_cast<thrust::complex<double>*>(potential), energy);
    check_launch();
}

} // namespace elecstate
