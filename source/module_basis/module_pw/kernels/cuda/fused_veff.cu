#include "module_basis/module_pw/kernels/fused_veff.h"

#include <base/macros/macros.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <thrust/complex.h>

namespace ModulePW
{
namespace
{

void check_launch()
{
    const cudaError_t status = cudaGetLastError();
    cudaErrcheck(status);
    cudaCheckOnDebug();
}

__device__ int compact_grid_index(int index, int plane, int left, int right)
{
    return index - (index / plane >= right ? (right - left) * plane : 0);
}

template <bool Compact>
__global__ void pack_spinor(int npw, int stride, int plane, int left, int right, const int* index,
                            const double2* input, double2* up, double2* down, int field_size)
{
    const int ig = blockIdx.x * blockDim.x + threadIdx.x;
    up += blockIdx.y * field_size;
    down += blockIdx.y * field_size;
    input += blockIdx.y * 2 * stride;
    if (ig < npw)
    {
        const int compact = Compact ? compact_grid_index(index[ig], plane, left, right) : index[ig];
        up[compact] = input[ig];
        down[compact] = input[stride + ig];
    }
}

template <bool Compact>
__global__ void add_spinor(int npw, int stride, int plane, int left, int right, int nxyz, const int* index,
                           const thrust::complex<double>* up, const thrust::complex<double>* down,
                           thrust::complex<double>* output, int field_size)
{
    const int ig = blockIdx.x * blockDim.x + threadIdx.x;
    up += blockIdx.y * field_size;
    down += blockIdx.y * field_size;
    output += blockIdx.y * 2 * stride;
    if (ig < npw)
    {
        const int compact = Compact ? compact_grid_index(index[ig], plane, left, right) : index[ig];
        output[ig] += 1.0 / nxyz * up[compact];
        output[stride + ig] += 1.0 / nxyz * down[compact];
    }
}

} // namespace

bool fused_veff_supported()
{
    return fused_fft(16) != nullptr;
}

void pack_spinor_grid(int npw, int stride, int plane, int left, int right, const int* index,
                     const std::complex<double>* input, std::complex<double>* up, std::complex<double>* down)
{
    pack_spinor<true><<<(npw + 255) / 256, 256>>>(npw, stride, plane, left, right, index,
        reinterpret_cast<const double2*>(input), reinterpret_cast<double2*>(up), reinterpret_cast<double2*>(down), 0);
    check_launch();
}

void add_spinor_grid(int npw, int stride, int plane, int left, int right, int nxyz, const int* index,
                    const std::complex<double>* up, const std::complex<double>* down, std::complex<double>* output)
{
    add_spinor<true><<<(npw + 255) / 256, 256>>>(npw, stride, plane, left, right, nxyz, index,
        reinterpret_cast<const thrust::complex<double>*>(up), reinterpret_cast<const thrust::complex<double>*>(down),
        reinterpret_cast<thrust::complex<double>*>(output), 0);
    check_launch();
}

void pack_spinor_sticks(int npw, int stride, const int* index, const std::complex<double>* input,
                        std::complex<double>* up, std::complex<double>* down, int bands, int stick_size)
{
    pack_spinor<false><<<dim3((npw + 255) / 256, bands), 256>>>(npw, stride, 0, 0, 0, index,
        reinterpret_cast<const double2*>(input), reinterpret_cast<double2*>(up), reinterpret_cast<double2*>(down), stick_size);
    check_launch();
}

void add_spinor_sticks(int npw, int stride, int nxyz, const int* index,
                       const std::complex<double>* up, const std::complex<double>* down, std::complex<double>* output, int bands, int stick_size)
{
    add_spinor<false><<<dim3((npw + 255) / 256, bands), 256>>>(npw, stride, 0, 0, 0, nxyz, index,
        reinterpret_cast<const thrust::complex<double>*>(up), reinterpret_cast<const thrust::complex<double>*>(down),
        reinterpret_cast<thrust::complex<double>*>(output), stick_size);
    check_launch();
}

void transform_stick_slabs(std::complex<double>* sticks, std::complex<double>* slabs,
                           const int* reverse, int nx, int ny, int nz, int left, int right,
                           bool inverse, int bands, int stick_size)
{
    fused_fft(ny)->sticks(sticks, slabs, reverse, nx, nz, left, right, inverse, bands, stick_size);
    check_launch();
}

void fused_veff_x(int nx, int plane, int left, int right, const double* potential,
                  std::complex<double>* up, std::complex<double>* down, int bands)
{
    fused_fft(nx)->potential(plane, left, right, potential, up, down, bands);
    check_launch();
}

void fused_density_x(int nx, int plane, int left, int right, double weight, double* density,
                     const std::complex<double>* up, const std::complex<double>* down)
{
    fused_density_x(nx, plane, left, right, &weight, 1, density, up, down);
}

void fused_density_x(int nx, int plane, int left, int right, const double* weights, int bands,
                     double* density, const std::complex<double>* up, const std::complex<double>* down)
{
    fused_fft(nx)->density(plane, left, right, weights, bands, density, up, down);
    check_launch();
}

} // namespace ModulePW
