#include "module_basis/module_pw/kernels/fused_veff.h"

#include <base/macros/macros.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cufftdx.hpp>
#include <thrust/complex.h>

#include <algorithm>

namespace ModulePW
{
namespace
{

using FFTBase = decltype(cufftdx::Size<200>() + cufftdx::Precision<double>()
                         + cufftdx::Type<cufftdx::fft_type::c2c>() + cufftdx::ElementsPerThread<10>()
                         + cufftdx::FFTsPerBlock<4>() + cufftdx::SM<700>() + cufftdx::Block());
using FFTInverse = decltype(FFTBase() + cufftdx::Direction<cufftdx::fft_direction::inverse>());
using FFTForward = decltype(FFTBase() + cufftdx::Direction<cufftdx::fft_direction::forward>());

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

__global__ void pack_spinor(int npw, int stride, int plane, int left, int right, const int* index,
                            const double2* input, double2* up, double2* down)
{
    const int ig = blockIdx.x * blockDim.x + threadIdx.x;
    if (ig < npw)
    {
        const int compact = compact_grid_index(index[ig], plane, left, right);
        up[compact] = input[ig];
        down[compact] = input[stride + ig];
    }
}

__global__ void add_spinor(int npw, int stride, int plane, int left, int right, int nxyz, const int* index,
                           const thrust::complex<double>* up, const thrust::complex<double>* down,
                           thrust::complex<double>* output)
{
    const int ig = blockIdx.x * blockDim.x + threadIdx.x;
    if (ig < npw)
    {
        const int compact = compact_grid_index(index[ig], plane, left, right);
        output[ig] += 1.0 / nxyz * up[compact];
        output[stride + ig] += 1.0 / nxyz * down[compact];
    }
}

// Cooperative transfers coalesce the adjacent y-z lines. Only the retained
// x slabs occupy the shared tile; absent Fourier coefficients are zero.
template <class FFT, bool Load>
__device__ void transfer_spinor_field(const cufftDoubleComplex* input, cufftDoubleComplex* output,
                              typename FFT::value_type* values, typename FFT::value_type* tile,
                              int plane, int left, int right)
{
    constexpr int nx = cufftdx::size_of<FFT>::value;
    constexpr int batches = FFT::ffts_per_block;
    const int tid = threadIdx.y * blockDim.x + threadIdx.x;
    if (!Load)
    {
        __syncthreads();
        for (int i = 0; i < FFT::elements_per_thread; ++i)
        {
            const int x = threadIdx.x + i * FFT::stride;
            if (x < left || (x >= right && x < nx))
            {
                const int compact_x = x < left ? x : x - right + left;
                tile[compact_x * batches + threadIdx.y] = values[i];
            }
        }
        __syncthreads();
    }
    for (int i = tid; i < (left + nx - right) * batches; i += blockDim.x * blockDim.y)
    {
        const int compact_x = i / batches;
        const int line = blockIdx.x * batches + i % batches;
        const bool valid = line < plane;
        if (Load)
        {
            const auto value = valid ? input[compact_x * plane + line] : make_cuDoubleComplex(0, 0);
            tile[i] = typename FFT::value_type(value.x, value.y);
        }
        else if (valid)
        {
            output[compact_x * plane + line] = make_cuDoubleComplex(tile[i].x, tile[i].y);
        }
    }
    if (Load)
    {
        __syncthreads();
        for (int i = 0; i < FFT::elements_per_thread; ++i)
        {
            const int x = threadIdx.x + i * FFT::stride;
            const bool retained = x < left || (x >= right && x < nx);
            const int compact_x = x < left ? x : x - right + left;
            values[i] = retained ? tile[compact_x * batches + threadIdx.y] : typename FFT::value_type(0, 0);
        }
    }
    __syncthreads();
}

template <class FFTInverse, class FFTForward>
__global__ void fused_spinor_potential(const cufftDoubleComplex* input_up, const cufftDoubleComplex* input_down,
                               cufftDoubleComplex* output_up, cufftDoubleComplex* output_down,
                               const double* potential, int plane, int left, int right)
{
    using Complex = typename FFTInverse::value_type;
    constexpr int nx = cufftdx::size_of<FFTInverse>::value;
    const int line = blockIdx.x * FFTInverse::ffts_per_block + threadIdx.y;
    const int size = nx * plane;
    extern __shared__ __align__(16) Complex shared[];
    Complex up[FFTInverse::storage_size], down[FFTInverse::storage_size];
    transfer_spinor_field<FFTInverse, true>(input_up, nullptr, up, shared, plane, left, right);
    FFTInverse().execute(up, shared);
    __syncthreads();
    transfer_spinor_field<FFTInverse, true>(input_down, nullptr, down, shared, plane, left, right);
    FFTInverse().execute(down, shared);
    for (int j = 0; j < FFTInverse::elements_per_thread; ++j)
    {
        const int x = threadIdx.x + j * FFTInverse::stride;
        if (line < plane && x < nx)
        {
            const int i = x * plane + line;
            const double v0 = potential[i], v1 = potential[size + i];
            const double v2 = potential[2 * size + i], v3 = potential[3 * size + i];
            const Complex a = up[j], b = down[j];
            const Complex cross_up {b.x * v1 + b.y * v2, b.y * v1 - b.x * v2};
            const Complex cross_down {a.x * v1 - a.y * v2, a.y * v1 + a.x * v2};
            up[j] = Complex(a.x * (v0 + v3) + cross_up.x, a.y * (v0 + v3) + cross_up.y);
            down[j] = Complex(b.x * (v0 - v3) + cross_down.x, b.y * (v0 - v3) + cross_down.y);
        }
    }
    __syncthreads();
    FFTForward().execute(up, shared);
    __syncthreads();
    FFTForward().execute(down, shared);
    transfer_spinor_field<FFTForward, false>(nullptr, output_up, up, shared, plane, left, right);
    transfer_spinor_field<FFTForward, false>(nullptr, output_down, down, shared, plane, left, right);
}

} // namespace

bool fused_veff_supported()
{
    int device = 0;
    cudaError_t status = cudaGetDevice(&device);
    cudaErrcheck(status);
    cudaDeviceProp properties;
    status = cudaGetDeviceProperties(&properties, device);
    cudaErrcheck(status);
    return properties.major == 7 && properties.minor == 0;
}

void pack_spinor_grid(int npw, int stride, int plane, int left, int right, const int* index,
                     const std::complex<double>* input, std::complex<double>* up, std::complex<double>* down)
{
    pack_spinor<<<(npw + 255) / 256, 256>>>(npw, stride, plane, left, right, index,
        reinterpret_cast<const double2*>(input), reinterpret_cast<double2*>(up), reinterpret_cast<double2*>(down));
    check_launch();
}

void add_spinor_grid(int npw, int stride, int plane, int left, int right, int nxyz, const int* index,
                    const std::complex<double>* up, const std::complex<double>* down, std::complex<double>* output)
{
    add_spinor<<<(npw + 255) / 256, 256>>>(npw, stride, plane, left, right, nxyz, index,
        reinterpret_cast<const thrust::complex<double>*>(up), reinterpret_cast<const thrust::complex<double>*>(down),
        reinterpret_cast<thrust::complex<double>*>(output));
    check_launch();
}

void fused_veff_x(int plane, int left, int right, const double* potential,
                  std::complex<double>* up, std::complex<double>* down)
{
    static_assert(!FFTInverse::requires_workspace && !FFTForward::requires_workspace,
                  "The fused transform requires no global workspace.");
    const size_t tile_bytes = (left + 200 - right) * FFTInverse::ffts_per_block * sizeof(cufftDoubleComplex);
    const size_t fft_bytes = std::max(FFTInverse::shared_memory_size, FFTForward::shared_memory_size);
    const size_t shared_bytes = std::max(tile_bytes, fft_bytes);
    auto* up_grid = reinterpret_cast<cufftDoubleComplex*>(up);
    auto* down_grid = reinterpret_cast<cufftDoubleComplex*>(down);
    fused_spinor_potential<FFTInverse, FFTForward>
        <<<(plane + FFTInverse::ffts_per_block - 1) / FFTInverse::ffts_per_block,
           FFTInverse::block_dim, shared_bytes>>>(up_grid, down_grid, up_grid, down_grid, potential, plane, left, right);
    check_launch();
}

} // namespace ModulePW
