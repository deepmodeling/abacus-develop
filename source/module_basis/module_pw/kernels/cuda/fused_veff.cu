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

template <unsigned int Batches>
using FFTBase = decltype(cufftdx::Size<200>() + cufftdx::Precision<double>()
                         + cufftdx::Type<cufftdx::fft_type::c2c>() + cufftdx::ElementsPerThread<10>()
                         + cufftdx::FFTsPerBlock<Batches>() + cufftdx::SM<700>() + cufftdx::Block());
using FFTInverse = decltype(FFTBase<6>() + cufftdx::Direction<cufftdx::fft_direction::inverse>());
using FFTForward = decltype(FFTBase<6>() + cufftdx::Direction<cufftdx::fft_direction::forward>());

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
                            const double2* input, double2* up, double2* down)
{
    const int ig = blockIdx.x * blockDim.x + threadIdx.x;
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
                           thrust::complex<double>* output)
{
    const int ig = blockIdx.x * blockDim.x + threadIdx.x;
    if (ig < npw)
    {
        const int compact = Compact ? compact_grid_index(index[ig], plane, left, right) : index[ig];
        output[ig] += 1.0 / nxyz * up[compact];
        output[stride + ig] += 1.0 / nxyz * down[compact];
    }
}

// The inverse transform expands occupied x-y sticks into complete y-z slabs.
// The forward transform retains only the sticks required by the plane-wave basis.
template <class FFT, bool Inverse>
__global__ void transform_y_sticks(cufftDoubleComplex* sticks, cufftDoubleComplex* compact,
                             const int* reverse, int nz, int left, int right)
{
    using Complex = typename FFT::value_type;
    constexpr int ny = 200, batches = FFT::ffts_per_block;
    const int cx = blockIdx.y;
    const int x = cx < left ? cx : cx + right - left;
    const int tid = threadIdx.y * blockDim.x + threadIdx.x;
    extern __shared__ __align__(16) Complex shared[];
    Complex values[FFT::storage_size];
    for (int i = tid; i < ny * batches; i += blockDim.x * blockDim.y)
    {
        const int y = i / batches, iz = blockIdx.x * batches + i % batches;
        const int stick = reverse[x * ny + y];
        const int index = Inverse ? stick * nz + iz : (cx * ny + y) * nz + iz;
        const bool valid = iz < nz && (!Inverse || stick >= 0);
        const auto v = valid ? (Inverse ? sticks[index] : compact[index]) : make_cuDoubleComplex(0, 0);
        shared[i] = Complex(v.x, v.y);
    }
    __syncthreads();
    for (int i = 0; i < FFT::elements_per_thread; ++i)
    {
        const int y = threadIdx.x + i * FFT::stride;
        values[i] = y < ny ? shared[y * batches + threadIdx.y] : Complex(0, 0);
    }
    __syncthreads();
    FFT().execute(values, shared);
    __syncthreads();
    for (int i = 0; i < FFT::elements_per_thread; ++i)
    {
        const int y = threadIdx.x + i * FFT::stride;
        if (y < ny) { shared[y * batches + threadIdx.y] = values[i]; }
    }
    __syncthreads();
    for (int i = tid; i < ny * batches; i += blockDim.x * blockDim.y)
    {
        const int y = i / batches, iz = blockIdx.x * batches + i % batches;
        const int stick = reverse[x * ny + y];
        if (iz < nz && (Inverse || stick >= 0))
        {
            const auto v = make_cuDoubleComplex(shared[i].x, shared[i].y);
            if (Inverse) { compact[(cx * ny + y) * nz + iz] = v; }
            else { sticks[stick * nz + iz] = v; }
        }
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

template <class FFT>
__global__ void fused_spinor_density(const cufftDoubleComplex* input_up, const cufftDoubleComplex* input_down,
                              double* rho, double weight, int plane, int left, int right)
{
    using Complex = typename FFT::value_type;
    constexpr int nx = cufftdx::size_of<FFT>::value;
    const int size = nx * plane;
    extern __shared__ __align__(16) Complex shared[];
    Complex up[FFT::storage_size], down[FFT::storage_size];
    transfer_spinor_field<FFT, true>(input_up, nullptr, up, shared, plane, left, right);
    FFT().execute(up, shared);
    __syncthreads();
    transfer_spinor_field<FFT, true>(input_down, nullptr, down, shared, plane, left, right);
    FFT().execute(down, shared);
    auto* tile = reinterpret_cast<double*>(shared);
    for (int spin = 0; spin < 4; ++spin)
    {
        __syncthreads();
        for (int j = 0; j < FFT::elements_per_thread; ++j)
        {
            const int x = threadIdx.x + j * FFT::stride;
            if (x < nx)
            {
                const Complex a = up[j], b = down[j];
                const double norm_up = a.x * a.x + a.y * a.y;
                const double norm_down = b.x * b.x + b.y * b.y;
                double value = weight * (norm_up + norm_down);
                if (spin == 1) { value = weight * 2.0 * (a.x * b.x + a.y * b.y); }
                else if (spin == 2) { value = weight * 2.0 * (a.x * b.y - b.x * a.y); }
                else if (spin == 3) { value = weight * (norm_up - norm_down); }
                tile[x * FFT::ffts_per_block + threadIdx.y] = value;
            }
        }
        __syncthreads();
        const int tid = threadIdx.y * blockDim.x + threadIdx.x;
        for (int i = tid; i < nx * FFT::ffts_per_block; i += blockDim.x * blockDim.y)
        {
            const int x = i / FFT::ffts_per_block;
            const int output_line = blockIdx.x * FFT::ffts_per_block + i % FFT::ffts_per_block;
            if (output_line < plane) { rho[spin * size + x * plane + output_line] += tile[i]; }
        }
    }
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
    pack_spinor<true><<<(npw + 255) / 256, 256>>>(npw, stride, plane, left, right, index,
        reinterpret_cast<const double2*>(input), reinterpret_cast<double2*>(up), reinterpret_cast<double2*>(down));
    check_launch();
}

void add_spinor_grid(int npw, int stride, int plane, int left, int right, int nxyz, const int* index,
                    const std::complex<double>* up, const std::complex<double>* down, std::complex<double>* output)
{
    add_spinor<true><<<(npw + 255) / 256, 256>>>(npw, stride, plane, left, right, nxyz, index,
        reinterpret_cast<const thrust::complex<double>*>(up), reinterpret_cast<const thrust::complex<double>*>(down),
        reinterpret_cast<thrust::complex<double>*>(output));
    check_launch();
}

void pack_spinor_sticks(int npw, int stride, const int* index, const std::complex<double>* input,
                        std::complex<double>* up, std::complex<double>* down)
{
    pack_spinor<false><<<(npw + 255) / 256, 256>>>(npw, stride, 0, 0, 0, index,
        reinterpret_cast<const double2*>(input), reinterpret_cast<double2*>(up), reinterpret_cast<double2*>(down));
    check_launch();
}

void add_spinor_sticks(int npw, int stride, int nxyz, const int* index,
                       const std::complex<double>* up, const std::complex<double>* down, std::complex<double>* output)
{
    add_spinor<false><<<(npw + 255) / 256, 256>>>(npw, stride, 0, 0, 0, nxyz, index,
        reinterpret_cast<const thrust::complex<double>*>(up), reinterpret_cast<const thrust::complex<double>*>(down),
        reinterpret_cast<thrust::complex<double>*>(output));
    check_launch();
}

void transform_stick_slabs(std::complex<double>* sticks, std::complex<double>* slabs,
                           const int* reverse, int nz, int left, int right, bool inverse)
{
    using Inverse = decltype(FFTBase<4>() + cufftdx::Direction<cufftdx::fft_direction::inverse>());
    using Forward = decltype(FFTBase<4>() + cufftdx::Direction<cufftdx::fft_direction::forward>());
    constexpr size_t tile_bytes = 200 * Inverse::ffts_per_block * sizeof(cufftDoubleComplex);
    const size_t shared_bytes = std::max(tile_bytes, size_t(Inverse::shared_memory_size));
    const dim3 blocks((nz + Inverse::ffts_per_block - 1) / Inverse::ffts_per_block, left + 200 - right);
    auto* stick_data = reinterpret_cast<cufftDoubleComplex*>(sticks);
    auto* slab_data = reinterpret_cast<cufftDoubleComplex*>(slabs);
    if (inverse)
    {
        transform_y_sticks<Inverse, true><<<blocks, Inverse::block_dim, shared_bytes>>>(
            stick_data, slab_data, reverse, nz, left, right);
    }
    else
    {
        transform_y_sticks<Forward, false><<<blocks, Forward::block_dim, shared_bytes>>>(
            stick_data, slab_data, reverse, nz, left, right);
    }
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

void fused_density_x(int plane, int left, int right, double weight, double* density,
                     const std::complex<double>* up, const std::complex<double>* down)
{
    using Inverse = decltype(FFTBase<8>() + cufftdx::Direction<cufftdx::fft_direction::inverse>());
    const size_t transfer_bytes = (left + 200 - right) * Inverse::ffts_per_block * sizeof(cufftDoubleComplex);
    const size_t density_bytes = 200 * Inverse::ffts_per_block * sizeof(double);
    const size_t shared_bytes = std::max(std::max(transfer_bytes, density_bytes), size_t(Inverse::shared_memory_size));
    fused_spinor_density<Inverse><<<(plane + Inverse::ffts_per_block - 1) / Inverse::ffts_per_block,
                                     Inverse::block_dim, shared_bytes>>>(
        reinterpret_cast<const cufftDoubleComplex*>(up), reinterpret_cast<const cufftDoubleComplex*>(down),
        density, weight, plane, left, right);
    check_launch();
}

} // namespace ModulePW
