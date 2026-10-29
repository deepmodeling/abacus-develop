#ifndef MODULEPW_FUSED_FFT_TEMPLATE_CUH
#define MODULEPW_FUSED_FFT_TEMPLATE_CUH

#include "module_basis/module_pw/kernels/fused_veff.h"
#include <base/macros/macros.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cufftdx.hpp>
#include <algorithm>

namespace ModulePW
{
namespace
{

// The inverse transform expands occupied x-y sticks into complete y-z slabs.
// The forward transform retains only the sticks required by the plane-wave basis.
template <class FFT, bool Inverse>
__global__ void transform_y_sticks(cufftDoubleComplex* sticks, cufftDoubleComplex* compact,
                             const int* reverse, int nx, int nz, int left, int right, int stick_size)
{
    using Complex = typename FFT::value_type;
    constexpr int ny = cufftdx::size_of<FFT>::value, batches = FFT::ffts_per_block;
    sticks += blockIdx.z * stick_size;
    compact += blockIdx.z * (left + nx - right) * ny * nz;
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
                              int plane, int left, int right, int spatial_block)
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
        const int line = spatial_block * batches + i % batches;
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
                               const double* potential, int plane, int left, int right, int bands)
{
    const int band = blockIdx.x % bands;
    const int spatial_block = blockIdx.x / bands;
    const int count = (left + cufftdx::size_of<FFTInverse>::value - right) * plane;
    input_up += band * count;
    input_down += band * count;
    output_up += band * count;
    output_down += band * count;
    using Complex = typename FFTInverse::value_type;
    constexpr int nx = cufftdx::size_of<FFTInverse>::value;
    const int line = spatial_block * FFTInverse::ffts_per_block + threadIdx.y;
    const int size = nx * plane;
    extern __shared__ __align__(16) Complex shared[];
    Complex up[FFTInverse::storage_size], down[FFTInverse::storage_size];
    transfer_spinor_field<FFTInverse, true>(input_up, nullptr, up, shared, plane, left, right, spatial_block);
    FFTInverse().execute(up, shared);
    __syncthreads();
    transfer_spinor_field<FFTInverse, true>(input_down, nullptr, down, shared, plane, left, right, spatial_block);
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
    transfer_spinor_field<FFTForward, false>(nullptr, output_up, up, shared, plane, left, right, spatial_block);
    transfer_spinor_field<FFTForward, false>(nullptr, output_down, down, shared, plane, left, right, spatial_block);
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
    transfer_spinor_field<FFT, true>(input_up, nullptr, up, shared, plane, left, right, blockIdx.x);
    FFT().execute(up, shared);
    __syncthreads();
    transfer_spinor_field<FFT, true>(input_down, nullptr, down, shared, plane, left, right, blockIdx.x);
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


template <class FFT, int Limit>
constexpr int fft_elements_at_most()
{
    using Candidate = decltype(FFT() + cufftdx::ElementsPerThread<Limit>());
    if constexpr (cufftdx::is_supported<Candidate, 700>::value)
    {
        return Limit;
    }
    else if constexpr (Limit > 1)
    {
        return fft_elements_at_most<FFT, Limit - 1>();
    }
    return FFT::elements_per_thread;
}

template <int Length>
struct FusedTransform
{
    using Common = decltype(cufftdx::Size<Length>() + cufftdx::Precision<double>()
                            + cufftdx::Type<cufftdx::fft_type::c2c>()
                            + cufftdx::SM<700>() + cufftdx::Block());
    using Default = decltype(Common() + cufftdx::Direction<cufftdx::fft_direction::inverse>());
    static constexpr int suggested = Default::elements_per_thread;
    using Half = decltype(Default() + cufftdx::ElementsPerThread<suggested / 2>());
    // Two simultaneous spin fields retain twice the FFT register state. Select
    // the largest legal count up to ten, retaining the library default only
    // when no smaller descriptor exists. Single-field stick transforms use
    // their own register budget.
    static constexpr int elements = fft_elements_at_most<Default, std::min(10, suggested)>();
    static constexpr int single_elements = suggested > 16 && suggested % 2 == 0
                                          && cufftdx::is_supported<Half, 700>::value ? suggested / 2 : suggested;
    using Base = decltype(Common() + cufftdx::ElementsPerThread<elements>());
    using Geometry = decltype(Base() + cufftdx::Direction<cufftdx::fft_direction::inverse>());
    static constexpr int stride = Geometry::stride;
    static constexpr int batches = std::min({16, std::max(1, 128 / stride), 49152 / (Length * 16)});
    using Inverse = decltype(Base() + cufftdx::Direction<cufftdx::fft_direction::inverse>()
                            + cufftdx::FFTsPerBlock<batches>());
    using Forward = decltype(Base() + cufftdx::Direction<cufftdx::fft_direction::forward>()
                            + cufftdx::FFTsPerBlock<batches>());

    static void potential(int plane, int left, int right, const double* values,
                           std::complex<double>* up, std::complex<double>* down, int bands)
    {
        const size_t tile = size_t(left + Length - right) * batches * sizeof(cufftDoubleComplex);
        const size_t bytes = std::max(tile, size_t(Inverse::shared_memory_size));
        auto* a = reinterpret_cast<cufftDoubleComplex*>(up);
        auto* b = reinterpret_cast<cufftDoubleComplex*>(down);
        fused_spinor_potential<Inverse, Forward>
            <<<bands * ((plane + batches - 1) / batches), Inverse::block_dim, bytes>>>(
                a, b, a, b, values, plane, left, right, bands);
    }

    static void density(int plane, int left, int right, double weight, double* rho,
                         const std::complex<double>* up, const std::complex<double>* down)
    {
        const size_t transfer = size_t(left + Length - right) * batches * sizeof(cufftDoubleComplex);
        const size_t density = size_t(Length) * batches * sizeof(double);
        const size_t bytes = std::max({transfer, density, size_t(Inverse::shared_memory_size)});
        fused_spinor_density<Inverse><<<(plane + batches - 1) / batches, Inverse::block_dim, bytes>>>(
            reinterpret_cast<const cufftDoubleComplex*>(up), reinterpret_cast<const cufftDoubleComplex*>(down),
            rho, weight, plane, left, right);
    }

    static void sticks(std::complex<double>* input, std::complex<double>* output, const int* reverse,
                        int nx, int nz, int left, int right, bool inverse, int bands, int stick_size)
    {
        using YBase = decltype(Common() + cufftdx::ElementsPerThread<single_elements>());
        using YGeometry = decltype(YBase() + cufftdx::Direction<cufftdx::fft_direction::inverse>());
        constexpr int y_stride = YGeometry::stride;
        constexpr int y_batches = std::min({16, std::max(1, 256 / y_stride), 49152 / (Length * 16)});
        using YInverse = decltype(YBase() + cufftdx::Direction<cufftdx::fft_direction::inverse>()
                                 + cufftdx::FFTsPerBlock<y_batches>());
        using YForward = decltype(YBase() + cufftdx::Direction<cufftdx::fft_direction::forward>()
                                 + cufftdx::FFTsPerBlock<y_batches>());
        const size_t tile = size_t(Length) * y_batches * sizeof(cufftDoubleComplex);
        const size_t bytes = std::max(tile, size_t(YInverse::shared_memory_size));
        const dim3 blocks((nz + y_batches - 1) / y_batches, left + nx - right, bands);
        auto* a = reinterpret_cast<cufftDoubleComplex*>(input);
        auto* b = reinterpret_cast<cufftDoubleComplex*>(output);
        if (inverse)
        {
            transform_y_sticks<YInverse, true><<<blocks, YInverse::block_dim, bytes>>>(
                a, b, reverse, nx, nz, left, right, stick_size);
        }
        else
        {
            transform_y_sticks<YForward, false><<<blocks, YForward::block_dim, bytes>>>(
                a, b, reverse, nx, nz, left, right, stick_size);
        }
    }
};

template <int Length>
const FusedFFT* make_fused_fft()
{
    using Transform = FusedTransform<Length>;
    if constexpr (!Transform::Default::requires_workspace)
    {
        static const FusedFFT functions = {Transform::potential, Transform::density, Transform::sticks};
        return &functions;
    }
    return nullptr;
}

} // namespace
} // namespace ModulePW

#endif
