#include "source_hsolver/kernels/ppcg_kernel_op.h"
#include "source_hsolver/ppcg_workspace.h"
#include <base/macros/macros.h>
#include <cuda_runtime.h>
#include <thrust/complex.h>

namespace hsolver
{
namespace
{
enum
{
    block_threads = 256
};

template <typename Real>
__global__ void column_products_kernel(
    const int dim,
    const int ld,
    const thrust::complex<Real> *x,
    const thrust::complex<Real> *y,
    double *products)
{
    __shared__ double sums[block_threads];
    const int column = blockIdx.x;
    double sum = 0.0;
    for (int row = threadIdx.x; row < dim; row += blockDim.x)
    {
        const int offset = row + column * ld;
        const auto product = thrust::conj(x[offset]) * y[offset];
        sum += double(product.real());
    }
    sums[threadIdx.x] = sum;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2)
    {
        if (threadIdx.x < stride)
        {
            sums[threadIdx.x] += sums[threadIdx.x + stride];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0)
    {
        products[column] = sums[0];
    }
}

template <typename Real>
__global__ void residual_kernel(
    const int dim,
    const int ld,
    const int count,
    const thrust::complex<Real> *hx,
    const thrust::complex<Real> *sx,
    const Real *eigenvalues,
    thrust::complex<Real> *output)
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < dim * count)
    {
        const int row = index % dim;
        const int column = index / dim;
        const int offset = row + column * ld;
        const Real real_part = fma(-eigenvalues[column], sx[offset].real(), hx[offset].real());
        const Real imaginary_part = fma(-eigenvalues[column], sx[offset].imag(), hx[offset].imag());
        output[offset] = thrust::complex<Real>(real_part, imaginary_part);
    }
}

template <typename Real>
__global__ void precondition_kernel(
    const int dim,
    const int ld,
    const int count,
    const Real *diagonal,
    thrust::complex<Real> *values)
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < dim * count)
    {
        const int row = index % dim;
        const int column = index / dim;
        const Real floor = Real(1.0e-12);
        const Real divisor = diagonal[row] > floor ? diagonal[row] : floor;
        values[row + column * ld] /= divisor;
    }
}

template <typename Real>
__global__ void scale_columns_kernel(
    const int dim,
    const int ld,
    const int count,
    const double *scales,
    thrust::complex<Real> *x,
    thrust::complex<Real> *sx,
    thrust::complex<Real> *hx)
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < dim * count)
    {
        const int row = index % dim;
        const int column = index / dim;
        const int offset = row + column * ld;
        const Real scale = Real(scales[column]);
        x[offset] *= scale;
        sx[offset] *= scale;
        hx[offset] *= scale;
    }
}
} // namespace

template <typename T>
void PPCGVectorOps<T, base_device::DEVICE_GPU>::column_products(
    const int dim,
    const int ld,
    const int count,
    const T *x,
    const T *y,
    double *host_products)
{
    if (count == 0)
    {
        return;
    }
    PPCGDeviceBuffer<double> products(count, 0.0);
    const auto *device_x = reinterpret_cast<const thrust::complex<Real> *>(x);
    const auto *device_y = reinterpret_cast<const thrust::complex<Real> *>(y);
    column_products_kernel<<<count, block_threads>>>(dim, ld, device_x, device_y, products.data());
    CHECK_CUDA(cudaGetLastError());
    base_device::memory::synchronize_memory_op<double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(
        host_products, products.data(), count);
}

template <typename T>
void PPCGVectorOps<T, base_device::DEVICE_GPU>::residual(
    const int dim,
    const int ld,
    const int count,
    const T *hx,
    const T *sx,
    const Real *eigenvalues,
    T *output)
{
    if (count == 0)
    {
        return;
    }
    const int blocks = 1 + (dim * count - 1) / block_threads;
    const auto *device_hx = reinterpret_cast<const thrust::complex<Real> *>(hx);
    const auto *device_sx = reinterpret_cast<const thrust::complex<Real> *>(sx);
    auto *device_output = reinterpret_cast<thrust::complex<Real> *>(output);
    residual_kernel<<<blocks, block_threads>>>(dim, ld, count, device_hx, device_sx, eigenvalues, device_output);
    CHECK_CUDA(cudaGetLastError());
}

template <typename T>
void PPCGVectorOps<T, base_device::DEVICE_GPU>::precondition(
    const int dim,
    const int ld,
    const int count,
    const Real *diagonal,
    T *values)
{
    if (count == 0)
    {
        return;
    }
    const int blocks = 1 + (dim * count - 1) / block_threads;
    auto *device_values = reinterpret_cast<thrust::complex<Real> *>(values);
    precondition_kernel<<<blocks, block_threads>>>(dim, ld, count, diagonal, device_values);
    CHECK_CUDA(cudaGetLastError());
}

template <typename T>
void PPCGVectorOps<T, base_device::DEVICE_GPU>::scale_columns(
    const int dim,
    const int ld,
    const int count,
    const double *scales,
    T *x,
    T *sx,
    T *hx)
{
    if (count == 0)
    {
        return;
    }
    const int blocks = 1 + (dim * count - 1) / block_threads;
    auto *device_x = reinterpret_cast<thrust::complex<Real> *>(x);
    auto *device_sx = reinterpret_cast<thrust::complex<Real> *>(sx);
    auto *device_hx = reinterpret_cast<thrust::complex<Real> *>(hx);
    scale_columns_kernel<<<blocks, block_threads>>>(dim, ld, count, scales, device_x, device_sx, device_hx);
    CHECK_CUDA(cudaGetLastError());
}

template struct PPCGVectorOps<std::complex<float>, base_device::DEVICE_GPU>;
template struct PPCGVectorOps<std::complex<double>, base_device::DEVICE_GPU>;
} // namespace hsolver
