#include "source_hsolver/linear_algebra.h"

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_device.h"
#include "source_hsolver/kernels/linear_op.h"

#include <cmath>
#include <limits>
#include <type_traits>

namespace hsolver
{
std::vector<std::complex<double>> linear_gram_basis(const std::vector<std::complex<double>>& gram, int n, double cutoff, int* rank)
{
    using Wide = std::complex<double>;
    *rank = 0;
    std::vector<Wide> coefficients, images, c(n), image(n);
    std::vector<double> remaining(n);
    for (int j = 0; j < n; ++j)
    {
        if (!std::isfinite(std::abs(gram[j + j * n])))
        {
            return {};
        }
        remaining[j] = std::max(0.0, gram[j + j * n].real());
    }
    const double largest = n == 0 ? 0 : *std::max_element(remaining.begin(), remaining.end());
    const double minimum = std::max(cutoff * cutoff, 1e-8 * largest);
    for (int step = 0; step < n; ++step)
    {
        const int pivot = std::max_element(remaining.begin(), remaining.end()) - remaining.begin();
        if (remaining[pivot] <= minimum)
        {
            break;
        }
        std::fill(c.begin(), c.end(), Wide(0));
        c[pivot] = 1;
        for (int pass = 0; pass < 2; ++pass)
        {
            for (int j = 0; j < *rank; ++j)
            {
                Wide overlap = 0;
                for (int i = 0; i < n; ++i)
                {
                    overlap += std::conj(images[i + j * n]) * c[i];
                }
                for (int i = 0; i < n; ++i)
                {
                    c[i] -= coefficients[i + j * n] * overlap;
                }
            }
        }
        std::fill(image.begin(), image.end(), Wide(0));
        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                image[i] += gram[i + j * n] * c[j];
            }
        }
        Wide norm2 = 0;
        for (int i = 0; i < n; ++i)
        {
            norm2 += std::conj(c[i]) * image[i];
        }
        remaining[pivot] = 0;
        if (!std::isfinite(std::abs(norm2)) || norm2.real() <= minimum)
        {
            continue;
        }
        const double norm = std::sqrt(norm2.real());
        for (int i = 0; i < n; ++i)
        {
            c[i] /= norm;
            image[i] /= norm;
            remaining[i] = std::max(0.0, remaining[i] - std::norm(image[i]));
        }
        coefficients.insert(coefficients.end(), c.begin(), c.end());
        images.insert(images.end(), image.begin(), image.end());
        ++*rank;
    }
    return coefficients;
}

bool LinearSmallLU::factor(const std::vector<std::complex<double>>& matrix, int n)
{
    size_ = 0;
    if (n <= 0 || matrix.size() != static_cast<size_t>(n) * n)
    {
        return false;
    }
    lu_ = matrix;
    pivots_.resize(n);
    double scale = 0.0;
    for (const std::complex<double>& value: matrix)
    {
        if (!std::isfinite(std::abs(value)))
        {
            return false;
        }
        scale = std::max(scale, std::abs(value));
    }
    const double cutoff = 100 * n * std::numeric_limits<double>::epsilon() * scale;
    for (int k = 0; k < n; ++k)
    {
        int pivot = k;
        for (int i = k + 1; i < n; ++i)
        {
            if (std::abs(lu_[i + k * n]) > std::abs(lu_[pivot + k * n]))
            {
                pivot = i;
            }
        }
        if (std::abs(lu_[pivot + k * n]) <= cutoff)
        {
            return false;
        }
        pivots_[k] = pivot;
        for (int j = 0; j < n; ++j)
        {
            std::swap(lu_[k + j * n], lu_[pivot + j * n]);
        }
        for (int i = k + 1; i < n; ++i)
        {
            lu_[i + k * n] /= lu_[k + k * n];
            for (int j = k + 1; j < n; ++j)
            {
                lu_[i + j * n] -= lu_[i + k * n] * lu_[k + j * n];
            }
        }
    }
    size_ = n;
    return true;
}

bool LinearSmallLU::solve(std::vector<std::complex<double>>* rhs, int columns) const
{
    if (size_ == 0 || rhs->size() != static_cast<size_t>(size_) * columns)
    {
        return false;
    }
    for (int j = 0; j < columns; ++j)
    {
        std::complex<double>* b = rhs->data() + j * size_;
        for (int k = 0; k < size_; ++k)
        {
            std::swap(b[k], b[pivots_[k]]);
        }
        for (int k = 0; k < size_; ++k)
        {
            for (int i = k + 1; i < size_; ++i)
            {
                b[i] -= lu_[i + k * size_] * b[k];
            }
        }
        for (int k = size_ - 1; k >= 0; --k)
        {
            b[k] /= lu_[k + k * size_];
            if (!std::isfinite(std::abs(b[k])))
            {
                return false;
            }
            for (int i = 0; i < k; ++i)
            {
                b[i] -= lu_[i + k * size_] * b[k];
            }
        }
    }
    return true;
}
template <typename T, typename Device>
void LinearAlgebra<T, Device>::reduce(std::vector<Wide>* values) const
{
#ifdef __MPI
    if (comm_.nproc > 1)
    {
        Parallel_Common::reduce_data(values->data(), values->size(), comm_.comm);
    }
#endif
}

template <typename T, typename Device>
std::vector<std::complex<double>> LinearAlgebra<T, Device>::projection_cross(int ld, int dim, int nx, int ny, const T* x, const T* y)
{
    if (std::is_same<T, Wide>::value)
    {
        return cross(ld, dim, nx, ny, x, y);
    }
    std::vector<T> native(nx * ny, T(0));
    if (dim > 0)
    {
        linear_buffer<T, Device>(&native_products_, nx * ny);
        const T one(1);
        const T zero(0);
        T* products = native_products_.template data<T>();
        ModuleBase::gemm_op<T, Device>()('C', 'N', nx, ny, dim, &one, x, ld, y, ld, &zero, products, nx);
        base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>()(native.data(), products, native.size());
    }
#ifdef __MPI
    if (comm_.nproc > 1)
    {
        Parallel_Common::reduce_data(native.data(), native.size(), comm_.comm);
    }
#endif
    return std::vector<Wide>(native.begin(), native.end());
}

template <typename T, typename Device>
std::vector<std::complex<double>> LinearAlgebra<T, Device>::arnoldi_dots(int ld,
                                                                         int dim,
                                                                         int bands,
                                                                         int count,
                                                                         int stride,
                                                                         const T* basis,
                                                                         const T* x)
{
    if (std::is_same<T, Wide>::value)
    {
        return dots(ld, dim, bands, count, stride, basis, x);
    }
    const int size = bands * count;
    linear_buffer<T, Device>(&native_products_, size);
    T* products = native_products_.template data<T>();
    linear_op<T, Device>().native_dots(ld, dim, bands, count, stride, basis, x, products);
    std::vector<T> native(size);
    base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>()(native.data(), products, size);
#ifdef __MPI
    if (comm_.nproc > 1)
    {
        Parallel_Common::reduce_data(native.data(), size, comm_.comm);
    }
#endif
    return std::vector<Wide>(native.begin(), native.end());
}

template <typename T, typename Device>
std::vector<std::complex<double>> LinearAlgebra<T, Device>::dots(int ld,
                                                                 int dim,
                                                                 int bands,
                                                                 int count,
                                                                 int stride,
                                                                 const T* basis,
                                                                 const T* x)
{
    const int size = bands * count;
    linear_buffer<Wide, Device>(&products_, size);
    linear_op<T, Device>().wide_dots(ld, dim, bands, count, stride, basis, x, products_.template data<Wide>());
    std::vector<Wide> result(size);
    base_device::memory::synchronize_memory_op<Wide, base_device::DEVICE_CPU, Device>()(result.data(),
                                                                                        products_.template data<Wide>(),
                                                                                        size);
    reduce(&result);
    return result;
}

template <typename T, typename Device>
std::vector<std::complex<double>> LinearAlgebra<T, Device>::cross(int ld, int dim, int nx, int ny, const T* x, const T* y)
{
    std::vector<Wide> result(nx * ny, Wide(0));
    if (dim > 0)
    {
        const Wide* a = reinterpret_cast<const Wide*>(x);
        const Wide* b = reinterpret_cast<const Wide*>(y);
        if (!std::is_same<T, Wide>::value)
        {
            linear_buffer<Wide, Device>(&left_, static_cast<int64_t>(ld) * nx);
            linear_buffer<Wide, Device>(&right_, static_cast<int64_t>(ld) * ny);
            // Copy valid rows only: padding may be uninitialized.
            for (int j = 0; j < nx; ++j)
            {
                base_device::memory::cast_memory_op<Wide, T, Device, Device>()(left_.template data<Wide>() + j * ld, x + j * ld, dim);
            }
            for (int j = 0; j < ny; ++j)
            {
                base_device::memory::cast_memory_op<Wide, T, Device, Device>()(right_.template data<Wide>() + j * ld, y + j * ld, dim);
            }
            a = left_.template data<Wide>();
            b = right_.template data<Wide>();
        }
        linear_buffer<Wide, Device>(&products_, nx * ny);
        const Wide one(1), zero(0);
        ModuleBase::gemm_op<Wide, Device>()('C', 'N', nx, ny, dim, &one, a, ld, b, ld, &zero, products_.template data<Wide>(), nx);
        base_device::memory::synchronize_memory_op<Wide, base_device::DEVICE_CPU, Device>()(result.data(),
                                                                                            products_.template data<Wide>(),
                                                                                            result.size());
    }
    reduce(&result);
    return result;
}

template <typename T, typename Device>
void LinearAlgebra<T, Device>::expand(int ld, int dim, int nx, int ny, const T* x, const std::vector<Wide>& c, T* y, T beta)
{
    if (dim == 0 || nx == 0 || ny == 0)
    {
        return;
    }
    std::vector<T> narrow(c.begin(), c.end());
    linear_buffer<T, Device>(&coefficients_, narrow.size());
    T* coeff = coefficients_.template data<T>();
    base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(coeff, narrow.data(), narrow.size());
    const T one(1);
    ModuleBase::gemm_op<T, Device>()('N', 'N', dim, ny, nx, &one, x, ld, coeff, nx, &beta, y, ld);
}

template class LinearAlgebra<std::complex<float>, base_device::DEVICE_CPU>;
template class LinearAlgebra<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class LinearAlgebra<std::complex<float>, base_device::DEVICE_GPU>;
template class LinearAlgebra<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
