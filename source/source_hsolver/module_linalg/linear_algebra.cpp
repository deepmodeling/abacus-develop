#include "source_hsolver/module_linalg/linear_algebra.h"

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_base/memory_recorder.h"
#include "source_hsolver/module_linalg/kernels/linear_op.h"

#include <cmath>
#include <limits>
#include <type_traits>

namespace hsolver
{
void linear_complete_hermitian_upper(int n, std::complex<double>* matrix)
{
    for (int j = 0; j < n; ++j)
    {
        matrix[j + j * n].imag(0.0);
        for (int i = j + 1; i < n; ++i)
        {
            matrix[i + j * n] = std::conj(matrix[j + i * n]);
        }
    }
}

std::vector<std::complex<double>> linear_gram_basis(const std::vector<std::complex<double>>& gram, int n, double cutoff, int* rank)
{
    using Wide = std::complex<double>;
    *rank = 0;
    std::vector<Wide> coefficients;
    std::vector<Wide> images;
    std::vector<Wide> c(n);
    std::vector<Wide> image(n);
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
    const double absolute_cutoff = cutoff * cutoff;
    const double relative_cutoff = 1e-8 * largest;
    const double minimum = std::max(absolute_cutoff, relative_cutoff);
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
            const double remainder = remaining[i] - std::norm(image[i]);
            remaining[i] = std::max(0.0, remainder);
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
        const double magnitude = std::abs(value);
        if (!std::isfinite(magnitude))
        {
            return false;
        }
        scale = std::max(scale, magnitude);
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
void LinearAlgebra<T, Device>::reduce(Wide* values, const int64_t count) const
{
#ifdef __MPI
    if (comm_.nproc > 1)
    {
        Parallel_Common::reduce_data(values, count, comm_.comm);
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
    const int64_t elements = static_cast<int64_t>(nx) * ny;
    std::vector<T> native(elements, T(0));
    if (dim > 0)
    {
        linear_buffer<T, Device>(&native_products_, elements);
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
ct::Tensor LinearAlgebra<T, Device>::dots_tensor(int ld, int dim, int bands, int count,
                                               int stride, const T* basis, const T* x)
{
    ModuleBase::timer::Scope function_timer("LinearAlgebra", "dots_tensor");
    const int size = bands * count;
    ct::Tensor result(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {size});
    linear_buffer<Wide, Device>(&products_, size);
    linear_op<T, Device>().wide_dots(ld, dim, bands, count, stride, basis, x, products_.template data<Wide>());
    base_device::memory::synchronize_memory_op<Wide, base_device::DEVICE_CPU, Device>()(
        result.data<Wide>(), products_.template data<Wide>(), size);
    reduce(result.data<Wide>(), size);
    return result;
}

template <typename T, typename Device>
std::vector<std::complex<double>> LinearAlgebra<T, Device>::cross(int ld, int dim, int nx, int ny, const T* x, const T* y)
{
    const int64_t elements = static_cast<int64_t>(nx) * ny;
    std::vector<Wide> result(elements, Wide(0));
    if (dim > 0)
    {
        const Wide* a = reinterpret_cast<const Wide*>(x);
        const Wide* b = reinterpret_cast<const Wide*>(y);
        if (!std::is_same<T, Wide>::value)
        {
            const int64_t left_elements = static_cast<int64_t>(ld) * nx;
            const int64_t right_elements = static_cast<int64_t>(ld) * ny;
            linear_buffer<Wide, Device>(&left_, left_elements);
            linear_buffer<Wide, Device>(&right_, right_elements);
            // Copy valid rows only: padding may be uninitialized.
            for (int j = 0; j < nx; ++j)
            {
                const int64_t offset = static_cast<int64_t>(j) * ld;
                Wide* destination = left_.template data<Wide>() + offset;
                const T* source = x + offset;
                base_device::memory::cast_memory_op<Wide, T, Device, Device>()(destination, source, dim);
            }
            for (int j = 0; j < ny; ++j)
            {
                const int64_t offset = static_cast<int64_t>(j) * ld;
                Wide* destination = right_.template data<Wide>() + offset;
                const T* source = y + offset;
                base_device::memory::cast_memory_op<Wide, T, Device, Device>()(destination, source, dim);
            }
            a = left_.template data<Wide>();
            b = right_.template data<Wide>();
        }
        linear_buffer<Wide, Device>(&products_, elements);
        const Wide one(1);
        const Wide zero(0);
        ModuleBase::gemm_op<Wide, Device>()('C', 'N', nx, ny, dim, &one, a, ld, b, ld, &zero, products_.template data<Wide>(), nx);
        base_device::memory::synchronize_memory_op<Wide, base_device::DEVICE_CPU, Device>()(result.data(),
                                                                                            products_.template data<Wide>(),
                                                                                            result.size());
    }
    reduce(&result);
    return result;
}

template <typename T, typename Device>
ct::Tensor LinearAlgebra<T, Device>::cross_tensor(int ld, int dim, int nx, int ny, const T* x, const T* y)
{
    ModuleBase::timer::Scope function_timer("LinearAlgebra", "cross_tensor");
    const int64_t count = static_cast<int64_t>(nx) * ny;
    ct::Tensor result(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {count});
    cross_local_into(ld, dim, nx, ny, x, y, result.data<Wide>());
    reduce(result.data<Wide>(), count);
    return result;
}

template <typename T, typename Device>
ct::Tensor LinearAlgebra<T, Device>::cross_root_tensor(int ld, int dim, int nx, int ny, const T* x, const T* y)
{
    ModuleBase::timer::Scope function_timer("LinearAlgebra", "cross_root_tensor");
    const int64_t elements = static_cast<int64_t>(nx) * ny;
    if (elements > std::numeric_limits<int>::max())
    {
        throw std::length_error("Root projection exceeds the MPI count range.");
    }
    ct::Tensor result(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {elements});
    // Root's result is the K/G owner, registered by the caller. Only non-root
    // ranks discard their local contribution after the reduction.
    const std::size_t local_bytes = comm_.rank == 0 ? 0 : result.buffer().GetAllocatedBytes();
    ModuleBase::Memory::record("LOBPCG::projection_local", local_bytes, false);
    cross_local_into(ld, dim, nx, ny, x, y, result.data<Wide>());
#ifdef __MPI
    if (comm_.nproc > 1)
    {
        const int count = elements;
        Parallel_Common::reduce_data_to_root(result.data<Wide>(), count, comm_.comm, 0);
    }
#endif
    if (comm_.rank != 0)
    {
        return ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {0});
    }
    return result;
}

template <typename T, typename Device>
void LinearAlgebra<T, Device>::cross_local_into(int ld, int dim, int nx, int ny, const T* x, const T* y, Wide* result)
{
    const int64_t elements = static_cast<int64_t>(nx) * ny;
    if (dim == 0)
    {
        std::fill_n(result, elements, Wide(0));
    }
    if (dim > 0)
    {
        const Wide* a = reinterpret_cast<const Wide*>(x);
        const Wide* b = reinterpret_cast<const Wide*>(y);
        if (!std::is_same<T, Wide>::value)
        {
            const int64_t left_elements = static_cast<int64_t>(ld) * nx;
            const int64_t right_elements = static_cast<int64_t>(ld) * ny;
            linear_buffer<Wide, Device>(&left_, left_elements);
            linear_buffer<Wide, Device>(&right_, right_elements);
            // Copy valid rows only: padding may be uninitialized.
            for (int j = 0; j < nx; ++j)
            {
                const int64_t offset = static_cast<int64_t>(j) * ld;
                Wide* destination = left_.template data<Wide>() + offset;
                const T* source = x + offset;
                base_device::memory::cast_memory_op<Wide, T, Device, Device>()(destination, source, dim);
            }
            for (int j = 0; j < ny; ++j)
            {
                const int64_t offset = static_cast<int64_t>(j) * ld;
                Wide* destination = right_.template data<Wide>() + offset;
                const T* source = y + offset;
                base_device::memory::cast_memory_op<Wide, T, Device, Device>()(destination, source, dim);
            }
            a = left_.template data<Wide>();
            b = right_.template data<Wide>();
        }
        Wide* products = result;
        if (!std::is_same<Device, base_device::DEVICE_CPU>::value)
        {
            linear_buffer<Wide, Device>(&products_, elements);
            products = products_.template data<Wide>();
        }
        const Wide one(1);
        const Wide zero(0);
        ModuleBase::gemm_op<Wide, Device>()('C', 'N', nx, ny, dim, &one, a, ld, b, ld, &zero, products, nx);
        if (!std::is_same<Device, base_device::DEVICE_CPU>::value)
        {
            base_device::memory::synchronize_memory_op<Wide, base_device::DEVICE_CPU, Device>()(
                result, products, elements);
        }
    }
}

template <typename T, typename Device>
std::vector<std::complex<double>> LinearAlgebra<T, Device>::gram(int ld, int dim, int bands, const T* input)
{
    ModuleBase::timer::start("LinearAlgebra", "gram");
    const int64_t elements = static_cast<int64_t>(bands) * bands;
    std::vector<Wide> result(elements, Wide(0));
    if (dim > 0 && bands > 0)
    {
        const Wide* wide = reinterpret_cast<const Wide*>(input);
        if (!std::is_same<T, Wide>::value)
        {
            const int64_t input_elements = static_cast<int64_t>(ld) * bands;
            linear_buffer<Wide, Device>(&left_, input_elements);
            Wide* converted = left_.template data<Wide>();
            // The existing conversion kernels use int indices, even though their wrappers accept size_t.
            if (ld == dim && input_elements <= std::numeric_limits<int>::max())
            {
                base_device::memory::cast_memory_op<Wide, T, Device, Device>()(converted, input, input_elements);
            }
            else
            {
                // Padding may be uninitialized; convert only the active rows of each column.
                for (int band = 0; band < bands; ++band)
                {
                    const int64_t offset = static_cast<int64_t>(band) * ld;
                    Wide* destination = converted + offset;
                    const T* source = input + offset;
                    base_device::memory::cast_memory_op<Wide, T, Device, Device>()(destination, source, dim);
                }
            }
            wide = converted;
        }
        linear_buffer<Wide, Device>(&products_, elements);
        const Wide one(1);
        const Wide zero(0);
        Wide* products = products_.template data<Wide>();
        ModuleBase::gemm_op<Wide, Device>()('C', 'N', bands, bands, dim, &one, wide, ld, wide, ld, &zero, products, bands);
        base_device::memory::synchronize_memory_op<Wide, base_device::DEVICE_CPU, Device>()(result.data(), products, elements);
    }
    reduce(&result);
    ModuleBase::timer::end("LinearAlgebra", "gram");
    return result;
}

template <typename T, typename Device>
ct::Tensor LinearAlgebra<T, Device>::gram_tensor(int ld, int dim, int bands, const T* input)
{
    ModuleBase::timer::Scope function_timer("LinearAlgebra", "gram_tensor");
    const int64_t count = static_cast<int64_t>(bands) * bands;
    ct::Tensor result(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {count});
    gram_local_into(ld, dim, bands, input, result.data<Wide>());
    reduce(result.data<Wide>(), count);
    return result;
}

template <typename T, typename Device>
ct::Tensor LinearAlgebra<T, Device>::gram_root_tensor(int ld, int dim, int bands, const T* input)
{
    ModuleBase::timer::Scope function_timer("LinearAlgebra", "gram_root_tensor");
    const int64_t elements = static_cast<int64_t>(bands) * bands;
    if (elements > std::numeric_limits<int>::max())
    {
        throw std::length_error("Root Gram exceeds the MPI count range.");
    }
    ct::Tensor result(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {elements});
    const std::size_t local_bytes = comm_.rank == 0 ? 0 : result.buffer().GetAllocatedBytes();
    ModuleBase::Memory::record("LOBPCG::projection_local", local_bytes, false);
    gram_local_into(ld, dim, bands, input, result.data<Wide>());
#ifdef __MPI
    if (comm_.nproc > 1)
    {
        const int count = elements;
        Parallel_Common::reduce_data_to_root(result.data<Wide>(), count, comm_.comm, 0);
    }
#endif
    if (comm_.rank != 0)
    {
        return ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {0});
    }
    return result;
}

template <typename T, typename Device>
void LinearAlgebra<T, Device>::gram_local_into(int ld, int dim, int bands, const T* input, Wide* result)
{
    const bool use_herk = std::is_same<Device, base_device::DEVICE_CPU>::value;
    const int64_t elements = static_cast<int64_t>(bands) * bands;
    if (dim == 0)
    {
        std::fill_n(result, elements, Wide(0));
    }
    if (dim > 0 && bands > 0)
    {
        const Wide* wide = reinterpret_cast<const Wide*>(input);
        if (!std::is_same<T, Wide>::value)
        {
            const int64_t input_elements = static_cast<int64_t>(ld) * bands;
            linear_buffer<Wide, Device>(&left_, input_elements);
            Wide* converted = left_.template data<Wide>();
            // Conversion kernels use int indices, although their wrappers accept size_t.
            if (ld == dim && input_elements <= std::numeric_limits<int>::max())
            {
                base_device::memory::cast_memory_op<Wide, T, Device, Device>()(converted, input, input_elements);
            }
            else
            {
                // Padding may be uninitialized; convert only the active rows of each column.
                for (int band = 0; band < bands; ++band)
                {
                    const int64_t offset = static_cast<int64_t>(band) * ld;
                    Wide* destination = converted + offset;
                    const T* source = input + offset;
                    base_device::memory::cast_memory_op<Wide, T, Device, Device>()(destination, source, dim);
                }
            }
            wide = converted;
        }
        Wide* products = result;
        if (!std::is_same<Device, base_device::DEVICE_CPU>::value)
        {
            linear_buffer<Wide, Device>(&products_, elements);
            products = products_.template data<Wide>();
        }
        if (use_herk)
        {
            ModuleBase::gram_herk_cpu(ld, dim, bands, wide, products);
        }
        else
        {
            const Wide one(1);
            const Wide zero(0);
            ModuleBase::gemm_op<Wide, Device>()('C', 'N', bands, bands, dim, &one,
                                               wide, ld, wide, ld, &zero, products, bands);
        }
        if (!std::is_same<Device, base_device::DEVICE_CPU>::value)
        {
            base_device::memory::synchronize_memory_op<Wide, base_device::DEVICE_CPU, Device>()(
                result, products, elements);
        }
    }
    // HERK leaves the lower triangle untouched. Complete it before reduction;
    // empty local partitions retain their initialized zero contribution.
    if (use_herk)
    {
        linear_complete_hermitian_upper(bands, result);
    }
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

template <typename T, typename Device>
const T* LinearAlgebra<T, Device>::prepare_coefficients(const Wide* coefficients, const int64_t count)
{
    ModuleBase::timer::Scope function_timer("LinearAlgebra", "prepare_coefficients");
    if (count == 0)
    {
        return nullptr;
    }
    linear_buffer<T, Device>(&coefficients_, count);
    T* coeff = coefficients_.template data<T>();
    if (std::is_same<T, Wide>::value)
    {
        const T* input = reinterpret_cast<const T*>(coefficients);
        base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(coeff, input, count);
    }
    else
    {
        const ct::DataType type = ct::DataTypeToEnum<T>::value;
        ct::Tensor narrow(type, ct::DeviceType::CpuDevice, {count});
        std::copy_n(coefficients, count, narrow.data<T>());
        base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(
            coeff, narrow.data<T>(), count);
    }
    return coeff;
}

template <typename T, typename Device>
void LinearAlgebra<T, Device>::expand_prepared(int ld, int dim, int nx, int ny, const T* x, const T* c, T* y, T beta)
{
    ModuleBase::timer::Scope function_timer("LinearAlgebra", "expand_prepared");
    if (dim == 0 || nx == 0 || ny == 0)
    {
        return;
    }
    const T one(1);
    ModuleBase::gemm_op<T, Device>()('N', 'N', dim, ny, nx, &one, x, ld, c, nx, &beta, y, ld);
}

template <typename T, typename Device>
std::size_t LinearAlgebra<T, Device>::storage_bytes(const ct::DeviceType device) const
{
    const ct::Tensor* buffers[] = {&products_, &coefficients_, &left_, &right_, &native_products_};
    std::size_t bytes = 0;
    for (const ct::Tensor* buffer : buffers)
    {
        if (buffer->device_type() == device)
        {
            bytes += buffer->buffer().GetAllocatedBytes();
        }
    }
    return bytes;
}

template class LinearAlgebra<std::complex<float>, base_device::DEVICE_CPU>;
template class LinearAlgebra<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class LinearAlgebra<std::complex<float>, base_device::DEVICE_GPU>;
template class LinearAlgebra<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
