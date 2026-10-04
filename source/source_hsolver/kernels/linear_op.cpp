#include "source_hsolver/kernels/linear_op.h"

#include <complex>
#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace hsolver
{
namespace
{
#ifdef _OPENMP
bool parallel_work(const int dim, const int nvec)
{
    return static_cast<long long>(dim) * nvec >= 32768 && !omp_in_parallel() && omp_get_max_threads() > 1;
}
#endif
} // namespace

template <typename T, typename Device>
void linear_op<T, Device>::dot(const int ld, const int dim, const int nvec, const T* x, const T* y, T* out) const
{
    for (int band = 0; band < nvec; ++band)
    {
        T sum = T(0);
        for (int i = 0; i < dim; ++i)
        {
            sum += std::conj(x[band * ld + i]) * y[band * ld + i];
        }
        out[band] = sum;
    }
}

template <typename T, typename Device>
void linear_op<T, Device>::batch(int ld,
                                 int dim,
                                 int nvec,
                                 T* out,
                                 const T* x,
                                 const T* y,
                                 T a,
                                 T b,
                                 const T* ca,
                                 const T* cb,
                                 const int* skip) const
{
#pragma omp parallel for schedule(static) if (parallel_work(dim, nvec))
    for (int band = 0; band < nvec; ++band)
    {
        if (skip && skip[band])
        {
            continue;
        }
        const T aa = ca ? a * ca[band] : a;
        const T bb = cb ? b * cb[band] : b;
        for (int i = 0; i < dim; ++i)
        {
            const int j = band * ld + i;
            out[j] = aa * x[j] + bb * y[j];
        }
    }
}
template <typename T, typename Device>
void linear_op<T, Device>::gmres_update(int ld,
                                        int dim,
                                        int nvec,
                                        T* solution,
                                        T* residual,
                                        const T* direction,
                                        const T* image,
                                        const T* coefficients) const
{
#pragma omp parallel for schedule(static) if (parallel_work(dim, nvec))
    for (int band = 0; band < nvec; ++band)
    {
        const T coefficient = coefficients[band];
        if (coefficient == T(0))
        {
            continue;
        }
        const std::int64_t offset = static_cast<std::int64_t>(band) * ld;
        for (int i = 0; i < dim; ++i)
        {
            solution[offset + i] += coefficient * direction[offset + i];
            if (residual)
            {
                residual[offset + i] -= coefficient * image[offset + i];
            }
        }
    }
}

template <typename T, typename Device>
void linear_op<T, Device>::dots(int ld, int dim, int nvec, const T* x, const T* y, const T* z, const T* w, T* out, T* partial, int tiles)
    const
{
    // Independent tiles also expose parallelism when there are few bands.
#pragma omp parallel for collapse(2) schedule(static) if (parallel_work(dim, nvec))
    for (int band = 0; band < nvec; ++band)
    {
        for (int tile = 0; tile < tiles; ++tile)
        {
            T first = T(0);
            T second = T(0);
            const int begin = static_cast<long long>(dim) * tile / tiles;
            const int end = static_cast<long long>(dim) * (tile + 1) / tiles;
            for (int i = begin; i < end; ++i)
            {
                const int j = band * ld + i;
                first += std::conj(x[j]) * y[j];
                if (z)
                {
                    second += std::conj(z[j]) * w[j];
                }
            }
            partial[band * tiles + tile] = first;
            if (z)
            {
                partial[(nvec + band) * tiles + tile] = second;
            }
        }
    }
    const int count = (z ? 2 : 1) * nvec;
    for (int band = 0; band < count; ++band)
    {
        T sum = T(0);
        for (int tile = 0; tile < tiles; ++tile)
        {
            sum += partial[band * tiles + tile];
        }
        out[band] = sum;
    }
}
template <typename T, typename Device>
void linear_op<T, Device>::gather(int ld, int dim, int nvec, int slots, int stride, const T* in, T* out, const int* map) const
{
    for (int slot = 0; slot < slots; ++slot)
    {
        const std::int64_t slot_offset = static_cast<std::int64_t>(slot) * stride;
        for (int band = 0; band < nvec; ++band)
        {
            const std::int64_t source_offset = slot_offset + static_cast<std::int64_t>(map[band]) * ld;
            const std::int64_t destination_offset = slot_offset + static_cast<std::int64_t>(band) * ld;
            for (int i = 0; i < dim; ++i)
            {
                out[destination_offset + i] = in[source_offset + i];
            }
        }
    }
}

template <typename T, typename Device>
void linear_op<T, Device>::swaps(int ld, int dim, int slots, int stride, int count, T* vectors, const int* pairs) const
{
    for (int slot = 0; slot < slots; ++slot)
    {
        const std::int64_t slot_offset = static_cast<std::int64_t>(slot) * stride;
        for (int i = 0; i < dim; ++i)
        {
            for (int k = 0; k < count; ++k)
            {
                const std::int64_t first = slot_offset + static_cast<std::int64_t>(pairs[2 * k]) * ld + i;
                const std::int64_t second = slot_offset + static_cast<std::int64_t>(pairs[2 * k + 1]) * ld + i;
                std::swap(vectors[first], vectors[second]);
            }
        }
    }
}

template <typename T, typename Device>
void linear_op<T, Device>::bicg_update(int ld, int dim, int nvec, bool direction, T* out, const T* x, const T* y, const T* a, const T* b)
    const
{
#pragma omp parallel for collapse(2) schedule(static) if (parallel_work(dim, nvec))
    for (int band = 0; band < nvec; ++band)
    {
        for (int i = 0; i < dim; ++i)
        {
            const int j = band * ld + i;
            out[j] = direction ? x[j] + a[band] * (out[j] - b[band] * y[j]) : (out[j] + a[band] * x[j]) + b[band] * y[j];
        }
    }
}

template <typename T, typename Device>
void linear_op<T, Device>::diagonal(int ld, int dim, int nvec, const T* inverse, const T* x, T* y) const
{
#pragma omp parallel for collapse(2) schedule(static) if (parallel_work(dim, nvec))
    for (int band = 0; band < nvec; ++band)
    {
        for (int i = 0; i < dim; ++i)
        {
            y[band * ld + i] = inverse[i] * x[band * ld + i];
        }
    }
}

template <typename T, typename Device>
void linear_op<T, Device>::cgs_direction(int ld, int dim, int nvec, T* p, T* u, const T* r, const T* q, const T* beta) const
{
#pragma omp parallel for collapse(2) schedule(static) if (parallel_work(dim, nvec))
    for (int band = 0; band < nvec; ++band)
    {
        for (int i = 0; i < dim; ++i)
        {
            const int j = band * ld + i;
            u[j] = r[j] + beta[band] * q[j];
            p[j] = u[j] + beta[band] * (q[j] + beta[band] * p[j]);
        }
    }
}

template <typename T, typename Device>
void linear_op<T, Device>::cgs_alpha(int ld, int dim, int nvec, T* q, T* direction, const T* u, const T* v, const T* alpha) const
{
#pragma omp parallel for collapse(2) schedule(static) if (parallel_work(dim, nvec))
    for (int band = 0; band < nvec; ++band)
    {
        for (int i = 0; i < dim; ++i)
        {
            const int j = band * ld + i;
            q[j] = u[j] - alpha[band] * v[j];
            direction[j] = u[j] + q[j];
        }
    }
}

template <typename T, typename Device>
void linear_op<T, Device>::cgs_finish(int ld, int dim, int nvec, T* x, T* r, const T* direction, const T* ad, const T* alpha) const
{
#pragma omp parallel for collapse(2) schedule(static) if (parallel_work(dim, nvec))
    for (int band = 0; band < nvec; ++band)
    {
        for (int i = 0; i < dim; ++i)
        {
            const int j = band * ld + i;
            x[j] += alpha[band] * direction[j];
            r[j] -= alpha[band] * ad[j];
        }
    }
}
template <typename T, typename Device>
void linear_op<T, Device>::wide_dots(int ld,
                                     int dim,
                                     int nvec,
                                     int count,
                                     int stride,
                                     const T* basis,
                                     const T* x,
                                     std::complex<double>* out) const
{
#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (parallel_work(dim, nvec * count))
#endif
    for (int j = 0; j < count; ++j)
    {
        for (int b = 0; b < nvec; ++b)
        {
            const std::int64_t band_offset = static_cast<std::int64_t>(b) * ld;
            const std::int64_t basis_offset = static_cast<std::int64_t>(j) * stride + band_offset;
            std::complex<double> sum = 0.0;
            for (int i = 0; i < dim; ++i)
            {
                sum += std::conj(std::complex<double>(basis[basis_offset + i])) * std::complex<double>(x[band_offset + i]);
            }
            out[j * nvec + b] = sum;
        }
    }
}
template <typename T, typename Device>
void linear_op<T, Device>::native_dots(int ld, int dim, int nvec, int count, int stride, const T* basis, const T* x, T* out) const
{
#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (parallel_work(dim, nvec * count))
#endif
    for (int j = 0; j < count; ++j)
    {
        for (int band = 0; band < nvec; ++band)
        {
            const std::int64_t band_offset = static_cast<std::int64_t>(band) * ld;
            const std::int64_t basis_offset = static_cast<std::int64_t>(j) * stride + band_offset;
            T sum(0);
            for (int i = 0; i < dim; ++i)
            {
                sum += std::conj(basis[basis_offset + i]) * x[band_offset + i];
            }
            out[j * nvec + band] = sum;
        }
    }
}

template <typename T, typename Device>
void linear_op<T, Device>::synchronize() const
{
}
template class linear_op<std::complex<float>, base_device::DEVICE_CPU>;
template class linear_op<std::complex<double>, base_device::DEVICE_CPU>;

} // namespace hsolver
