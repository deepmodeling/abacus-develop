#include "ppcg_kernel_op.h"
#include <algorithm>

namespace hsolver
{
namespace
{
#ifdef _OPENMP
enum { openmp_work_threshold = 4096 };
#endif
}

template <typename T, typename Device>
void PPCGVectorOps<T, Device>::column_products(
    const int dim, const int ld, const int count, const T* x, const T* y, double* host_products)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (dim * count > openmp_work_threshold)
#endif
    for (int j = 0; j < count; ++j)
    {
        double sum = 0.0;
        for (int i = 0; i < dim; ++i)
        {
            const int offset = i + j * ld;
            sum += double(std::real(std::conj(x[offset]) * y[offset]));
        }
        host_products[j] = sum;
    }
}

template <typename T, typename Device>
void PPCGVectorOps<T, Device>::residual(
    const int dim, const int ld, const int count, const T* hx, const T* sx, const Real* eigenvalues, T* output)
{
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) if (dim * count > openmp_work_threshold)
#endif
    for (int j = 0; j < count; ++j)
    {
        for (int i = 0; i < dim; ++i)
        {
            const int offset = i + j * ld;
            output[offset] = hx[offset] - sx[offset] * eigenvalues[j];
        }
    }
}

template <typename T, typename Device>
void PPCGVectorOps<T, Device>::precondition(
    const int dim, const int ld, const int count, const Real* diagonal, T* values)
{
    const Real minimum_preconditioner = Real(1.0e-12);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) if (dim * count > openmp_work_threshold)
#endif
    for (int j = 0; j < count; ++j)
    {
        for (int i = 0; i < dim; ++i)
        {
            values[i + j * ld] /= std::max(diagonal[i], minimum_preconditioner);
        }
    }
}

template <typename T, typename Device>
void PPCGVectorOps<T, Device>::scale_columns(
    const int dim, const int ld, const int count, const double* scales, T* x, T* sx, T* hx)
{
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) if (dim * count > openmp_work_threshold)
#endif
    for (int j = 0; j < count; ++j)
    {
        for (int i = 0; i < dim; ++i)
        {
            const int offset = i + j * ld;
            const Real scale = Real(scales[j]);
            x[offset] *= scale;
            sx[offset] *= scale;
            hx[offset] *= scale;
        }
    }
}

template struct PPCGVectorOps<std::complex<float>, base_device::DEVICE_CPU>;
template struct PPCGVectorOps<std::complex<double>, base_device::DEVICE_CPU>;
} // namespace hsolver
