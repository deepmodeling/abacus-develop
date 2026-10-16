#include "source_hsolver/cn_subspace.h"

#include "source_hsolver/kernels/linear_op.h"
#include "source_hsolver/linear_workspace.h"

namespace hsolver
{
template <typename T, typename Device>
bool CNSubspace<T, Device>::prepare(LinearAlgebra<T, Device>& algebra, int ld, int dim, int bands, const T* u, const T* b)
{
    const int64_t elements = static_cast<int64_t>(ld) * bands;
    linear_buffer<T, Device>(&image_, elements);
    linear_buffer<T, Device>(&seed_, elements);
    linear_buffer<T, Device>(&residual_, elements);
    linear_op<T, Device>().batch(ld, dim, bands, image(), u, b, T(2), T(-1), nullptr, nullptr, nullptr);
    const std::vector<std::complex<double>> projected = algebra.projection_cross(ld, dim, bands, bands, u, image());
    std::vector<std::complex<double>> c = algebra.projection_cross(ld, dim, bands, bands, u, b);
    if (!factor_.factor(projected, bands) || !factor_.solve(&c, bands))
    {
        return false;
    }
    for (const std::complex<double>& value: c)
    {
        if (!std::isfinite(std::abs(T(value))))
        {
            return false;
        }
    }
    algebra.expand(ld, dim, bands, bands, u, c, seed(), T(0));
    algebra.expand(ld, dim, bands, bands, image(), c, residual(), T(0));
    linear_op<T, Device>().batch(ld, dim, bands, residual(), b, residual(), T(1), T(-1), nullptr, nullptr, nullptr);
    const int stride = ld * bands;
    const std::vector<std::complex<double>> norms = algebra.dots(ld, dim, bands, 1, stride, residual(), residual());
    for (const std::complex<double>& value: norms)
    {
        if (!std::isfinite(linear_norm(value)))
        {
            return false;
        }
    }
    return true;
}

template class CNSubspace<std::complex<float>, base_device::DEVICE_CPU>;
template class CNSubspace<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class CNSubspace<std::complex<float>, base_device::DEVICE_GPU>;
template class CNSubspace<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
