#include "source_hsolver/linear_workspace.h"

#include "source_base/module_device/memory_op.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_hsolver/kernels/linear_op.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace hsolver
{

LinearSolveTimer::LinearSolveTimer(const char* solver) : solver_(solver)
{
    ModuleBase::timer::start(solver_, "solve");
}

LinearSolveTimer::~LinearSolveTimer()
{
    ModuleBase::timer::end(solver_, "solve");
}

template <typename T, typename Device>
LinearWorkspace<T, Device>::LinearWorkspace(const diag_comm_info& comm, int slots) : slots_(slots), comm_(comm)
{
    if (slots_ <= 0)
    {
        throw std::invalid_argument("Linear workspace requires a positive slot count.");
    }
    // A default Tensor is a one-element float CPU tensor, not an empty buffer.
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device = ct::DeviceTypeToEnum<CtDevice>::value;
    dots_ = ct::Tensor(ct::DataTypeToEnum<T>::value, device, {0});
    partial_ = ct::Tensor(ct::DataTypeToEnum<T>::value, device, {0});
    coefficients_ = ct::Tensor(ct::DataTypeToEnum<T>::value, device, {0});
    mask_ = ct::Tensor(ct::DataTypeToEnum<int>::value, device, {0});
    permutation_ = ct::Tensor(ct::DataTypeToEnum<int>::value, device, {0});
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::prepare(const int ld,
                                         const int dim,
                                         const int nvec,
                                         const T* x,
                                         const T* b,
                                         const double tolerance,
                                         const int max_iter)
{
    if (dim < 0 || ld < dim || nvec <= 0 || !std::isfinite(tolerance) || tolerance <= 0.0 || max_iter < 0
        || ld > std::numeric_limits<int>::max() / nvec || (ld > 0 && (x == nullptr || b == nullptr)) || (ld > 0 && x == b)
        || comm_.nproc < 1)
    {
        throw std::invalid_argument("Invalid linear-solver layout, buffers, or convergence controls.");
    }
    const int vector_elements = ld * nvec;
    const int size = std::max(1, vector_elements);
    if (size > capacity_)
    {
        using CtDevice = typename ct::PsiToContainer<Device>::type;
        const int64_t slot_elements = size;
        vectors_ = ct::Tensor(ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<CtDevice>::value, {slots_, slot_elements});
        capacity_ = size;
    }
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::clear()
{
    const size_t elements = static_cast<size_t>(slots_) * capacity_;
    base_device::memory::set_memory_op<T, Device>()(vectors_.template data<T>(), 0, elements);
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::reset_statistics()
{
    operator_calls_ = 0;
    operator_columns_ = 0;
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::statistics(LinearSolveResult* result) const
{
    result->operator_calls = operator_calls_;
    result->operator_columns = operator_columns_;
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::apply(const LinearOperator<T, Device>& op, const T* x, T* y, int ld, int nvec)
{
    ModuleBase::timer::start("LinearWorkspace", "apply");
    ++operator_calls_;
    operator_columns_ += nvec;
    op.apply(x, y, ld, nvec);
    ModuleBase::timer::end("LinearWorkspace", "apply");
}

template <typename T, typename Device>
T* LinearWorkspace<T, Device>::data(const int slot)
{
    return vectors_.data<T>() + static_cast<size_t>(slot) * capacity_;
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::copy(const int ld, const int dim, const int nvec, const T* x, T* y) const
{
    if (dim > 0)
    {
        base_device::memory::synchronize_memory_2d_op<T, Device, Device>()(y, ld, x, ld, dim, nvec);
    }
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::dot(const int ld, const int dim, const int nvec, const T* x, const T* y, T* out)
{

    dot_pair(ld, dim, nvec, x, y, nullptr, nullptr, out, nullptr);
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::residual(const LinearOperator<T, Device>& op,
                                          const int ld,
                                          const int dim,
                                          const int nvec,
                                          const T* x,
                                          const T* b,
                                          T* r)
{
    apply(op, x, r, ld, nvec);
    linear_op<T, Device>().batch(ld, dim, nvec, r, b, r, T(1), T(-1), nullptr, nullptr, nullptr);
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::verify(const LinearOperator<T, Device>& op,
                                        const int ld,
                                        const int dim,
                                        const int nvec,
                                        const T* x,
                                        const T* b,
                                        const std::vector<double>& threshold,
                                        T* scratch,
                                        LinearSolveResult* result)
{
    residual(op, ld, dim, nvec, x, b, scratch);
    ++result->true_checks;
    result->operator_calls = operator_calls_;
    result->operator_columns = operator_columns_;
    std::vector<T> norms(nvec);
    dot(ld, dim, nvec, scratch, scratch, norms.data());
    result->max_residual = 0.0;
    result->failed_band = -1;
    bool converged = true;
    for (int band = 0; band < nvec; ++band)
    {
        const double error = linear_norm(norms[band]);
        result->max_residual = std::isfinite(error) ? std::max(result->max_residual, error) : std::numeric_limits<double>::infinity();
        if (!std::isfinite(error) || !std::isfinite(threshold[band]) || error > threshold[band])
        {
            if (result->failed_band < 0)
            {
                result->failed_band = band;
            }
            converged = false;
        }
    }
    if (converged)
    {
        result->status = LinearSolveStatus::converged;
        result->failed_band = -1;
    }
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::dot_pair(int ld, int dim, int nvec, const T* x, const T* y, const T* z, const T* w, T* first, T* second)
{
    ModuleBase::timer::start("LinearWorkspace", "dot_pair");
    if (nvec == 0)
    {
        ModuleBase::timer::end("LinearWorkspace", "dot_pair");
        return;
    }
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device = ct::DeviceTypeToEnum<CtDevice>::value;
    const int count = (z ? 2 : 1) * nvec;
    const int requested_tiles = (dim + 2047) / 2048;
    const int capped_tiles = std::min(32, requested_tiles);
    const int tiles = std::max(1, capped_tiles);
    if (dots_.NumElements() < count)
    {
        dots_ = ct::Tensor(ct::DataTypeToEnum<T>::value, device, {count});
    }
    const int64_t partial_elements = static_cast<int64_t>(count) * tiles;
    if (partial_.NumElements() < partial_elements)
    {
        partial_ = ct::Tensor(ct::DataTypeToEnum<T>::value, device, {partial_elements});
    }
    host_dots_.resize(count);
    linear_op<T, Device>().dots(ld, dim, nvec, x, y, z, w, dots_.template data<T>(), partial_.template data<T>(), tiles);
    base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>()(host_dots_.data(), dots_.template data<T>(), count);
#ifdef __MPI
    if (comm_.nproc > 1)
    {
        Parallel_Common::reduce_data(host_dots_.data(), count, comm_.comm);
    }
#endif
    const typename std::vector<T>::const_iterator second_begin = host_dots_.cbegin() + nvec;
    std::copy(host_dots_.cbegin(), second_begin, first);
    if (z)
    {
        std::copy(second_begin, host_dots_.cend(), second);
    }
    ModuleBase::timer::end("LinearWorkspace", "dot_pair");
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::batch(int ld,
                                       int dim,
                                       int nvec,
                                       T* out,
                                       const T* x,
                                       const T* y,
                                       T a,
                                       T b,
                                       const T* ca,
                                       const T* cb,
                                       const int* skip)
{
    ModuleBase::timer::start("LinearWorkspace", "batch");
    if (dim == 0 || nvec == 0)
    {
        ModuleBase::timer::end("LinearWorkspace", "batch");
        return;
    }
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device = ct::DeviceTypeToEnum<CtDevice>::value;
    const int64_t coefficient_elements = 2 * static_cast<int64_t>(nvec);
    if (coefficients_.NumElements() < coefficient_elements)
    {
        coefficients_ = ct::Tensor(ct::DataTypeToEnum<T>::value, device, {coefficient_elements});
    }
    T* coeff = coefficients_.template data<T>();
    if (ca)
    {
        base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(coeff, ca, nvec);
    }
    if (cb)
    {
        T* second_coefficients = coeff + nvec;
        base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(second_coefficients, cb, nvec);
    }
    const int* mask = nullptr;
    if (skip)
    {
        if (mask_.NumElements() < nvec)
        {
            mask_ = ct::Tensor(ct::DataTypeToEnum<int>::value, device, {nvec});
        }
        base_device::memory::synchronize_memory_op<int, Device, base_device::DEVICE_CPU>()(mask_.template data<int>(), skip, nvec);
        mask = mask_.template data<int>();
    }
    const T* first_coefficients = ca ? coeff : nullptr;
    const T* second_coefficients = cb ? coeff + nvec : nullptr;
    linear_op<T, Device>().batch(ld, dim, nvec, out, x, y, a, b, first_coefficients, second_coefficients, mask);
    ModuleBase::timer::end("LinearWorkspace", "batch");
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::gmres_update(int ld,
                                              int dim,
                                              int nvec,
                                              T* solution,
                                              T* residual,
                                              const T* direction,
                                              const T* image,
                                              const T* coefficients,
                                              bool has_zero_coefficients)
{
    ModuleBase::timer::start("LinearWorkspace", "gmres_update");
    if (dim == 0 || nvec == 0)
    {
        ModuleBase::timer::end("LinearWorkspace", "gmres_update");
        return;
    }
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device = ct::DeviceTypeToEnum<CtDevice>::value;
    if (coefficients_.NumElements() < nvec)
    {
        coefficients_ = ct::Tensor(ct::DataTypeToEnum<T>::value, device, {nvec});
    }
    T* coeff = coefficients_.template data<T>();
    base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(coeff, coefficients, nvec);
    if (std::is_same<Device, base_device::DEVICE_GPU>::value && residual == nullptr && !has_zero_coefficients)
    {
        // Keep the existing GPU path when neither zero-column skipping nor a paired update is needed.
        linear_op<T, Device>().batch(ld, dim, nvec, solution, solution, direction, T(1), T(1), nullptr, coeff, nullptr);
    }
    else
    {
        linear_op<T, Device>().gmres_update(ld, dim, nvec, solution, residual, direction, image, coeff);
    }
    ModuleBase::timer::end("LinearWorkspace", "gmres_update");
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::swap_columns(int ld, int dim, const std::vector<int>& pairs)
{
    swap_vectors(ld, dim, slots_, capacity_, vectors_.template data<T>(), pairs);
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::swap_vectors(int ld, int dim, int slots, int stride, T* vectors, const std::vector<int>& pairs)
{
    if (pairs.empty() || dim == 0)
    {
        return;
    }
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const int64_t permutation_elements = pairs.size();
    if (permutation_.NumElements() < permutation_elements)
    {
        permutation_ = ct::Tensor(ct::DataTypeToEnum<int>::value, ct::DeviceTypeToEnum<CtDevice>::value, {permutation_elements});
    }
    int* map = permutation_.template data<int>();
    base_device::memory::synchronize_memory_op<int, Device, base_device::DEVICE_CPU>()(map, pairs.data(), pairs.size());
    const int pair_count = pairs.size() / 2;
    linear_op<T, Device>().swaps(ld, dim, slots, stride, pair_count, vectors, map);
}
template <typename T, typename Device>
void LinearWorkspace<T, Device>::restore(int ld, int dim, int nvec, const std::vector<int>& order, const T* source, T* destination)
{
    if (dim == 0)
    {
        return;
    }
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    if (permutation_.NumElements() < nvec)
    {
        permutation_ = ct::Tensor(ct::DataTypeToEnum<int>::value, ct::DeviceTypeToEnum<CtDevice>::value, {nvec});
    }
    std::vector<int> inverse(nvec);
    for (int band = 0; band < nvec; ++band)
    {
        inverse[order[band]] = band;
    }
    int* map = permutation_.template data<int>();
    base_device::memory::synchronize_memory_op<int, Device, base_device::DEVICE_CPU>()(map, inverse.data(), nvec);
    const int stride = ld * nvec;
    linear_op<T, Device>().gather(ld, dim, nvec, 1, stride, source, destination, map);
}

template <typename T, typename Device>
void LinearWorkspace<T,
                     Device>::bicg_update(int ld, int dim, int nvec, bool direction, T* out, const T* x, const T* y, const T* a, const T* b)
{
    ModuleBase::timer::start("LinearWorkspace", "bicg_update");
    if (dim == 0 || nvec == 0)
    {
        ModuleBase::timer::end("LinearWorkspace", "bicg_update");
        return;
    }
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const int64_t coefficient_elements = 2 * static_cast<int64_t>(nvec);
    if (coefficients_.NumElements() < coefficient_elements)
    {
        coefficients_ = ct::Tensor(ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<CtDevice>::value, {coefficient_elements});
    }
    host_coefficients_.resize(coefficient_elements);
    const T* a_end = a + nvec;
    const T* b_end = b + nvec;
    const typename std::vector<T>::iterator second_begin = host_coefficients_.begin() + nvec;
    std::copy(a, a_end, host_coefficients_.begin());
    std::copy(b, b_end, second_begin);
    T* coeff = coefficients_.template data<T>();
    base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(coeff,
                                                                                     host_coefficients_.data(),
                                                                                     coefficient_elements);
    const T* second_coefficients = coeff + nvec;
    linear_op<T, Device>().bicg_update(ld, dim, nvec, direction, out, x, y, coeff, second_coefficients);
    ModuleBase::timer::end("LinearWorkspace", "bicg_update");
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::cgs_direction(int ld, int dim, int nvec, T* p, T* u, const T* r, const T* q, const T* beta)
{
    ModuleBase::timer::start("LinearWorkspace", "cgs_direction");
    if (dim == 0 || nvec == 0)
    {
        ModuleBase::timer::end("LinearWorkspace", "cgs_direction");
        return;
    }
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    if (coefficients_.NumElements() < nvec)
    {
        coefficients_ = ct::Tensor(ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<CtDevice>::value, {nvec});
    }
    T* coeff = coefficients_.template data<T>();
    base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(coeff, beta, nvec);
    linear_op<T, Device>().cgs_direction(ld, dim, nvec, p, u, r, q, coeff);
    ModuleBase::timer::end("LinearWorkspace", "cgs_direction");
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::cgs_alpha(int ld, int dim, int nvec, T* q, T* direction, const T* u, const T* v, const T* alpha)
{
    ModuleBase::timer::start("LinearWorkspace", "cgs_alpha");
    if (dim == 0 || nvec == 0)
    {
        ModuleBase::timer::end("LinearWorkspace", "cgs_alpha");
        return;
    }
    T* coeff = coefficients_.template data<T>();
    base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(coeff, alpha, nvec);
    linear_op<T, Device>().cgs_alpha(ld, dim, nvec, q, direction, u, v, coeff);
    ModuleBase::timer::end("LinearWorkspace", "cgs_alpha");
}

template <typename T, typename Device>
void LinearWorkspace<T, Device>::cgs_finish(int ld, int dim, int nvec, T* x, T* r, const T* direction, const T* ad)
{
    ModuleBase::timer::start("LinearWorkspace", "cgs_finish");
    if (dim == 0 || nvec == 0)
    {
        ModuleBase::timer::end("LinearWorkspace", "cgs_finish");
        return;
    }
    linear_op<T, Device>().cgs_finish(ld, dim, nvec, x, r, direction, ad, coefficients_.template data<T>());
    ModuleBase::timer::end("LinearWorkspace", "cgs_finish");
}
template class LinearWorkspace<std::complex<float>, base_device::DEVICE_CPU>;
template class LinearWorkspace<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class LinearWorkspace<std::complex<float>, base_device::DEVICE_GPU>;
template class LinearWorkspace<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hsolver
