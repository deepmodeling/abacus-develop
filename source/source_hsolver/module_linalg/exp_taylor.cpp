#include "source_hsolver/module_linalg/exp_taylor.h"

#include "source_base/timer.h"
#include "source_hsolver/module_linalg/kernels/linear_op.h"

namespace hsolver
{
template <typename T, typename Device>
ExpTaylor<T, Device>::ExpTaylor(const ExpControls& exp_ctrl, const diag_comm_info& comm) : ExpAction<T, Device>(exp_ctrl, comm)
{
}

template <typename T, typename Device>
ExpResult ExpTaylor<T, Device>::apply(const HSOperator<T, Device>& op, const OrbLayout& layout, double alpha, const T* input, T* output)
{
    return evaluate(op, layout, alpha, input, output, nullptr);
}

template <typename T, typename Device>
ExpResult ExpTaylor<T, Device>::apply_pair(const HSOperator<T, Device>& op,
                                           const OrbLayout& layout,
                                           double alpha,
                                           const T* input,
                                           T* half,
                                           T* full)
{
    return evaluate(op, layout, alpha, input, half, full);
}

template <typename T, typename Device>
ExpResult ExpTaylor<T, Device>::evaluate(const HSOperator<T, Device>& op,
                                         const OrbLayout& layout,
                                         double alpha,
                                         const T* input,
                                         T* half,
                                         T* full)
{
    ModuleBase::timer::start("ExpTaylor", "evaluate");
    const int ld = layout.ld;
    const int bands = layout.bands;
    const int stride = ld * bands;
    const int dim = layout.npol == 2 ? ld : layout.npw;
    const int64_t count = static_cast<int64_t>(ld) * bands;
    linear_buffer<T, Device>(&term_, count);
    linear_buffer<T, Device>(&image_, count);
    T* term = term_.template data<T>();
    T* image = image_.template data<T>();
    base_device::memory::synchronize_memory_op<T, Device, Device>()(term, input, count);
    zero_pad<T, Device>(term, layout);
    base_device::memory::synchronize_memory_op<T, Device, Device>()(half, term, count);
    if (full)
    {
        base_device::memory::synchronize_memory_op<T, Device, Device>()(full, term, count);
    }
    ExpResult result;
    result.order_or_dim = this->exp_ctrl_.taylor_order;
    // Scale at the longer requested time; half-step weights cannot overflow at high order.
    const double series_alpha = full ? 2.0 * alpha : alpha;
    for (int k = 1; k <= this->exp_ctrl_.taylor_order; ++k)
    {
        op.hpsi(term, image, ld, bands);
        zero_pad<T, Device>(image, layout);
        const T coefficient(0, -series_alpha / k);
        linear_op<T, Device>().batch(ld, dim, bands, term, image, image, coefficient, T(0), nullptr, nullptr, nullptr);
        const T half_weight = full ? T(std::ldexp(1.0, -k)) : T(1);
        linear_op<T, Device>().batch(ld, dim, bands, half, half, term, T(1), half_weight, nullptr, nullptr, nullptr);
        if (full)
        {
            linear_op<T, Device>().batch(ld, dim, bands, full, full, term, T(1), T(1), nullptr, nullptr, nullptr);
        }
        ++result.hpsi_calls;
        result.hpsi_cols += bands;
    }
    const std::vector<std::complex<double>> norms = this->algebra_.dots(ld, dim, bands, 1, stride, half, half);
    std::vector<std::complex<double>> full_norms;
    if (full)
    {
        full_norms = this->algebra_.dots(ld, dim, bands, 1, stride, full, full);
    }
    for (int b = 0; b < bands; ++b)
    {
        if (!std::isfinite(norms[b].real()) || (full && !std::isfinite(full_norms[b].real())))
        {
            result.state = ExpState::failed;
            result.failed_band = b;
            break;
        }
    }
    ModuleBase::timer::end("ExpTaylor", "evaluate");
    return result;
}

template class ExpTaylor<std::complex<float>, base_device::DEVICE_CPU>;
template class ExpTaylor<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class ExpTaylor<std::complex<float>, base_device::DEVICE_GPU>;
template class ExpTaylor<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
