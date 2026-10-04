#include "source_pw/module_pwdft/td_history_pw.h"

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"

#include <algorithm>

namespace pw
{
namespace
{

template <typename Real, typename Device>
void prepare_pair(const ModuleBase::matrix& potential, ct::Tensor* previous, ct::Tensor* midpoint)
{
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device = ct::DeviceTypeToEnum<CtDevice>::value;
    const int64_t size = std::max<int64_t>(1, static_cast<int64_t>(potential.nr) * potential.nc);
    if (previous->NumElements() < size || previous->data_type() != ct::DataTypeToEnum<Real>::value || previous->device_type() != device)
    {
        *previous = ct::Tensor(ct::DataTypeToEnum<Real>::value, device, {size});
        *midpoint = ct::Tensor(ct::DataTypeToEnum<Real>::value, device, {size});
    }
}

template <typename Real, typename Device>
const Real* average_potential(const ModuleBase::matrix& shape, const Real* current, const ct::Tensor& previous, ct::Tensor* midpoint)
{
    const int size = shape.nr * shape.nc;
    if (size == 0)
    {
        return current;
    }
    Real* buffer = midpoint->template data<Real>();
    ModuleBase::vector_add_vector_op<Real, Device, Real>()(size, buffer, previous.template data<Real>(), Real(0.5), current, Real(0.5));
    return buffer;
}

template <typename Real, typename Device>
void save_potential(const ModuleBase::matrix& shape, const Real* current, ct::Tensor* previous)
{
    const int size = shape.nr * shape.nc;
    if (size > 0)
    {
        base_device::memory::synchronize_memory_op<Real, Device, Device>()(previous->template data<Real>(), current, size);
    }
}

} // namespace

template <typename T, typename Device>
void TDHistoryPW<T, Device>::prepare(const elecstate::Potential& potential, const bool needs_ked)
{
    ModuleBase::timer::start("TDHistoryPW", "prepare");
    prepare_pair<Real, Device>(potential.get_veff_smooth(), &veff_prev_, &veff_mid_);
    if (needs_ked)
    {
        prepare_pair<Real, Device>(potential.get_vofk_smooth(), &vofk_prev_, &vofk_mid_);
    }
    ModuleBase::timer::end("TDHistoryPW", "prepare");
}

template <typename T, typename Device>
TDPotentialView<typename GetTypeReal<T>::type> TDHistoryPW<T, Device>::prop_pot(elecstate::Potential& potential,
                                                                                const bool needs_ked,
                                                                                const bool corrector)
{
    ModuleBase::timer::start("TDHistoryPW", "prop_pot");
    TDPotentialView<Real> view{potential.template get_veff_smooth_data<Real>(), potential.template get_vofk_smooth_data<Real>()};
    if (corrector)
    {
        view.veff = average_potential<Real, Device>(potential.get_veff_smooth(), view.veff, veff_prev_, &veff_mid_);
        if (needs_ked)
        {
            view.vofk = average_potential<Real, Device>(potential.get_vofk_smooth(), view.vofk, vofk_prev_, &vofk_mid_);
        }
    }
    ModuleBase::timer::end("TDHistoryPW", "prop_pot");
    return view;
}

template <typename T, typename Device>
void TDHistoryPW<T, Device>::save(const psi::Psi<T, Device>& current, elecstate::Potential& potential, const bool needs_ked)
{
    ModuleBase::timer::start("TDHistoryPW", "save");
    if (!previous_)
    {
        previous_.reset(new psi::Psi<T, Device>(current));
    }
    *previous_ = current;
    save_potential<Real, Device>(potential.get_veff_smooth(), potential.template get_veff_smooth_data<Real>(), &veff_prev_);
    if (needs_ked)
    {
        save_potential<Real, Device>(potential.get_vofk_smooth(), potential.template get_vofk_smooth_data<Real>(), &vofk_prev_);
    }
    ModuleBase::timer::end("TDHistoryPW", "save");
}

template <typename T, typename Device>
const psi::Psi<T, Device>& TDHistoryPW<T, Device>::previous() const
{
    if (!previous_)
    {
        ModuleBase::WARNING_QUIT("TDHistoryPW", "No converged state has been saved for propagation.");
    }
    return *previous_;
}

template class TDHistoryPW<std::complex<float>, base_device::DEVICE_CPU>;
template class TDHistoryPW<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class TDHistoryPW<std::complex<float>, base_device::DEVICE_GPU>;
template class TDHistoryPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace pw
