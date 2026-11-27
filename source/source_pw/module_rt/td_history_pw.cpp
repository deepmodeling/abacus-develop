#include "source_pw/module_rt/td_history_pw.h"

#include "source_base/module_device/memory_op.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"

#include <algorithm>

namespace pw
{
namespace
{

template <typename Real, typename Device>
void prepare_potential(const ModuleBase::matrix& potential, ct::Tensor* previous)
{
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device = ct::DeviceTypeToEnum<CtDevice>::value;
    const int64_t size = std::max<int64_t>(1, static_cast<int64_t>(potential.nr) * potential.nc);
    if (previous->NumElements() < size || previous->data_type() != ct::DataTypeToEnum<Real>::value || previous->device_type() != device)
    {
        *previous = ct::Tensor(ct::DataTypeToEnum<Real>::value, device, {size});
    }
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

TDElecSnap capture_td(const elecstate::fenergy& energy, const ModuleBase::matrix& eigenvalues, const Charge& charge)
{
    TDElecSnap state;
    state.energy = energy;
    state.eigenvalues = eigenvalues;
    state.density.resize(charge.nspin * charge.nrxx);
    for (int spin = 0; spin < charge.nspin; ++spin)
    {
        std::copy(charge.rho[spin], charge.rho[spin] + charge.nrxx, state.density.begin() + spin * charge.nrxx);
    }
    return state;
}

void restore_td(const TDElecSnap& state, elecstate::fenergy& energy, ModuleBase::matrix& eigenvalues, Charge& charge)
{
    ModuleBase::timer::start("pw", "restore_td");
    energy = state.energy;
    eigenvalues = state.eigenvalues;
    for (int spin = 0; spin < charge.nspin; ++spin)
    {
        const double* density = state.density.data() + spin * charge.nrxx;
        std::copy(density, density + charge.nrxx, charge.rho[spin]);
    }
    charge.save_rho_before_sum_band();
    ModuleBase::timer::end("pw", "restore_td");
}

template <typename T, typename Device>
void TDHistoryPW<T, Device>::prepare(const elecstate::Potential& potential)
{
    ModuleBase::timer::start("TDHistoryPW", "prepare");
    prepare_potential<Real, Device>(potential.get_veff_smooth(), &veff_prev_);
    ModuleBase::timer::end("TDHistoryPW", "prepare");
}

template <typename T, typename Device>
const typename GetTypeReal<T>::type* TDHistoryPW<T, Device>::left_potential() const
{
    return veff_prev_.template data<Real>();
}

template <typename T, typename Device>
void TDHistoryPW<T, Device>::save(const psi::Psi<T, Device>& current, elecstate::Potential& potential)
{
    ModuleBase::timer::start("TDHistoryPW", "save");
    if (!previous_)
    {
        previous_.reset(new psi::Psi<T, Device>(current));
    }
    *previous_ = current;
    ++generation_;
    save_potential<Real, Device>(potential.get_veff_smooth(), potential.template get_veff_smooth_data<Real>(), &veff_prev_);
    ModuleBase::timer::end("TDHistoryPW", "save");
}

template <typename T, typename Device>
const psi::Psi<T, Device>& TDHistoryPW<T, Device>::previous() const
{
    if (!previous_)
    {
        ModuleBase::WARNING_QUIT("TDHistoryPW", "No endpoint state has been saved for propagation.");
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
