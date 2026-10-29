#include "source_pw/module_rt/td_hamilt_pw.h"

namespace pw
{
template <typename T, typename Device>
TDHamiltPW<T, Device>::TDHamiltPW(hamilt::HamiltPW<T, Device>& hamilt, const ModulePW::PW_Basis_K& basis)
    : hamilt_(hamilt), op_(&hamilt, &basis)
{
}
template <typename T, typename Device>
void TDHamiltPW<T, Device>::bind(const Real* veff, const ModuleBase::Vector3<double>& vecpot_ha)
{
    veff_ = veff;
    vecpot_ha_ = vecpot_ha;
    hamilt_.bind_td_state(veff_, vecpot_ha_);
}
template <typename T, typename Device>
void TDHamiltPW<T, Device>::update_k(int ik)
{
    hamilt_.bind_td_state(veff_, vecpot_ha_);
    op_.update_k(ik);
}
template <typename T, typename Device>
void TDHamiltPW<T, Device>::hpsi(const T* x, T* hx, int ld, int bands) const
{
    op_.hpsi(x, hx, ld, bands);
}
template <typename T, typename Device>
void TDHamiltPW<T, Device>::spsi(const T* x, T* sx, int ld, int bands) const
{
    op_.spsi(x, sx, ld, bands);
}
template class TDHamiltPW<std::complex<float>, base_device::DEVICE_CPU>;
template class TDHamiltPW<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class TDHamiltPW<std::complex<float>, base_device::DEVICE_GPU>;
template class TDHamiltPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace pw
