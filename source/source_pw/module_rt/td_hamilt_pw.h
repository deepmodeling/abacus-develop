#ifndef PW_TD_HAMILT_PW_H
#define PW_TD_HAMILT_PW_H
#include "source_hamilt/hamilt_hs_adapter.h"
#include "source_pw/module_pwdft/hamilt_pw.h"

namespace pw
{
/** @brief A call-scoped PW sample; every k update refreshes shifted projectors. */
template <typename T, typename Device>
class TDHamiltPW final : public hsolver::HSOperator<T, Device>
{
  private:
    using Real = typename T::value_type;
    hamilt::HamiltPW<T, Device>& hamilt_;
    hamilt::HamiltHSOperator<T, Device> op_;
    const Real* veff_ = nullptr;
    ModuleBase::Vector3<double> vecpot_ha_;

  public:
    TDHamiltPW(hamilt::HamiltPW<T, Device>& hamilt, const ModulePW::PW_Basis_K& basis);
    /** @brief Select a potential sample whose storage outlives its applications. */
    void bind(const Real* veff, const ModuleBase::Vector3<double>& vecpot_ha);
    void update_k(int ik) override;
    void hpsi(const T* x, T* hx, int ld, int bands) const override;
    void spsi(const T* x, T* sx, int ld, int bands) const override;
};
} // namespace pw
#endif
