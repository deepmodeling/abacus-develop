#ifndef ELECSTATEPW_SDFT_H
#define ELECSTATEPW_SDFT_H
#include "elecstate_pw.h"
namespace elecstate
{
template <typename T, typename Device>
class ElecStatePW_SDFT : public ElecStatePW<T, Device>
{
  public:
    ElecStatePW_SDFT(ModulePW::PW_Basis_K* wfc_basis_in,
                     Charge* chg_in,
                     K_Vectors* pkv_in,
                     UnitCell* ucell_in,
                     pseudopot_cell_vnl* ppcell_in,
                     ModulePW::PW_Basis* rhodpw_in,
                     ModulePW::PW_Basis* rhopw_in,
                     ModulePW::PW_Basis_Big* bigpw_in)
        : ElecStatePW<T,
                      Device>(wfc_basis_in, chg_in, pkv_in, ucell_in, ppcell_in, rhodpw_in, rhopw_in, bigpw_in)
    {
        this->classname = "ElecStatePW_SDFT";
    }
    virtual void psiToRho(const psi::Psi<T, Device>& psi) override;
  private:
    using Real = typename GetTypeReal<T>::type;
    using setmem_var_op = base_device::memory::set_memory_op<Real, Device>;
    using castmem_var_d2h_op = base_device::memory::cast_memory_op<double, Real, base_device::DEVICE_CPU, Device>;
};
} // namespace elecstate
#endif