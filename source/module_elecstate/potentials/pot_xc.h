#ifndef POTXC_H
#define POTXC_H

#include "module_hamilt_general/module_xc/xc_functional.h"
#include "pot_base.h"

namespace elecstate
{

class PotXC : public PotBase
{
  public:
    // constructor for exchange-correlation potential
    // meta-GGA should input matrix of kinetic potential, it is optional
    PotXC(const ModulePW::PW_Basis* rho_basis_in,
          double* etxc_in,
          double* vtxc_in,
          ModuleBase::matrix* vofk_in = nullptr)
        : etxc_(etxc_in), vtxc_(vtxc_in), vofk(vofk_in)
    {
        this->rho_basis_ = rho_basis_in;
        this->dynamic_mode = true;
        this->fixed_mode = false;
    }

    void cal_v_eff(const Charge* chg, const UnitCell* ucell, ModuleBase::matrix& v_eff) override;

    ModuleBase::matrix* vofk = nullptr;
    double* etxc_ = nullptr;
    double* vtxc_ = nullptr;
};

} // namespace elecstate

#endif