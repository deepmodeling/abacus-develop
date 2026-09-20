#ifndef POTXC_H
#define POTXC_H

#include "pot_base.h"

namespace elecstate
{

class PotXC : public PotBase
{
  public:
    // constructor for exchange-correlation potential
    // meta-GGA should input matrix of kinetic potential, it is optional
    PotXC(const bool domag,
          const bool domag_z,
          const int gga_grad,
          const ModulePW::PW_Basis* rho_basis_in,
          double* etxc_in,
          double* vtxc_in,
          ModuleBase::matrix* vofk_in = nullptr)
        : vofk(vofk_in), etxc_(etxc_in), vtxc_(vtxc_in),
          domag_(domag), domag_z_(domag_z), gga_grad_(gga_grad)
    {
        this->rho_basis_ = rho_basis_in;
        this->dynamic_mode = true;
        this->fixed_mode = false;
    }

    void cal_v_eff(const Charge*const chg, const UnitCell*const ucell, ModuleBase::matrix& v_eff) override;

    ModuleBase::matrix* vofk = nullptr;
    double* etxc_ = nullptr;
    double* vtxc_ = nullptr;
  private:
    const bool domag_;
    const bool domag_z_;
    const int gga_grad_;
};

} // namespace elecstate

#endif
