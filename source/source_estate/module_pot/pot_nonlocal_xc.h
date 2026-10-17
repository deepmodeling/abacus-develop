#ifndef ABACUS_POT_NONLOCAL_XC_H
#define ABACUS_POT_NONLOCAL_XC_H
#include "pot_base.h"

#include <memory>
#include <string>

namespace vdw
{
class NonlocalXC;
}

namespace elecstate
{
// Additive density-grid nonlocal XC correction, registered after PotXC on the
// CPU/PW SCF path. The numerical method lives behind vdw::NonlocalXC.
class PotNonlocalXC final : public PotBase
{
  public:
    PotNonlocalXC(const ModulePW::PW_Basis* density_basis,
                  double* etxc,
                  double* vtxc,
                  const std::string& method,
                  double b,
                  double c);
    ~PotNonlocalXC() override;
    void cal_v_eff(const Charge* chg, const UnitCell* cell, ModuleBase::matrix& v_eff) override;

  private:
    double* etxc_;
    double* vtxc_;
    std::unique_ptr<vdw::NonlocalXC> backend_;
};
} // namespace elecstate
#endif // ABACUS_POT_NONLOCAL_XC_H
