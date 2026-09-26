#ifndef ABACUS_POT_RVV10_H
#define ABACUS_POT_RVV10_H
#include "pot_base.h"
#include "source_hamilt/module_xc/xc_rvv10_pw.h"

namespace elecstate
{
// Independent rVV10 nonlocal correction for the CPU/PW SCF path.
class PotRvv10 final : public PotBase
{
  public:
    PotRvv10(const ModulePW::PW_Basis* density_basis,
             double* etxc,
             double* vtxc,
             double b = Rvv10::rvv10_b_default,
             double c = Rvv10::rvv10_c_default);
    void cal_v_eff(const Charge* chg, const UnitCell* cell, ModuleBase::matrix& v_eff) override;

  private:
    double* etxc_;
    double* vtxc_;
    Rvv10::Evaluator evaluator_;
};
} // namespace elecstate
#endif
