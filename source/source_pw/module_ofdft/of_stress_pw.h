#ifndef OF_STRESS_PW_H
#define OF_STRESS_PW_H

#include "source_estate/elecstate.h"
#include "source_pw/module_pwdft/vl_pw.h"
#include "source_pw/module_pwdft/stress_func.h"

namespace vdw
{
struct VdwResult;
}

class OF_Stress_PW : public Stress_Func<double>
{
  public:
    OF_Stress_PW(const elecstate::ElecState* pelec_in, ModulePW::PW_Basis* rhopw_in)
        : pelec(pelec_in), rhopw(rhopw_in){};

    // calculate the stress in OFDFT
    void cal_stress(ModuleBase::matrix& sigmatot,
                    ModuleBase::matrix& kinetic_stress,
                    UnitCell& ucell,
                    const vdw::VdwResult* vdw_result,
                    ModuleSymmetry::Symmetry* p_symm,
                    const pseudopot_cell_vl& locpp,
                    Structure_Factor* p_sf,
                    K_Vectors* p_kv);

  protected:
    const elecstate::ElecState* pelec = nullptr;
    ModulePW::PW_Basis* rhopw = nullptr;
};
#endif
