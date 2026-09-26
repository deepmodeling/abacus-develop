#ifndef ESOLVER_DFTB_NATIVE_H
#define ESOLVER_DFTB_NATIVE_H

#include "esolver.h"
#include "source_dftb/periodic_scc.h"
#include "source_dftb/skf_data.h"

#include <vector>

namespace ModuleESolver
{
class ESolver_DFTBNative : public ESolver
{
  public:
    ESolver_DFTBNative() { classname = "ESolver_DFTBNative"; }

    void before_all_runners(BaseCell& cell, const Input_para& inp) override;
    void runner(BaseCell& cell, const int istep) override;
    void after_all_runners(BaseCell& cell) override;
    double cal_energy() override;
    void cal_force(BaseCell& cell, ModuleBase::matrix& force) override;
    void cal_stress(BaseCell& cell, ModuleBase::matrix& stress) override;

  private:
    void load_model(const UnitCell& ucell, const Input_para& inp);
    ModuleDFTB::DftbPeriodicInput make_geometry(const UnitCell& ucell) const;

    std::vector<ModuleDFTB::SkfData> skfiles_;
    ModuleDFTB::DftbPeriodicInput template_;
    ModuleDFTB::DftbPeriodicResult result_;
    double energy_ry_ = 0.0;
    int output_precision_ = 12;
};
} // namespace ModuleESolver

#endif // ESOLVER_DFTB_NATIVE_H
