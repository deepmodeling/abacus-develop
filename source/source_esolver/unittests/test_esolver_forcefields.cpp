#include "source_esolver/esolver_forcefields.h"

#include "source_base/constants.h"
#include "source_cell/mdcell.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/parameter.h"
#include "source_md/forcefield.h"

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <memory>
#include <utility>

namespace
{
class FixedStressField : public ModuleMD::ForceField
{
  public:
    void initialize(BaseCell&, const Input_para&) override {}

    void compute(BaseCell&, double& energy, ModuleBase::matrix& force, ModuleBase::matrix& stress) override
    {
        energy = 6.0;
        force.create(0, 3);
        stress.create(3, 3);
        stress.zero_out();
        stress(0, 0) = 8.0;
        stress(1, 1) = 8.0;
        stress(2, 2) = 8.0;
    }
};

void expect_external_pressure(BaseCell& cell)
{
    Input_para input;
    input.press1 = 1.0;
    input.press2 = 2.0;
    input.press3 = 3.0;
    std::ofstream running_log;
    std::unique_ptr<ModuleMD::ForceField> field(new FixedStressField());
    ModuleESolver::esolver_forcefields solver(std::move(field), running_log);
    solver.before_all_runners(cell, input);
    const int step = 0;
    solver.runner(cell, step);

    ModuleBase::matrix stress;
    solver.cal_stress(cell, stress);
    const double unit_transform = ModuleBase::RYDBERG_SI / std::pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    EXPECT_NEAR(stress(0, 0), 8.0 - input.press1 / unit_transform, 1.0e-12);
    EXPECT_NEAR(stress(1, 1), 8.0 - input.press2 / unit_transform, 1.0e-12);
    EXPECT_NEAR(stress(2, 2), 8.0 - input.press3 / unit_transform, 1.0e-12);
}

TEST(ESolverForceFieldsTest, AppliesExternalPressureToUnitCellStress)
{
    UnitCell cell;
    expect_external_pressure(cell);
}

TEST(ESolverForceFieldsTest, AppliesExternalPressureToMDCellStress)
{
    MDCell cell;
    expect_external_pressure(cell);
}
} // namespace
