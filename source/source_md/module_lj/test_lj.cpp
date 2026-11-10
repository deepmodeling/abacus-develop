#include "source_md/module_lj/lj.h"
#include "source_esolver/esolver_potential.h"
#include "source_cell/dom_decomp.h"
#include "source_md/test/setcell.h"

#include "gtest/gtest.h"

#include <memory>
#include <fstream>
#include <utility>

TEST(LJPotential, UnitCellEnergyForceAndStress)
{
    UnitCell cell;
    Setcell::setupcell(cell);
    Input_para input;
    Setcell::parameters(input);
    const std::string output_dir = testing::TempDir();
    std::ofstream running_log;
    std::unique_ptr<ModuleMD::Potential> potential(new ModuleMD::LJ(output_dir, running_log));
    ModuleESolver::ESolver_Potential solver(std::move(potential), running_log);
    solver.before_all_runners(cell, input);
    solver.runner(cell, 0);

    EXPECT_NEAR(solver.cal_energy(), -0.023915637247068761, 1e-12);
    ModuleBase::matrix force;
    solver.cal_force(cell, force);
    EXPECT_NEAR(force(0, 0), 0.00099635466178755387, 1e-12);
    ModuleBase::matrix stress;
    solver.cal_stress(cell, stress);
    EXPECT_NEAR(stress(0, 0), 1.6072044445526372e-6, 1e-10);
}

TEST(LJPotential, MDCellWritesOwnedForces)
{
    UnitCell unitcell;
    Setcell::setupcell(unitcell);
    MDCell cell = Setcell::setup_mdcell(unitcell);
    Input_para input;
    Setcell::parameters(input);
    const std::string output_dir = testing::TempDir();
    std::ofstream running_log;
    std::unique_ptr<ModuleMD::Potential> potential(new ModuleMD::LJ(output_dir, running_log));
    ModuleESolver::ESolver_Potential solver(std::move(potential), running_log);
    solver.before_all_runners(cell, input);
    DomainDecomposition decomp;
    decomp.init(ModuleBase::world_comm_domain(), cell.latvec(), cell.lat0(), 0.0, 0.0);
    decomp.prepare_neighbors(cell);
    solver.runner(cell, 0);

    EXPECT_NEAR(solver.cal_energy(), -0.023915637247068761, 1e-12);
    ModuleBase::matrix force;
    solver.cal_force(cell, force);
    for (std::size_t i = 0; i < cell.owned_atoms().size(); ++i)
    {
        EXPECT_DOUBLE_EQ(force(i, 0), cell.owned_atoms()[i].force.x);
        EXPECT_DOUBLE_EQ(force(i, 1), cell.owned_atoms()[i].force.y);
        EXPECT_DOUBLE_EQ(force(i, 2), cell.owned_atoms()[i].force.z);
    }
}
