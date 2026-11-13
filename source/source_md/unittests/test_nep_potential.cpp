#include "source_md/module_nep/nep_potential.h"

#include "source_cell/dom_decomp.h"
#include "source_md/test/setcell.h"

#include "gtest/gtest.h"

#include <cmath>
#include <fstream>
#include <string>

TEST(NEPPotential, UnitCellCalculationOrMissingBackend)
{
    const std::string source_file = __FILE__;
    const std::string marker = "/source/source_md/";
    const std::size_t marker_pos = source_file.find(marker);
    ASSERT_NE(marker_pos, std::string::npos);
    const std::string model_file = source_file.substr(0, marker_pos) + "/tests/PP_ORB/nep_hfo2.txt";

    UnitCell cell;
    Setcell::setupcell(cell);
    cell.atoms[0].label = "Hf";
    Input_para input;
    Setcell::parameters(input);
    const std::string output_dir = testing::TempDir();
    std::ofstream running_log;
    ModuleMD::NEPPotential potential(model_file, output_dir, running_log);
    potential.initialize(cell, input);

    double energy = 0.0;
    ModuleBase::matrix force;
    ModuleBase::matrix stress(3, 3);
#ifdef __NEP
    potential.compute(cell, energy, force, stress);
    EXPECT_TRUE(std::isfinite(energy));
    for (int i = 0; i < cell.nat; ++i)
    {
        for (int direction = 0; direction < 3; ++direction)
        {
            EXPECT_TRUE(std::isfinite(force(i, direction)));
        }
    }

    MDCell mdcell = Setcell::setup_mdcell(cell);
    ModuleMD::NEPPotential md_potential(model_file, output_dir, running_log);
    md_potential.initialize(mdcell, input);
    DomainDecomposition decomp;
    decomp.init(ModuleBase::world_comm_domain(), mdcell.latvec(), mdcell.lat0(), 0.0, 0.0);
    decomp.prepare_neighbors(mdcell);
    double md_energy = 0.0;
    ModuleBase::matrix unused_force;
    ModuleBase::matrix md_stress(3, 3);
    md_potential.compute(mdcell, md_energy, unused_force, md_stress);
    decomp.accumulate_ghost_forces(mdcell);
    EXPECT_NEAR(md_energy, energy, 1e-8);
    for (int row = 0; row < 3; ++row)
    {
        for (int column = 0; column < 3; ++column)
        {
            EXPECT_NEAR(md_stress(row, column), stress(row, column), 1e-8);
        }
    }
    for (std::size_t i = 0; i < mdcell.owned_atoms().size(); ++i)
    {
        EXPECT_NEAR(mdcell.owned_atoms()[i].force.x, force(i, 0), 1e-8);
        EXPECT_NEAR(mdcell.owned_atoms()[i].force.y, force(i, 1), 1e-8);
        EXPECT_NEAR(mdcell.owned_atoms()[i].force.z, force(i, 2), 1e-8);
    }
#else
    EXPECT_EXIT(potential.compute(cell, energy, force, stress), ::testing::ExitedWithCode(1), "");
#endif
}
