#include "source_md/module_snap/snap_potential.h"

#include "source_cell/dom_decomp.h"
#include "source_io/module_parameter/parameter.h"
#include "source_md/test/setcell.h"

#include "gtest/gtest.h"

#include <cmath>
#include <fstream>
#include <string>

TEST(SnapPotential, UnitCellAndMDCellAgree)
{
    const std::string source_file = __FILE__;
    const std::string marker = "/source/source_md/";
    const std::size_t marker_pos = source_file.find(marker);
    ASSERT_NE(marker_pos, std::string::npos);
    const std::string model_dir = source_file.substr(0, marker_pos) + "/tests/SIMULATION_MODELS/";
    const std::string coefficient_file = model_dir + "snap_102.coeff";
    const std::string parameter_file = model_dir + "snap_102.param";

    UnitCell unitcell;
    Setcell::setupcell(unitcell);
    unitcell.atoms[0].label = "Si";
    Input_para input;
    Setcell::parameters(input);
    std::ofstream running_log;
    ModuleMD::SnapPotential unit_potential(coefficient_file, parameter_file, running_log, 0);
    unit_potential.initialize(unitcell, input);
    double unit_energy = 0.0;
    ModuleBase::matrix unit_force;
    ModuleBase::matrix unit_stress(3, 3);
    unit_potential.compute(unitcell, unit_energy, unit_force, unit_stress);

    MDCell mdcell = Setcell::setup_mdcell(unitcell);
    ModuleMD::SnapPotential md_potential(coefficient_file, parameter_file, running_log, 0);
    md_potential.initialize(mdcell, input);
    DomainDecomposition decomp;
    decomp.init(ModuleBase::world_comm_domain(), mdcell.latvec(), mdcell.lat0(), 0.0, 0.0);
    decomp.prepare_neighbors(mdcell);
    double md_energy = 0.0;
    ModuleBase::matrix unused_force;
    ModuleBase::matrix md_stress(3, 3);
    md_potential.compute(mdcell, md_energy, unused_force, md_stress);
    decomp.accumulate_ghost_forces(mdcell);

    EXPECT_TRUE(std::isfinite(unit_energy));
    EXPECT_NEAR(md_energy, unit_energy, 1e-10);
    for (int direction = 0; direction < 3; ++direction)
    {
        for (int column = 0; column < 3; ++column)
        {
            EXPECT_NEAR(md_stress(direction, column), unit_stress(direction, column), 1e-10);
        }
    }
    for (std::size_t i = 0; i < mdcell.owned_atoms().size(); ++i)
    {
        EXPECT_NEAR(mdcell.owned_atoms()[i].force.x, unit_force(i, 0), 1e-10);
        EXPECT_NEAR(mdcell.owned_atoms()[i].force.y, unit_force(i, 1), 1e-10);
        EXPECT_NEAR(mdcell.owned_atoms()[i].force.z, unit_force(i, 2), 1e-10);
    }
}
