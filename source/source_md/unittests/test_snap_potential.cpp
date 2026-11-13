#include "source_md/module_snap/snap_potential.h"

#include "source_cell/dom_decomp.h"
#include "source_io/module_parameter/parameter.h"
#include "source_md/test/setcell.h"

#include "gtest/gtest.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

TEST(SnapPotential, UnitCellAndMDCellAgree)
{
    const std::string coefficient_file = testing::TempDir() + "snap_potential_coeff.txt";
    const std::string parameter_file = testing::TempDir() + "snap_potential_param.txt";
    {
        std::ofstream output(coefficient_file.c_str());
        ASSERT_TRUE(output.is_open());
        output << "1 56\nSi 0.5 1\n0.125\n";
        for (int index = 2; index <= 56; ++index)
        {
            const double coefficient = (index % 2 == 0 ? -1.0 : 1.0) * index * 0.001;
            output << coefficient << '\n';
        }
        ASSERT_TRUE(output.good());
    }
    {
        std::ofstream output(parameter_file.c_str());
        ASSERT_TRUE(output.is_open());
        output << "rcutfac 4.9\ntwojmax 8\nrfac0 0.99363\nrmin0 0\nquadraticflag 0\nbzeroflag 0\n";
        ASSERT_TRUE(output.good());
    }

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
    std::remove(coefficient_file.c_str());
    std::remove(parameter_file.c_str());
}
