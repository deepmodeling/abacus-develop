#include "source_estate/module_pot/sccs_fixed_src.h"
#include "source_estate/module_pot/pot_sccs.h"
#include "source_cell/cell_tools.h"
#include "source_hamilt/module_sccs/sccs_ion_charge.h"
#include "source_hamilt/module_sccs/test/sccs_test.h"
#include <sstream>
#include <string>

namespace
{
const elecstate::SccsResume no_resume = elecstate::SccsResume();

// The value after a key in an sccs_debug report.
std::string report_value(const std::string& report, const std::string& key)
{
    std::istringstream tokens(report);
    std::string token;
    while (tokens >> token)
    {
        if (token == key && tokens >> token)
        {
            return token;
        }
    }
    return std::string();
}
}

class SccsFixedSourcesTest : public SccsTest::PwTest
{
protected:
    void SetUp() override
    {
        SccsTest::PwTest::SetUp();
        cell.lat0 = length;
        cell.tpiba = tpiba;
        cell.omega = basis.omega;
        unitcell::AtomData atom;
        atom.position = {2.5, 2.5, 2.5};
        atom.mass = 1.0;
        atom.valence_charge = 1.0;
        atom.atomic_number = 1;
        atoms.push_back(atom);
        atom.position.x = 4.0;
        atom.mass = 16.0;
        atom.valence_charge = 6.0;
        atom.atomic_number = 8;
        atoms.push_back(atom);
        config.core_electrons = true;
        config.core_spreads = {0.5, 0.7};
    }
    void rebuild_and_compare()
    {
        EXPECT_FALSE(cache.update(cell, basis, atoms, config));
        std::vector<double> expected;
        ModuleSccs::gauss_ion_rho(atoms, basis, cell.tpiba,
                                  ModuleSccs::gaussian_ion_spread, expected);
        EXPECT_EQ(cache.ionic_density(), expected);
        ModuleSccs::gauss_core_rho(atoms, basis, cell.tpiba, config.core_spreads, expected);
        EXPECT_EQ(cache.core_density(), expected);
        EXPECT_TRUE(cache.update(cell, basis, atoms, config));
    }
    UnitCell cell;
    std::vector<unitcell::AtomData> atoms;
    ModuleSccs::SccsConfig config;
    elecstate::SccsFixedSources cache;
};

TEST_F(SccsFixedSourcesTest, ReusesSourcesUntilTheAtomsChange)
{
    EXPECT_TRUE(cache.ionic_density().empty());
    rebuild_and_compare();
    const double* stored_density = cache.ionic_density().data();
    EXPECT_TRUE(cache.update(cell, basis, atoms, config));
    EXPECT_EQ(cache.ionic_density().data(), stored_density);
    EXPECT_EQ(cache.positions().size(), static_cast<std::size_t>(basis.nrxx));
    atoms[1].mass = 17.0; // Not an input to the Gaussians.
    EXPECT_TRUE(cache.update(cell, basis, atoms, config));
    atoms[0].position.x += 0.1;
    rebuild_and_compare();
    atoms[1].valence_charge = 5.0;
    rebuild_and_compare();
    atoms[0].atomic_number = 2; // Changes which atoms a single sccs_corespread skips (Z == zv).
    rebuild_and_compare();
}

TEST_F(SccsFixedSourcesTest, CachedSourcesMatchFreshFullCavityEnergyPotentialAndForce)
{
    Atom atom;
    cell.ntype = 1;
    cell.nat = 1;
    cell.atoms = &atom;
    atom.na = 1;
    atom.mass = 1.0;
    atom.ncpp.zv = 1.0;
    atom.tau = {ModuleBase::Vector3<double>(0.25, 0.25, 0.25)};
    const double density_value = 0.6 / basis.omega;
    std::vector<double> density(basis.nrxx, density_value);
    double* channel = density.data();
    Charge charge;
    charge.nspin = 1;
    charge.rho = &channel;
    config = ModuleSccs::preset_config(ModuleSccs::Preset::WaterNeutral);
    config.boundary = ModuleSccs::Boundary::Pcc0d;
    config.cavity.epsilon_bulk = 5.0;
    config.core_electrons = true;
    config.core_spreads = {0.5};
    config.collect_report = true;
    ModuleSccs::PolarizationSolverParameters solver;
    solver.tolerance_rms = 1e-12;
    solver.tolerance_max = 1e-12;
    elecstate::PotSccs cached(&basis, config, solver, 0.6, no_resume);
    for (int step = 0; step < 3; ++step)
    {
        if (step == 1) { for (double& value : density) { value *= 0.95; } }
        if (step == 2)
        {
            atom.tau[0].x += 0.01;
        }
        elecstate::PotSccs fresh(&basis, config, solver, 0.6, no_resume);
        ModuleBase::matrix first(1, basis.nrxx);
        ModuleBase::matrix second(1, basis.nrxx);
        cached.cal_v_eff(&charge, &cell, first);
        fresh.cal_v_eff(&charge, &cell, second);
        EXPECT_NEAR(cached.get_energy(), fresh.get_energy(), 1e-11);
        for (int i = 0; i < basis.nrxx; ++i)
        {
            EXPECT_NEAR(first(0, i), second(0, i), 1e-9);
        }
        ModuleBase::matrix force1(1, 3);
        ModuleBase::matrix force2(1, 3);
        cached.add_sol_force(cell, force1);
        fresh.add_sol_force(cell, force2);
        for (int i = 0; i < 3; ++i)
        {
            EXPECT_NEAR(force1(0, i), force2(0, i), 1e-9);
        }
        std::ostringstream output;
        cached.write_iter(output, 2, 1e-5, 0.0);
        const std::string reused = step == 1 ? "1" : "0";
        EXPECT_EQ(report_value(output.str(), "CACHED_SOURCES"), reused);
    }
}

TEST_F(SccsFixedSourcesTest, NextIonicStepResumesTheWarmStart)
{
    Atom atom;
    cell.ntype = 1;
    cell.nat = 1;
    cell.atoms = &atom;
    atom.na = 1;
    atom.mass = 1.0;
    atom.ncpp.zv = 1.0;
    atom.tau = {ModuleBase::Vector3<double>(0.25, 0.25, 0.25)};
    const double density_value = 0.6 / basis.omega;
    std::vector<double> density(basis.nrxx, density_value);
    double* channel = density.data();
    Charge charge;
    charge.nspin = 1;
    charge.rho = &channel;
    config = ModuleSccs::preset_config(ModuleSccs::Preset::WaterNeutral);
    config.boundary = ModuleSccs::Boundary::Pcc0d;
    config.cavity.epsilon_bulk = 5.0;
    config.core_electrons = true;
    config.core_spreads = {0.5};
    ModuleSccs::PolarizationSolverParameters solver;
    solver.tolerance_rms = 1e-12;
    solver.tolerance_max = 1e-12;
    elecstate::PotSccs first_step(&basis, config, solver, 0.6, no_resume);
    ModuleBase::matrix first_potential(1, basis.nrxx);
    first_step.cal_v_eff(&charge, &cell, first_potential);
    const elecstate::SccsResume resume = first_step.resume_state();
    EXPECT_TRUE(resume.active);
    EXPECT_EQ(resume.restart_potential.size(), static_cast<std::size_t>(basis.nrxx));

    // The rebuilt component of the next ionic step starts from the previous
    // solution; moved atoms do not discard it.
    atom.tau[0].x += 0.002;
    elecstate::PotSccs resumed(&basis, config, solver, 0.6, resume);
    elecstate::PotSccs cold(&basis, config, solver, 0.6, no_resume);
    ModuleBase::matrix resumed_potential(1, basis.nrxx);
    ModuleBase::matrix cold_potential(1, basis.nrxx);
    resumed.cal_v_eff(&charge, &cell, resumed_potential);
    cold.cal_v_eff(&charge, &cell, cold_potential);
    EXPECT_NEAR(resumed.get_energy(), cold.get_energy(), 1e-10);
    for (int i = 0; i < basis.nrxx; ++i)
    {
        EXPECT_NEAR(resumed_potential(0, i), cold_potential(0, i), 1e-8);
    }
    std::ostringstream resumed_output;
    std::ostringstream cold_output;
    resumed.write_iter(resumed_output, 2, 1e-5, 0.0);
    cold.write_iter(cold_output, 2, 1e-5, 0.0);
    EXPECT_EQ(report_value(resumed_output.str(), "WARM_START"), "1");
    EXPECT_EQ(report_value(cold_output.str(), "WARM_START"), "0");
}
