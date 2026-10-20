#include "source_estate/module_pot/pot_sccs.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_hamilt/module_sccs/test/sccs_test.h"
#include "source_base/global_variable.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

namespace
{
const elecstate::SccsResume no_resume = elecstate::SccsResume();
}

// Public storage fixtures: no INPUT initialization or privately owned atom maps.
UnitCell::UnitCell() {}
UnitCell::~UnitCell() {}
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}
SepPot::SepPot() {}
SepPot::~SepPot() {}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}
Charge::Charge() {}
Charge::~Charge() {}

class PotSccsTest : public SccsTest::PwTest
{
protected:
    // Two atoms in a uniform electron density: the solvation force is the
    // finite difference of the Rydberg energy and is added once to an existing
    // force.
    void check_solvation_force(const ModuleSccs::SccsConfig& config,
                               const ModuleSccs::PolarizationSolverParameters& solver)
    {
        UnitCell cell;
        Atom atom;
        cell.lat0 = length;
        cell.tpiba = tpiba;
        cell.omega = basis.omega;
        cell.ntype = 1;
        cell.nat = 2;
        cell.atoms = &atom;
        atom.na = 2;
        atom.ncpp.zv = 1.0;
        atom.tau = {ModuleBase::Vector3<double>(0.21, 0.32, 0.43),
                    ModuleBase::Vector3<double>(0.64, 0.51, 0.27)};
        const double density_value = 2.0 / basis.omega;
        std::vector<double> density(basis.nrxx, density_value);
        double* channels[] = {density.data()};
        Charge charge;
        charge.nspin = 1;
        charge.rho = channels;
        elecstate::PotSccs component(&basis, config, solver, 2.0, no_resume);
        ModuleBase::matrix potential(1, basis.nrxx);
        component.cal_v_eff(&charge, &cell, potential);
        ModuleBase::matrix force(2, 3);
        component.add_sol_force(cell, force);
        const double step = 1e-4;
        for (int ia = 0; ia < 2; ++ia)
        {
            for (int axis = 0; axis < 3; ++axis)
            {
                const double original = atom.tau[ia][axis];
                atom.tau[ia][axis] = original + step / length;
                potential.zero_out();
                component.cal_v_eff(&charge, &cell, potential);
                const double positive = component.get_energy();
                atom.tau[ia][axis] = original - step / length;
                potential.zero_out();
                component.cal_v_eff(&charge, &cell, potential);
                const double negative = component.get_energy();
                atom.tau[ia][axis] = original;
                const double finite_difference = -(positive - negative) / (2.0 * step);
                EXPECT_NEAR(force(ia, axis), finite_difference, 1e-8);
            }
        }
        potential.zero_out();
        component.cal_v_eff(&charge, &cell, potential);
        ModuleBase::matrix accumulated(2, 3);
        for (int ia = 0; ia < 2; ++ia)
        {
            for (int axis = 0; axis < 3; ++axis) { accumulated(ia, axis) = 7.0; }
        }
        component.add_sol_force(cell, accumulated);
        for (int ia = 0; ia < 2; ++ia)
        {
            for (int axis = 0; axis < 3; ++axis)
            {
                const double expected_force = 7.0 + force(ia, axis);
                EXPECT_NEAR(accumulated(ia, axis), expected_force, 1e-10);
            }
        }
    }
};

TEST_F(PotSccsTest, IndependentInstancesAndSpinChannels)
{
    UnitCell cell;
    Atom atom;
    cell.lat0 = length;
    cell.tpiba = tpiba;
    cell.omega = basis.omega;
    cell.ntype = 1;
    cell.nat = 1;
    cell.atoms = &atom;
    atom.na = 1;
    atom.ncpp.zv = 1.0;
    atom.tau = {ModuleBase::Vector3<double>(0.25, 0.25, 0.25)};
    const double density_value = 1.0 / basis.omega;
    std::vector<double> density(basis.nrxx, density_value);
    double* channels[] = {density.data(), density.data()};
    Charge charge;
    charge.nspin = 1;
    charge.rho = channels;
    Input_para input;
    ModuleSccs::SccsConfig config;
    ModuleSccs::PolarizationSolverParameters solver;
    elecstate::sccs_cfg_input(input, config, solver);
    config.cavity.epsilon_bulk = 5.0;
    elecstate::PotSccs first(&basis, config, solver, 1.0, no_resume);
    config.cavity.epsilon_bulk = 1.0;
    elecstate::PotSccs vacuum(&basis, config, solver, 1.0, no_resume);
    EXPECT_DOUBLE_EQ(first.get_energy(), 0.0);
    ModuleBase::matrix potential(1, basis.nrxx);
    first.cal_v_eff(&charge, &cell, potential);
    const double first_energy = first.get_energy();
    EXPECT_LT(first_energy, 0.0);
    double el = 0.0;
    double cav = 0.0;
    first.sol_energy(el, cav);
    EXPECT_DOUBLE_EQ(el, first_energy);
    EXPECT_DOUBLE_EQ(cav, 0.0);
    ModuleBase::matrix zero(1, basis.nrxx);
    vacuum.cal_v_eff(&charge, &cell, zero);
    EXPECT_NEAR(vacuum.get_energy(), 0.0, 1e-12);
    EXPECT_DOUBLE_EQ(first.get_energy(), first_energy);
    for (int ir = 0; ir < basis.nrxx; ++ir) { EXPECT_NEAR(zero(0, ir), 0.0, 1e-11); }
    for (double& value : density) { value *= 0.5; }
    charge.nspin = 2;
    ModuleBase::matrix spin_potential(2, basis.nrxx);
    first.cal_v_eff(&charge, &cell, spin_potential);
    EXPECT_NEAR(first.get_energy(), first_energy, 1e-12);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_NEAR(spin_potential(0, ir), potential(0, ir), 1e-11);
        EXPECT_DOUBLE_EQ(spin_potential(0, ir), spin_potential(1, ir));
    }
}

TEST_F(PotSccsTest, InputMapsOntoSccsConfig)
{
    // The custom defaults follow Environ environ_type input:
    // env_static_permittivity 1, rhomin 1e-4, rhomax 5e-3, no surface tension
    // or pressure.
    Input_para input;
    ModuleSccs::SccsConfig defaults;
    ModuleSccs::PolarizationSolverParameters default_solver;
    elecstate::sccs_cfg_input(input, defaults, default_solver);
    EXPECT_EQ(input.sccs_preset, "custom");
    EXPECT_DOUBLE_EQ(defaults.cavity.epsilon_bulk, 1.0);
    EXPECT_DOUBLE_EQ(defaults.cavity.density_min, 1e-4);
    EXPECT_DOUBLE_EQ(defaults.cavity.density_max, 5e-3);
    EXPECT_DOUBLE_EQ(defaults.surface_tension, 0.0);
    EXPECT_DOUBLE_EQ(defaults.pressure, 0.0);
    EXPECT_EQ(default_solver.max_iterations, 200);
    input.sccs_preset = "water-neutral";
    input.sccs_epsilon = 2.0; // replaced by the preset
    input.sccs_maxiter = 42;
    input.sccs_surface_eta = 1e-6;
    ModuleSccs::SccsConfig config;
    ModuleSccs::PolarizationSolverParameters solver;
    elecstate::sccs_cfg_input(input, config, solver);
    EXPECT_DOUBLE_EQ(config.cavity.epsilon_bulk, 78.3);
    EXPECT_DOUBLE_EQ(config.cavity.density_min, 1e-4);
    EXPECT_LT(config.pressure, 0.0);
    EXPECT_DOUBLE_EQ(config.surface_regularization, 1e-6);
    EXPECT_EQ(solver.max_iterations, 42);
    EXPECT_EQ(config.boundary, ModuleSccs::Boundary::Periodic);
    input.assume_isolated = "pcc_0d";
    elecstate::sccs_cfg_input(input, config, solver);
    EXPECT_EQ(config.boundary, ModuleSccs::Boundary::Pcc0d);
    input.assume_isolated = "pcc_2d";
    for (int axis = 0; axis < 3; ++axis)
    {
        input.pcc_2d_axis = axis;
        elecstate::sccs_cfg_input(input, config, solver);
        EXPECT_EQ(config.boundary, ModuleSccs::Boundary::Pcc2d);
        EXPECT_EQ(config.pcc_2d_axis, axis);
    }
}

TEST_F(PotSccsTest, InputMapsCavityModelKeys)
{
    Input_para input;
    ModuleSccs::SccsConfig config;
    ModuleSccs::PolarizationSolverParameters solver;
    elecstate::sccs_cfg_input(input, config, solver);
    EXPECT_DOUBLE_EQ(config.cavity.lowpass_p1, -1.0);
    EXPECT_DOUBLE_EQ(config.cavity.lowpass_p2, -1.0);
    EXPECT_FALSE(config.core_electrons);
    input.sccs_preset = "water-neutral";
    input.assume_isolated = "pcc_0d";
    input.sccs_lowpass_p1 = 10.0;
    input.sccs_lowpass_p2 = 5.0;
    input.sccs_solvent_mode = "full";
    input.sccs_corespread = {0.8, 0.0, -1.0};
    elecstate::sccs_cfg_input(input, config, solver);
    EXPECT_DOUBLE_EQ(config.cavity.lowpass_p1, 10.0);
    EXPECT_DOUBLE_EQ(config.cavity.lowpass_p2, 5.0);
    EXPECT_TRUE(config.core_electrons);
    EXPECT_EQ(config.core_spreads, input.sccs_corespread);
    // The cavity model keys leave the preset dielectric unchanged.
    EXPECT_DOUBLE_EQ(config.cavity.epsilon_bulk, 78.3);
    EXPECT_DOUBLE_EQ(config.start_drho, 0.0);
    EXPECT_FALSE(config.collect_report);
    input.sccs_start_drho = 1e-3;
    input.sccs_start_nmax = 3;
    input.sccs_debug = 2;
    elecstate::sccs_cfg_input(input, config, solver);
    EXPECT_DOUBLE_EQ(config.start_drho, 1e-3);
    EXPECT_EQ(config.start_nmax, 3);
    EXPECT_TRUE(config.collect_report);
    input.sccs_debug = 1;
    elecstate::sccs_cfg_input(input, config, solver);
    EXPECT_FALSE(config.collect_report);
}

TEST_F(PotSccsTest, IonicForceAddsRydbergDerivativeOnce)
{
    Input_para input;
    input.sccs_epsilon = 5.0;
    input.sccs_tol_rms = 1e-13;
    input.sccs_tol_max = 1e-12;
    ModuleSccs::SccsConfig config;
    ModuleSccs::PolarizationSolverParameters solver;
    elecstate::sccs_cfg_input(input, config, solver);
    check_solvation_force(config, solver);
}

// The cavity terms of core Gaussians, with surface tension and pressure.
TEST_F(PotSccsTest, FullCavityForceMatchesTheEnergyDerivative)
{
    Input_para input;
    input.sccs_epsilon = 1.0;
    input.sccs_tol_rms = 1e-13;
    input.sccs_tol_max = 1e-12;
    ModuleSccs::SccsConfig config;
    ModuleSccs::PolarizationSolverParameters solver;
    elecstate::sccs_cfg_input(input, config, solver);
    config.core_electrons = true;
    config.core_spreads = {2.0};
    config.surface_tension = 1e-5;
    config.pressure = 1e-6;
    check_solvation_force(config, solver);
}

TEST_F(PotSccsTest, DelayedActivationHasOneIrreversibleTransitionAndResumes)
{
    ModuleSccs::SccsConfig config;
    ModuleSccs::PolarizationSolverParameters solver;
    config.start_drho = 1e-3;
    config.start_nmax = 3;
    elecstate::PotSccs threshold(&basis, config, solver, 1.0, no_resume);
    EXPECT_FALSE(threshold.is_active());
    EXPECT_FALSE(threshold.try_activate(1, 0.1));
    EXPECT_TRUE(threshold.try_activate(2, 1e-3));
    EXPECT_TRUE(threshold.is_active());
    EXPECT_FALSE(threshold.try_activate(3, 0.2));
    EXPECT_FALSE(threshold.try_activate(1, 0.2));
    elecstate::PotSccs forced(&basis, config, solver, 1.0, no_resume);
    EXPECT_FALSE(forced.try_activate(2, 0.1));
    EXPECT_TRUE(forced.try_activate(3, 0.1));
    config.start_drho = 0.0;
    elecstate::PotSccs immediate(&basis, config, solver, 1.0, no_resume);
    EXPECT_TRUE(immediate.is_active());
    EXPECT_FALSE(immediate.try_activate(1, 0.0));
    // A later ionic step resumes an activation already reached.
    config.start_drho = 1e-3;
    elecstate::PotSccs first_ionic_step(&basis, config, solver, 1.0, no_resume);
    EXPECT_FALSE(first_ionic_step.is_active());
    elecstate::SccsResume resumed;
    resumed.active = true;
    elecstate::PotSccs next_ionic_step(&basis, config, solver, 1.0, resumed);
    EXPECT_TRUE(next_ionic_step.is_active());
    EXPECT_FALSE(next_ionic_step.try_activate(1, 0.1));
    EXPECT_DOUBLE_EQ(next_ionic_step.get_energy(), 0.0);
    const std::vector<double>* potential = next_ionic_step.sol_elec_pot();
    EXPECT_TRUE(potential->empty());
}

TEST_F(PotSccsTest, WaitingSccsContributesNoPotentialEnergyOrForce)
{
    UnitCell cell;
    Atom atom;
    cell.lat0 = length;
    cell.tpiba = tpiba;
    cell.omega = basis.omega;
    cell.ntype = 1;
    cell.nat = 1;
    cell.atoms = &atom;
    atom.na = 1;
    atom.ncpp.zv = 1.0;
    atom.tau = {ModuleBase::Vector3<double>(0.25, 0.25, 0.25)};
    const double value = 1.0 / basis.omega;
    std::vector<double> density(basis.nrxx, value);
    double* channels[] = {density.data()};
    Charge charge;
    charge.nspin = 1;
    charge.rho = channels;
    ModuleSccs::SccsConfig config = ModuleSccs::preset_config(ModuleSccs::Preset::WaterNeutral);
    config.start_drho = 1e-3;
    ModuleSccs::PolarizationSolverParameters solver;
    elecstate::PotSccs delayed(&basis, config, solver, 1.0, no_resume);
    ModuleBase::matrix potential(1, basis.nrxx);
    delayed.cal_v_eff(&charge, &cell, potential);
    EXPECT_DOUBLE_EQ(delayed.get_energy(), 0.0);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_DOUBLE_EQ(potential(0, ir), 0.0);
    }
    ModuleBase::matrix force(1, 3);
    force(0, 0) = 42.0;
    delayed.add_sol_force(cell, force);
    EXPECT_DOUBLE_EQ(force(0, 0), 42.0);
    EXPECT_TRUE(delayed.try_activate(2, 1e-4));
    delayed.cal_v_eff(&charge, &cell, potential);
    EXPECT_LT(delayed.get_energy(), 0.0);
}

TEST_F(PotSccsTest, FullCavityStructureChecksWarnOnceAndRejectWidthCount)
{
    UnitCell cell;
    Atom atom;
    cell.ntype = 1;
    cell.nat = 2;
    cell.atoms = &atom;
    atom.na = 2;
    atom.label = "X";
    atom.ncpp.psd = "Xx";
    ModuleSccs::SccsConfig config;
    config.core_electrons = true;
    config.core_spreads = {0.5};
    auto warnings = [&]() { return elecstate::check_sccs_cell(config, cell); };
    const std::string unknown = warnings();
    EXPECT_NE(unknown.find("unknown"), std::string::npos);
    atom.ncpp.psd = "O";
    EXPECT_EQ(warnings(), "");
    atom.ncpp.psd = "Xx";
    config.core_spreads = {0.5, 0.6};
    EXPECT_EQ(warnings(), "");
    config.core_electrons = false;
    config.core_spreads = {0.5};
    EXPECT_EQ(warnings(), "");
    config.core_electrons = true;
    config.core_spreads = {0.5, 0.6, 0.7};
    if (SccsTest::pool_size == 1)
    {
        testing::internal::CaptureStdout();
        EXPECT_EXIT(elecstate::check_sccs_cell(config, cell), ::testing::ExitedWithCode(1), "");
        testing::internal::GetCapturedStdout();
    }
}

TEST_F(PotSccsTest, ChargedDielectricCellWarnsAndRuns)
{
    UnitCell cell;
    Atom atom;
    cell.lat0 = length;
    cell.tpiba = tpiba;
    cell.omega = basis.omega;
    cell.ntype = 1;
    cell.nat = 1;
    cell.atoms = &atom;
    atom.na = 1;
    atom.ncpp.zv = 1.0;
    atom.tau = {ModuleBase::Vector3<double>(0.25, 0.25, 0.25)};
    Input_para input;
    ModuleSccs::SccsConfig config;
    ModuleSccs::PolarizationSolverParameters solver;
    elecstate::sccs_cfg_input(input, config, solver);
    config.cavity.epsilon_bulk = 5.0;
    EXPECT_TRUE(elecstate::check_charge(config, cell, 1.0).empty());
    EXPECT_NE(elecstate::check_charge(config, cell, 0.5).find("G = 0"), std::string::npos);
    ModuleSccs::SccsConfig vacuum_config = config;
    vacuum_config.cavity.epsilon_bulk = 1.0;
    EXPECT_TRUE(elecstate::check_charge(vacuum_config, cell, 0.5).empty());
    ModuleSccs::SccsConfig pcc_config = config;
    pcc_config.boundary = ModuleSccs::Boundary::Pcc0d;
    EXPECT_TRUE(elecstate::check_charge(pcc_config, cell, 0.5).empty());

    const double density_value = 0.5 / basis.omega;
    std::vector<double> density(basis.nrxx, density_value);
    double* channels[] = {density.data()};
    Charge charge;
    charge.nspin = 1;
    charge.rho = channels;
    // One log per rank: the MPI variant runs this test on several ranks at once.
    const std::string log_name = "pot_sccs_charge_warning_" + std::to_string(basis.poolrank) + ".log";
    std::ofstream& warning_log = GlobalV::ofs_warning;
    warning_log.open(log_name.c_str());
    elecstate::PotSccs matched(&basis, config, solver, 0.5, no_resume);
    ModuleBase::matrix matched_potential(1, basis.nrxx);
    matched.cal_v_eff(&charge, &cell, matched_potential);
    warning_log.flush();
    std::ifstream matched_log(log_name.c_str());
    std::stringstream matched_text;
    matched_text << matched_log.rdbuf();
    EXPECT_EQ(matched_text.str().find("SCCS grid electron count"), std::string::npos);
    EXPECT_TRUE(std::isfinite(matched.get_energy()));

    elecstate::PotSccs mismatched(&basis, config, solver, 1.0, no_resume);
    ModuleBase::matrix mismatched_potential(1, basis.nrxx);
    mismatched.cal_v_eff(&charge, &cell, mismatched_potential);
    warning_log.close();
    std::ifstream mismatched_log(log_name.c_str());
    std::stringstream mismatched_text;
    mismatched_text << mismatched_log.rdbuf();
    EXPECT_NE(mismatched_text.str().find("SCCS grid electron count differs"), std::string::npos);
    EXPECT_DOUBLE_EQ(mismatched.get_energy(), matched.get_energy());
    std::remove(log_name.c_str());
}
