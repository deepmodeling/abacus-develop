#include "source_estate/module_pot/pot_sccs.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_hamilt/module_sccs/test/sccs_test.h"

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

using PotSccsTest = SccsTest::PwTest;

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
    elecstate::make_sccs_config_from_input(input, config, solver);
    config.cavity.epsilon_bulk = 5.0;
    elecstate::PotSccs first(&basis, config, solver);
    config.cavity.epsilon_bulk = 1.0;
    elecstate::PotSccs vacuum(&basis, config, solver);
    EXPECT_DOUBLE_EQ(first.get_energy(), 0.0);
    ModuleBase::matrix potential(1, basis.nrxx);
    first.cal_v_eff(&charge, &cell, potential);
    const double first_energy = first.get_energy();
    EXPECT_LT(first_energy, 0.0);
    double el = 0.0;
    double cav = 0.0;
    first.get_solvation_energy(el, cav);
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
    Input_para input;
    input.sccs_preset = "water-neutral";
    input.sccs_epsilon = 2.0; // replaced by the preset
    input.sccs_maxiter = 42;
    input.sccs_surface_eta = 1e-6;
    ModuleSccs::SccsConfig config;
    ModuleSccs::PolarizationSolverParameters solver;
    elecstate::make_sccs_config_from_input(input, config, solver);
    EXPECT_DOUBLE_EQ(config.cavity.epsilon_bulk, 78.3);
    EXPECT_DOUBLE_EQ(config.cavity.density_min, 1e-4);
    EXPECT_LT(config.pressure, 0.0);
    EXPECT_DOUBLE_EQ(config.surface_regularization, 1e-6);
    EXPECT_EQ(solver.max_iterations, 42);
}

TEST_F(PotSccsTest, IonicForceAddsRydbergDerivativeOnce)
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
    Input_para input;
    input.sccs_epsilon = 5.0;
    input.sccs_tol_rms = 1e-13;
    input.sccs_tol_max = 1e-12;
    ModuleSccs::SccsConfig config;
    ModuleSccs::PolarizationSolverParameters solver;
    elecstate::make_sccs_config_from_input(input, config, solver);
    elecstate::PotSccs component(&basis, config, solver);
    ModuleBase::matrix potential(1, basis.nrxx);
    component.cal_v_eff(&charge, &cell, potential);
    ModuleBase::matrix force(2, 3);
    component.add_solvation_force(cell, force);
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
    component.add_solvation_force(cell, accumulated);
    for (int ia = 0; ia < 2; ++ia)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            const double expected_force = 7.0 + force(ia, axis);
            EXPECT_NEAR(accumulated(ia, axis), expected_force, 1e-10);
        }
    }
}
