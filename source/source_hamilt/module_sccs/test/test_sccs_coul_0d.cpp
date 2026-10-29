#include "sccs_test.h"
#include "../sccs_coul_0d.h"
#include "../sccs_cav_deriv.h"
#include "../sccs_functional.h"
#include "../sccs_parameters.h"
#include "../sccs_response.h"

#include "source_base/parallel_reduce.h"
#include "source_basis/module_pw/pw_grid_geometry.h"
#include "source_cell/cell_geometry.h"


class SccsPcc0dCoulombTest : public SccsTest::PwTest
{
protected:
    void SetUp() override
    {
        SccsTest::PwTest::SetUp();
        const ModuleBase::Matrix3 lattice;
        unitcell::OrthogonalCell cell;
        // Same cell tolerance as PotPcc.
        ASSERT_TRUE(unitcell::make_orthogonal_cell(lattice, length, 1e-10, cell));
        ASSERT_TRUE(elecstate::make_pcc_0d_parameters(cell, 1e-10, parameters));
        const double half_length = 0.5 * length;
        cell.origin = ModuleBase::Vector3<double>(half_length, half_length, half_length);
        ModulePW::grid_positions(basis, lattice, length, positions);
        for (ModuleBase::Vector3<double>& position : positions)
        {
            position = unitcell::relative_position(position, cell);
        }
    }

    double integrate_product(const std::vector<double>& left, const std::vector<double>& right) const
    {
        double value = 0.0;
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            value += left[ir] * right[ir];
        }
        Parallel_Reduce::reduce_pool(value);
        return value * basis.omega / basis.nxyz;
    }

    elecstate::Pcc0dParameters parameters;
    std::vector<ModuleBase::Vector3<double>> positions;
};

TEST_F(SccsPcc0dCoulombTest, CorrectionMatchesExistingPccEnergyAndIsSymmetric)
{
    std::vector<double> charge = cosine_mode(0);
    std::vector<double> direction = cosine_mode(1);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        charge[ir] = SccsTest::charge_mean + SccsTest::charge_amplitude * charge[ir];
        direction[ir] = SccsTest::probe_mean + SccsTest::probe_amplitude * direction[ir];
    }
    ModuleSccs::PeriodicCoulombOperator periodic(basis, tpiba);
    ModuleSccs::Pcc0dCoulombOperator corrected(basis, tpiba, positions, parameters);
    EXPECT_TRUE(corrected.has_correction());
    std::vector<double> periodic_potential;
    std::vector<double> potential;
    std::vector<double> direction_potential;
    periodic.apply_potential(charge, periodic_potential);
    corrected.apply_potential(charge, potential);
    corrected.apply_potential(direction, direction_potential);
    const double charge_direction = integrate_product(charge, direction_potential);
    const double direction_charge = integrate_product(direction, potential);
    EXPECT_NEAR(charge_direction, direction_charge, 1e-13);

    const double dv = basis.omega / basis.nxyz;
    const double* charge_data = charge.data();
    const ModuleBase::Vector3<double>* position_data = positions.data();
    elecstate::ChargeMoments moments = elecstate::charge_moments(charge_data, position_data, basis.nrxx, dv);
    Parallel_Reduce::reduce_pool(moments.charge);
    Parallel_Reduce::reduce_pool(moments.dipole.x);
    Parallel_Reduce::reduce_pool(moments.dipole.y);
    Parallel_Reduce::reduce_pool(moments.dipole.z);
    Parallel_Reduce::reduce_pool(moments.second_moment);
    std::vector<double> correction(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        correction[ir] = potential[ir] - periodic_potential[ir];
        const double expected = elecstate::pcc_0d_potential(moments, positions[ir], parameters);
        EXPECT_NEAR(correction[ir], expected, 1e-14);
    }
    const double grid_energy = 0.5 * integrate_product(charge, correction);
    const double moment_energy = elecstate::pcc_0d_energy(moments, parameters);
    EXPECT_NEAR(grid_energy, moment_energy, 1e-13);
}

TEST_F(SccsPcc0dCoulombTest, ChargedUniformDielectricPreservesGaugeAndVacuumSubtraction)
{
    std::vector<double> charge = cosine_mode(0);
    std::vector<double> direction = cosine_mode(1);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        charge[ir] = SccsTest::charge_mean + SccsTest::charge_amplitude * charge[ir];
        // Perturbation direction for the finite-difference check, not proportional to charge.
        direction[ir] = 1e-3 + 3e-4 * direction[ir];
    }
    ModuleSccs::Pcc0dCoulombOperator coulomb(basis, tpiba, positions, parameters);
    const std::vector<double> density(basis.nrxx, 0.0);
    const std::vector<double> cold;
    std::vector<double> vacuum;
    coulomb.apply_potential(charge, vacuum);
    const std::vector<double> ones(basis.nrxx, 1.0);
    const double vacuum_mean = integrate_product(vacuum, ones) / basis.omega;
    EXPECT_GT(std::abs(vacuum_mean), 1e-3);
    ModuleSccs::PolarizationSolverParameters solver;
    const double dielectrics[] = {1.0, 5.0};
    for (double epsilon : dielectrics)
    {
        ModuleSccs::SccsConfig config = ModuleSccs::preset_config(ModuleSccs::Preset::Vacuum);
        config.cavity.epsilon_bulk = epsilon;
        ModuleSccs::SccsResponse response;
        ModuleSccs::solve_response(density, charge, config.cavity, solver, cold,
                                   basis, tpiba, coulomb, response);
        EXPECT_EQ(response.polarization.iterations, 1);
        const double solute_sum = integrate_product(charge, ones);
        const double screening = (1.0 / epsilon - 1.0) * solute_sum;
        EXPECT_NEAR(response.far_field_polarization_charge, screening, 1e-13);
        ModuleSccs::FunctionalResult result;
        ModuleSccs::eval_functional(charge, response, config, basis, tpiba,
                                    coulomb, result);
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            const double dielectric_potential = vacuum[ir] / epsilon;
            const double reaction = dielectric_potential - vacuum[ir];
            EXPECT_NEAR(response.polarization.potential[ir], dielectric_potential, 1e-13);
            EXPECT_NEAR(result.reaction_potential[ir], reaction, 1e-13);
            const double expected_electron_potential = -reaction;
        EXPECT_NEAR(result.electron_potential[ir], expected_electron_potential, 1e-13);
        }
        const double expected_energy = 0.5 * (1.0 / epsilon - 1.0) * integrate_product(charge, vacuum);
        EXPECT_NEAR(result.reaction_energy, expected_energy, 1e-13);
        // Electron potential carries the opposite sign to signed solute charge.
        const double predicted = -integrate_product(result.electron_potential, direction);
        // Relative perturbation of 1e-4 x 1e-3: the central difference error is below 1e-10.
        const double step = 1e-4;
        const double signs[] = {-1.0, 1.0};
        double energy[2];
        int index = 0;
        for (double sign : signs)
        {
            std::vector<double> shifted = charge;
            for (int ir = 0; ir < basis.nrxx; ++ir)
            {
                shifted[ir] += sign * step * direction[ir];
            }
            ModuleSccs::SccsResponse shifted_response;
            ModuleSccs::FunctionalResult shifted_result;
            ModuleSccs::solve_response(density, shifted, config.cavity, solver, cold,
                                       basis, tpiba, coulomb, shifted_response);
            ModuleSccs::eval_functional(shifted, shifted_response, config, basis, tpiba,
                                        coulomb, shifted_result);
            energy[index++] = shifted_result.reaction_energy;
        }
        const double derivative = (energy[1] - energy[0]) / (2.0 * step);
        EXPECT_NEAR(derivative, predicted, 1e-10);
    }
}

TEST_F(SccsPcc0dCoulombTest, NonuniformResponseSatisfiesFixedPointAndFarFieldCharge)
{
    ModuleSccs::CavityParameters cavity;
    // Water-neutral switching range.
    cavity.density_min = 1e-4;
    cavity.density_max = 5e-3;
    cavity.epsilon_bulk = 5.0;
    const std::vector<double> mode = cosine_mode(0);
    std::vector<double> density(basis.nrxx);
    std::vector<double> charge(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        density[ir] = SccsTest::cavity_mean + SccsTest::cavity_amplitude * mode[ir];
        charge[ir] = SccsTest::charge_mean + SccsTest::charge_amplitude * mode[ir];
    }
    ModuleSccs::Pcc0dCoulombOperator coulomb(basis, tpiba, positions, parameters);
    ModuleSccs::PolarizationSolverParameters solver;
    solver.tolerance_rms = 1e-12;
    solver.tolerance_max = 1e-12;
    const std::vector<double> cold;
    ModuleSccs::SccsResponse response;
    ModuleSccs::solve_response(density, charge, cavity, solver, cold,
                               basis, tpiba, coulomb, response);
    EXPECT_GT(response.polarization.iterations, 1);
    const bool open_boundary = true;
    ModuleSccs::SccsResponse derivatives;
    ModuleSccs::CavityDerivatives cavity_derivatives;
    const std::vector<double>& coefficient = cavity_derivatives.coefficient;
    ModuleSccs::make_cav_derivs(density, cavity, basis, tpiba, open_boundary,
                                derivatives, cavity_derivatives);
    std::vector<double> source(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double invsqrt = 1.0 / std::sqrt(response.epsilon[ir]);
        source[ir] = (charge[ir] - coefficient[ir] * response.polarization.potential[ir]) * invsqrt;
    }
    std::vector<double> reconstructed;
    coulomb.apply_potential(source, reconstructed);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double expected = reconstructed[ir] / std::sqrt(response.epsilon[ir]);
        EXPECT_NEAR(response.polarization.potential[ir], expected, 1e-10);
    }
    const std::vector<double> ones(basis.nrxx, 1.0);
    const double source_sum = integrate_product(source, ones);
    const double solute_sum = integrate_product(charge, ones);
    const double expected_screening = source_sum / std::sqrt(cavity.epsilon_bulk) - solute_sum;
    EXPECT_NEAR(response.far_field_polarization_charge, expected_screening, 1e-13);
    ModuleSccs::SccsResponse restarted;
    ModuleSccs::solve_response(density, charge, cavity, solver, response.restart_potential,
                               basis, tpiba, coulomb, restarted);
    EXPECT_TRUE(restarted.polarization.warm_started);
    EXPECT_LT(restarted.polarization.iterations, response.polarization.iterations);
    EXPECT_NEAR(restarted.far_field_polarization_charge, response.far_field_polarization_charge, 1e-10);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_NEAR(restarted.polarization.potential[ir], response.polarization.potential[ir], 1e-10);
    }
}
