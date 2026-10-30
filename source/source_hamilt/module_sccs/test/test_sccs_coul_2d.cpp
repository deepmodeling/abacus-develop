#include "sccs_test.h"
#include "../sccs_coul_2d.h"
#include "../sccs_cav_deriv.h"
#include "../sccs_functional.h"
#include "../sccs_parameters.h"
#include "../sccs_response.h"

#include "source_base/parallel_reduce.h"
#include "source_basis/module_pw/pw_grid_geometry.h"
#include "source_cell/cell_geometry.h"


class SccsPcc2dCoulombTest : public SccsTest::PwTest
{
protected:
    void SetUp() override
    {
        // Oblique periodic plane, with a perpendicular open vector along x.
        const ModuleBase::Matrix3 lattice(0.0, 1.0, 0.0,
                                          0.0, 0.2, 1.0,
                                          1.5, 0.0, 0.0);
        set_up_basis(lattice);
        unitcell::SlabCell cell;
        ASSERT_TRUE(unitcell::make_slab_cell(lattice, length, 2, 1e-10, cell));
        cell.origin = 0.5 * cell.length;
        normal = cell.normal;
        parameters = elecstate::make_pcc_2d_parameters(cell);
        ModulePW::grid_positions(basis, lattice, length, positions);
        for (ModuleBase::Vector3<double>& position : positions)
        {
            const double coordinate = unitcell::relative_coordinate(position, cell);
            position = normal * coordinate;
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

    elecstate::Pcc2dParameters parameters;
    ModuleBase::Vector3<double> normal;
    std::vector<ModuleBase::Vector3<double>> positions;
};

TEST_F(SccsPcc2dCoulombTest, ChargedPlanarModeMatchesAnalyticGauge)
{
    // One cosine mode along the open axis on a uniform part: the mode has the
    // analytic potential 4 pi amplitude / k^2 and the uniform part carries the
    // monopole, so the charged gauge -pi L Q / (3 A) is tested as well.
    const double density_mean = SccsTest::charge_mean;
    const double amplitude = SccsTest::charge_amplitude;
    double monopole = 0.0;
    double dipole = 0.0;
    double second = 0.0;
    const double plane_volume = basis.omega / basis.nz;
    for (int iz = 0; iz < basis.nz; ++iz)
    {
        const double angle = ModuleBase::TWO_PI * iz / basis.nz;
        const double weight = (density_mean + amplitude * std::cos(angle)) * plane_volume;
        const double coordinate = parameters.length * iz / basis.nz - 0.5 * parameters.length;
        monopole += weight;
        dipole += weight * coordinate;
        second += weight * coordinate * coordinate;
    }
    std::vector<double> charge(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const int iz = ir % basis.nplane + basis.startz_current;
        const double angle = ModuleBase::TWO_PI * iz / basis.nz;
        charge[ir] = density_mean + amplitude * std::cos(angle);
    }
    ModuleSccs::Pcc2dCoulombOperator coulomb(basis, tpiba, positions, normal, parameters);
    EXPECT_TRUE(coulomb.has_correction());
    std::vector<double> potential;
    coulomb.apply_potential(charge, potential);
    const double wavevector = ModuleBase::TWO_PI / parameters.length;
    const double mode_factor = ModuleBase::FOUR_PI * amplitude / (wavevector * wavevector);
    const double constant = -ModuleBase::PI * parameters.length / (3.0 * parameters.area);
    const double factor = -2.0 * ModuleBase::PI / basis.omega;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const int iz = ir % basis.nplane + basis.startz_current;
        const double angle = ModuleBase::TWO_PI * iz / basis.nz;
        const double coordinate = positions[ir] * normal;
        const double multipole = monopole * coordinate * coordinate - 2.0 * dipole * coordinate + second;
        const double expected = mode_factor * std::cos(angle) + constant * monopole + factor * multipole;
        // The operator is exact for this charge; the tolerance is round-off.
        EXPECT_NEAR(potential[ir], expected, 1e-12);
    }
    std::vector<double> direction = cosine_mode(0);
    for (double& value : direction)
    {
        value = SccsTest::probe_mean + SccsTest::probe_amplitude * value;
    }
    std::vector<double> direction_potential;
    coulomb.apply_potential(direction, direction_potential);
    const double left = integrate_product(charge, direction_potential);
    const double right = integrate_product(direction, potential);
    EXPECT_NEAR(left, right, 1e-13);
}

TEST_F(SccsPcc2dCoulombTest, UniformDielectricRetainsChargedGaugeAndReactionEnergy)
{
    std::vector<double> charge = cosine_mode(0);
    for (double& value : charge)
    {
        value = SccsTest::charge_mean + SccsTest::charge_amplitude * value;
    }
    ModuleSccs::Pcc2dCoulombOperator coulomb(basis, tpiba, positions, normal, parameters);
    std::vector<double> vacuum;
    coulomb.apply_potential(charge, vacuum);
    const std::vector<double> density(basis.nrxx, 0.0);
    const std::vector<double> ones(basis.nrxx, 1.0);
    const std::vector<double> cold;
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
        ModuleSccs::FunctionalResult result;
        ModuleSccs::eval_functional(charge, response, config, basis, tpiba,
                                    coulomb, result);
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            const double expected_potential = vacuum[ir] / epsilon;
            const double expected_reaction = expected_potential - vacuum[ir];
            EXPECT_NEAR(response.polarization.potential[ir], expected_potential, 1e-13);
            EXPECT_NEAR(result.reaction_potential[ir], expected_reaction, 1e-13);
            const double expected_electron_potential = -expected_reaction;
        EXPECT_NEAR(result.electron_potential[ir], expected_electron_potential, 1e-13);
        }
        const double expected_energy = 0.5 * (1.0 / epsilon - 1.0) * integrate_product(charge, vacuum);
        EXPECT_NEAR(result.reaction_energy, expected_energy, 1e-13);
        const double expected_screening = (1.0 / epsilon - 1.0) * integrate_product(charge, ones);
        EXPECT_NEAR(response.far_field_polarization_charge, expected_screening, 1e-13);
    }
}

TEST_F(SccsPcc2dCoulombTest, NonuniformResponseSatisfiesFixedPointAndWarmStart)
{
    const std::vector<double> mode = cosine_mode(0);
    std::vector<double> density(basis.nrxx);
    std::vector<double> charge(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        density[ir] = SccsTest::cavity_mean + SccsTest::cavity_amplitude * mode[ir];
        charge[ir] = SccsTest::charge_mean + SccsTest::charge_amplitude * mode[ir];
    }
    ModuleSccs::CavityParameters cavity;
    // Water-neutral switching range.
    cavity.density_min = 1e-4;
    cavity.density_max = 5e-3;
    cavity.epsilon_bulk = 5.0;
    ModuleSccs::PolarizationSolverParameters solver;
    solver.tolerance_rms = 1e-12;
    solver.tolerance_max = 1e-12;
    ModuleSccs::Pcc2dCoulombOperator coulomb(basis, tpiba, positions, normal, parameters);
    const std::vector<double> cold;
    ModuleSccs::SccsResponse response;
    ModuleSccs::solve_response(density, charge, cavity, solver, cold,
                               basis, tpiba, coulomb, response);
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
    ModuleSccs::SccsResponse restarted;
    ModuleSccs::solve_response(density, charge, cavity, solver, response.restart_potential,
                               basis, tpiba, coulomb, restarted);
    EXPECT_TRUE(restarted.polarization.warm_started);
    EXPECT_LT(restarted.polarization.iterations, response.polarization.iterations);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_NEAR(restarted.polarization.potential[ir], response.polarization.potential[ir], 1e-10);
    }
}
