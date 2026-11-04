#include "../snap_calculator.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>

namespace
{

ModuleMD::SNAP::ModelData make_model()
{
    ModuleMD::SNAP::ModelData model;
    model.parameters.rcut = 4.5;
    model.parameters.rcutfac = 1.0;
    model.parameters.rmin = 0.0;
    model.parameters.rfac0 = 0.99363;
    model.parameters.twojmax = 0;
    model.parameters.elements.push_back("Si");
    model.parameters.element_weights.push_back(1.0);
    model.parameters.element_radii.push_back(2.25);
    model.coefficients.push_back(std::vector<double>{2.0, 3.0});
    return model;
}

ModuleMD::SNAP::ModelData make_multi_element_model()
{
    ModuleMD::SNAP::ModelData model;
    model.parameters.rcut = 4.5;
    model.parameters.rcutfac = 1.0;
    model.parameters.rmin = 0.0;
    model.parameters.rfac0 = 0.99363;
    model.parameters.twojmax = 0;
    model.parameters.bzero_flag = false;
    model.parameters.elements.push_back("Si");
    model.parameters.elements.push_back("C");
    model.parameters.element_weights.push_back(1.0);
    model.parameters.element_weights.push_back(0.8);
    model.parameters.element_radii.push_back(2.25);
    model.parameters.element_radii.push_back(1.7);
    model.coefficients.push_back(std::vector<double>{2.0, 3.0});
    model.coefficients.push_back(std::vector<double>{5.0, 7.0});
    return model;
}

ModuleMD::SNAP::ModelData make_chemical_model()
{
    ModuleMD::SNAP::ModelData model;
    model.parameters.rcut = 4.5;
    model.parameters.rcutfac = 1.0;
    model.parameters.rmin = 0.0;
    model.parameters.rfac0 = 0.99363;
    model.parameters.twojmax = 0;
    model.parameters.bnorm_flag = true;
    model.parameters.chem_flag = true;
    model.parameters.wselfall_flag = true;
    model.parameters.elements.push_back("Si");
    model.parameters.elements.push_back("C");
    model.parameters.element_weights.push_back(1.0);
    model.parameters.element_weights.push_back(0.8);
    model.parameters.element_radii.push_back(2.25);
    model.parameters.element_radii.push_back(1.7);
    model.coefficients.push_back(std::vector<double>(9, 0.0));
    model.coefficients.push_back(std::vector<double>(9, 0.0));
    model.coefficients[0][1] = 2.0;
    model.coefficients[1][1] = 3.0;
    return model;
}

ModuleMD::SNAP::ModelData make_angular_quadratic_model()
{
    ModuleMD::SNAP::ModelData model = make_model();
    model.parameters.twojmax = 2;
    model.parameters.quadratic_flag = true;
    model.coefficients[0] = {0.2, 0.3, -0.1, 0.5, 0.2, -0.4,
                             0.25, 0.0, 0.0, 0.0, -0.15,
                             0.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0, 0.0};
    return model;
}

struct SmallCellResult
{
    double energy = 0.0;
    double virial_xx = 0.0;
    std::array<std::array<double, 3>, 3> forces = {};
};

SmallCellResult evaluate_small_cell(const ModuleMD::SNAP::Calculator& calculator,
                                    const std::array<std::array<double, 3>, 3>& positions)
{
    SmallCellResult total;
    ModuleMD::SNAP::CellAccumulator accumulator(positions.size());
    for (std::size_t center = 0; center < positions.size(); ++center)
    {
        std::vector<ModuleMD::SNAP::Neighbor> neighbors;
        std::vector<std::size_t> neighbor_indices;
        for (std::size_t atom = 0; atom < positions.size(); ++atom)
        {
            if (atom == center)
            {
                continue;
            }
            ModuleMD::SNAP::Neighbor neighbor;
            neighbor.x = positions[atom][0] - positions[center][0];
            neighbor.y = positions[atom][1] - positions[center][1];
            neighbor.z = positions[atom][2] - positions[center][2];
            neighbor.cutoff_radius = 4.5;
            neighbors.push_back(neighbor);
            neighbor_indices.push_back(atom);
        }
        const ModuleMD::SNAP::ForceResult result = calculator.compute_with_analytic_force(neighbors, 0);
        accumulator.add(center, neighbor_indices, neighbors, result);
    }
    total.energy = accumulator.energy();
    total.virial_xx = accumulator.virial()[0];
    for (std::size_t atom = 0; atom < positions.size(); ++atom)
    {
        total.forces[atom] = accumulator.forces()[atom];
    }
    return total;
}

TEST(SnapCalculator, ComputesLinearAtomicEnergy)
{
    const ModuleMD::SNAP::Calculator calculator(make_model());
    ModuleMD::SNAP::Neighbor neighbor;
    neighbor.x = 1.0;
    const ModuleMD::SNAP::EnergyResult result =
        calculator.compute(std::vector<ModuleMD::SNAP::Neighbor>(1, neighbor), 0);

    ASSERT_EQ(result.descriptor.size(), static_cast<std::size_t>(1));
    EXPECT_TRUE(std::isfinite(result.descriptor[0]));
    EXPECT_DOUBLE_EQ(result.energy, 2.0 + 3.0 * result.descriptor[0]);
}

TEST(SnapCalculator, SupportsMultipleCentralElements)
{
    const ModuleMD::SNAP::Calculator calculator(make_multi_element_model());
    const std::vector<ModuleMD::SNAP::Neighbor> neighbors;
    const ModuleMD::SNAP::EnergyResult silicon = calculator.compute(neighbors, 0);
    const ModuleMD::SNAP::EnergyResult carbon = calculator.compute(neighbors, 1);

    EXPECT_DOUBLE_EQ(silicon.energy, 2.0 + 3.0 * silicon.descriptor[0]);
    EXPECT_DOUBLE_EQ(carbon.energy, 5.0 + 7.0 * carbon.descriptor[0]);
    EXPECT_NE(silicon.energy, carbon.energy);
}

TEST(SnapCalculator, SupportsChemicalChannels)
{
    const ModuleMD::SNAP::Calculator calculator(make_chemical_model());
    ModuleMD::SNAP::Neighbor silicon;
    silicon.x = 1.0;
    silicon.cutoff_radius = 4.5;
    silicon.element_index = 0;
    ModuleMD::SNAP::Neighbor carbon = silicon;
    carbon.y = 1.0;
    carbon.element_index = 1;
    const std::vector<ModuleMD::SNAP::Neighbor> neighbors = {silicon, carbon};

    const ModuleMD::SNAP::ForceResult analytic = calculator.compute_with_analytic_force(neighbors, 0);
    const ModuleMD::SNAP::ForceResult finite_difference
        = calculator.compute_with_finite_difference_force(neighbors, 0, 1.0e-6);

    EXPECT_EQ(analytic.energy.descriptor.size(), static_cast<std::size_t>(8));
    EXPECT_NEAR(analytic.energy.energy, finite_difference.energy.energy, 1.0e-12);
    for (int direction = 0; direction < 3; ++direction)
    {
        EXPECT_NEAR(analytic.center_force[static_cast<std::size_t>(direction)],
                    finite_difference.center_force[static_cast<std::size_t>(direction)],
                    1.0e-5);
        EXPECT_NEAR(analytic.neighbor_forces[0][static_cast<std::size_t>(direction)],
                    finite_difference.neighbor_forces[0][static_cast<std::size_t>(direction)],
                    1.0e-5);
    }
}

TEST(SnapCalculator, ChemicalAngularChannelsFollowElementOrder)
{
    ModuleMD::SNAP::ModelData model = make_chemical_model();
    model.parameters.twojmax = 2;
    const std::vector<double> zero_coefficients(41, 0.0);
    model.coefficients.assign(2, zero_coefficients);
    // Bispectrum index 2 couples the Si and C angular channels.
    const std::size_t original_index = 1 + 2 * 5 + 2;
    model.coefficients[0][original_index] = 1.7;
    const ModuleMD::SNAP::Calculator calculator(model);

    ModuleMD::SNAP::Neighbor silicon;
    silicon.x = 1.0;
    silicon.y = 0.1;
    silicon.z = 0.2;
    silicon.cutoff_radius = 4.5;
    silicon.element_index = 0;
    ModuleMD::SNAP::Neighbor carbon = silicon;
    carbon.x = 0.3;
    carbon.y = 1.1;
    carbon.z = -0.4;
    carbon.weight = 0.8;
    carbon.element_index = 1;
    const std::vector<ModuleMD::SNAP::Neighbor> neighbors = {silicon, carbon};
    const ModuleMD::SNAP::ForceResult original = calculator.compute_with_analytic_force(neighbors, 0);
    const ModuleMD::SNAP::ForceResult finite = calculator.compute_with_finite_difference_force(neighbors, 0, 1.0e-6);
    ModuleMD::SNAP::Neighbor rotated_carbon = carbon;
    rotated_carbon.x = carbon.y;
    rotated_carbon.y = carbon.x;
    const std::vector<ModuleMD::SNAP::Neighbor> rotated_neighbors = {silicon, rotated_carbon};
    const double rotated_energy = calculator.compute(rotated_neighbors, 0).energy;
    EXPECT_GT(std::abs(original.energy.energy - rotated_energy), 1.0e-6);

    ModuleMD::SNAP::ModelData reordered = model;
    reordered.parameters.elements = {"C", "Si"};
    reordered.parameters.element_weights = {0.8, 1.0};
    reordered.parameters.element_radii = {1.7, 2.25};
    reordered.coefficients.assign(2, zero_coefficients);
    // The same physical channel triple is (Si, C, Si) = (1, 0, 1) after reordering.
    const std::size_t reordered_index = 1 + 5 * 5 + 2;
    reordered.coefficients[1][reordered_index] = 1.7;
    const ModuleMD::SNAP::Calculator reordered_calculator(reordered);
    silicon.element_index = 1;
    carbon.element_index = 0;
    const std::vector<ModuleMD::SNAP::Neighbor> reordered_neighbors = {silicon, carbon};
    const ModuleMD::SNAP::ForceResult permuted
        = reordered_calculator.compute_with_analytic_force(reordered_neighbors, 1);

    EXPECT_NEAR(original.energy.energy, permuted.energy.energy, 1.0e-10);
    for (std::size_t atom = 0; atom < neighbors.size(); ++atom)
    {
        double force_size = 0.0;
        for (int direction = 0; direction < 3; ++direction)
        {
            const std::size_t component = static_cast<std::size_t>(direction);
            const double force = original.neighbor_forces[atom][component];
            force_size += std::abs(force);
            EXPECT_NEAR(force, permuted.neighbor_forces[atom][component], 1.0e-10);
            EXPECT_NEAR(force, finite.neighbor_forces[atom][component], 1.0e-5);
        }
        EXPECT_GT(force_size, 1.0e-6);
    }
}

TEST(SnapCalculator, RejectsMismatchedCoefficientCount)
{
    ModuleMD::SNAP::ModelData model = make_model();
    model.coefficients[0].push_back(4.0);
    try
    {
        const ModuleMD::SNAP::Calculator calculator(model);
        FAIL() << "Expected a coefficient-count error";
    }
    catch (const std::invalid_argument& error)
    {
        EXPECT_NE(std::string(error.what()).find("expected 2, found 3"), std::string::npos);
    }
}

TEST(SnapCalculator, CachedSelfDescriptorMatchesEnergyAndForce)
{
    ModuleMD::SNAP::ModelData model = make_model();
    const ModuleMD::SNAP::Calculator calculator(model);
    ModuleMD::SNAP::Neighbor neighbor;
    neighbor.x = 1.0;
    neighbor.y = 0.4;
    const std::vector<ModuleMD::SNAP::Neighbor> neighbors(1, neighbor);
    const ModuleMD::SNAP::EnergyResult energy = calculator.compute(neighbors, 0);
    const ModuleMD::SNAP::ForceResult analytic = calculator.compute_with_analytic_force(neighbors, 0);
    const ModuleMD::SNAP::ForceResult finite = calculator.compute_with_finite_difference_force(neighbors, 0, 1.0e-6);
    EXPECT_DOUBLE_EQ(energy.energy, analytic.energy.energy);
    ASSERT_EQ(energy.descriptor.size(), analytic.energy.descriptor.size());
    EXPECT_DOUBLE_EQ(energy.descriptor[0], analytic.energy.descriptor[0]);
    for (int direction = 0; direction < 3; ++direction)
    {
        const std::size_t index = static_cast<std::size_t>(direction);
        EXPECT_NEAR(analytic.neighbor_forces[0][index], finite.neighbor_forces[0][index], 1.0e-6);
    }
}

TEST(SnapCalculator, CachedSelfDescriptorSelectsCentralElement)
{
    ModuleMD::SNAP::ModelData model = make_chemical_model();
    model.parameters.wselfall_flag = false;
    model.coefficients[0] = {2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0};
    model.coefficients[1] = {7.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.0};
    const ModuleMD::SNAP::Calculator calculator(model);
    const std::vector<ModuleMD::SNAP::Neighbor> neighbors;
    for (std::size_t element = 0; element < model.parameters.elements.size(); ++element)
    {
        const ModuleMD::SNAP::EnergyResult energy = calculator.compute(neighbors, element);
        const ModuleMD::SNAP::ForceResult force = calculator.compute_with_analytic_force(neighbors, element);
        EXPECT_DOUBLE_EQ(energy.energy, model.coefficients[element][0]);
        EXPECT_DOUBLE_EQ(force.energy.energy, model.coefficients[element][0]);
        for (std::size_t index = 0; index < energy.descriptor.size(); ++index)
        {
            EXPECT_DOUBLE_EQ(energy.descriptor[index], 0.0);
        }
    }
}

TEST(SnapCalculator, RejectsInvalidCentralElement)
{
    const ModuleMD::SNAP::Calculator calculator(make_model());
    EXPECT_THROW(calculator.compute(std::vector<ModuleMD::SNAP::Neighbor>(), 1), std::out_of_range);
}

TEST(SnapCalculator, FiniteDifferenceForcesConservePairContribution)
{
    const ModuleMD::SNAP::Calculator calculator(make_model());
    ModuleMD::SNAP::Neighbor neighbor;
    neighbor.x = 1.0;
    const ModuleMD::SNAP::ForceResult result = calculator.compute_with_finite_difference_force(
        std::vector<ModuleMD::SNAP::Neighbor>(1, neighbor), 0, 1.0e-5);

    for (int direction = 0; direction < 3; ++direction)
    {
        EXPECT_DOUBLE_EQ(result.center_force[static_cast<std::size_t>(direction)]
                             + result.neighbor_forces[0][static_cast<std::size_t>(direction)],
                         0.0);
    }
}

TEST(SnapCalculator, AnalyticForceMatchesFiniteDifference)
{
    const ModuleMD::SNAP::Calculator calculator(make_model());
    ModuleMD::SNAP::Neighbor neighbor;
    neighbor.x = 1.0;
    neighbor.y = 0.4;
    neighbor.z = -0.2;
    const std::vector<ModuleMD::SNAP::Neighbor> neighbors(1, neighbor);
    const ModuleMD::SNAP::ForceResult analytic = calculator.compute_with_analytic_force(neighbors, 0);
    const ModuleMD::SNAP::ForceResult finite = calculator.compute_with_finite_difference_force(neighbors, 0, 1.0e-6);
    for (int direction = 0; direction < 3; ++direction)
    {
        EXPECT_NEAR(analytic.neighbor_forces[0][static_cast<std::size_t>(direction)],
                    finite.neighbor_forces[0][static_cast<std::size_t>(direction)],
                    1.0e-6);
    }
}

TEST(SnapCalculator, InnerSwitchMatchesFiniteDifference)
{
    ModuleMD::SNAP::ModelData model = make_model();
    model.parameters.twojmax = 2;
    model.parameters.switchinner_flag = true;
    model.parameters.sinner.push_back(1.5);
    model.parameters.dinner.push_back(0.4);
    model.coefficients[0] = {0.2, 0.3, -0.1, 0.5, 0.2, -0.4};
    const ModuleMD::SNAP::Calculator calculator(model);
    ModuleMD::SNAP::Neighbor first;
    first.x = 1.4;
    first.y = 0.3;
    first.z = 0.2;
    first.sinner = 1.5;
    first.dinner = 0.4;
    ModuleMD::SNAP::Neighbor second = first;
    second.x = -0.3;
    second.y = 1.7;
    second.z = -0.2;
    const std::vector<ModuleMD::SNAP::Neighbor> neighbors = {first, second};
    const ModuleMD::SNAP::ForceResult analytic = calculator.compute_with_analytic_force(neighbors, 0);
    const ModuleMD::SNAP::ForceResult finite = calculator.compute_with_finite_difference_force(neighbors, 0, 1.0e-6);
    for (std::size_t i = 0; i < neighbors.size(); ++i)
    {
        for (int direction = 0; direction < 3; ++direction)
        {
            EXPECT_NEAR(analytic.neighbor_forces[i][static_cast<std::size_t>(direction)],
                        finite.neighbor_forces[i][static_cast<std::size_t>(direction)], 1.0e-5);
        }
    }

    model.parameters.quadratic_flag = true;
    model.coefficients[0].resize(21, 0.0);
    model.coefficients[0][6] = 0.25;
    model.coefficients[0][10] = -0.15;
    const ModuleMD::SNAP::Calculator quadratic(model);
    const ModuleMD::SNAP::ForceResult quadratic_analytic = quadratic.compute_with_analytic_force(neighbors, 0);
    const ModuleMD::SNAP::ForceResult quadratic_finite = quadratic.compute_with_finite_difference_force(neighbors, 0, 1.0e-6);
    for (std::size_t i = 0; i < neighbors.size(); ++i)
    {
        for (int direction = 0; direction < 3; ++direction)
        {
            EXPECT_NEAR(quadratic_analytic.neighbor_forces[i][static_cast<std::size_t>(direction)],
                        quadratic_finite.neighbor_forces[i][static_cast<std::size_t>(direction)], 1.0e-5);
        }
    }
}

TEST(SnapCalculator, InnerSwitchSkipsSingularRadius)
{
    ModuleMD::SNAP::ModelData model = make_model();
    model.parameters.rmin = 0.5;
    model.parameters.switchinner_flag = true;
    model.parameters.sinner.push_back(1.5);
    model.parameters.dinner.push_back(0.4);
    const ModuleMD::SNAP::Calculator calculator(model);
    ModuleMD::SNAP::Neighbor neighbor;
    neighbor.x = 0.5;
    neighbor.sinner = 1.5;
    neighbor.dinner = 0.4;
    const std::vector<ModuleMD::SNAP::Neighbor> neighbors(1, neighbor);
    EXPECT_NO_THROW(calculator.compute_with_analytic_force(neighbors, 0));
    model.parameters.switchinner_flag = false;
    const ModuleMD::SNAP::Calculator bare(model);
    EXPECT_THROW(bare.compute(neighbors, 0), std::domain_error);
    EXPECT_THROW(bare.compute_with_analytic_force(neighbors, 0), std::domain_error);

    model.parameters.twojmax = 2;
    model.coefficients[0] = {0.2, 0.3, -0.1, 0.5, 0.2, -0.4};
    const ModuleMD::SNAP::Calculator angular(model);
    neighbor.x = 0.55;
    neighbor.y = 0.12;
    const std::vector<ModuleMD::SNAP::Neighbor> near_neighbors(1, neighbor);
    const ModuleMD::SNAP::ForceResult near_analytic = angular.compute_with_analytic_force(near_neighbors, 0);
    const ModuleMD::SNAP::ForceResult near_finite = angular.compute_with_finite_difference_force(near_neighbors, 0, 1.0e-6);
    EXPECT_TRUE(std::isfinite(near_analytic.energy.energy));
    for (int direction = 0; direction < 3; ++direction)
    {
        EXPECT_TRUE(std::isfinite(near_analytic.neighbor_forces[0][static_cast<std::size_t>(direction)]));
        EXPECT_NEAR(near_analytic.neighbor_forces[0][static_cast<std::size_t>(direction)],
                    near_finite.neighbor_forces[0][static_cast<std::size_t>(direction)], 1.0e-4);
    }
}

TEST(SnapCalculator, ZeroWeightSkipsInvalidRadius)
{
    ModuleMD::SNAP::ModelData model = make_model();
    model.parameters.rmin = 0.5;
    const ModuleMD::SNAP::Calculator calculator(model);
    ModuleMD::SNAP::Neighbor neighbor;
    neighbor.x = 0.4;
    neighbor.weight = 0.0;
    const std::vector<ModuleMD::SNAP::Neighbor> neighbors(1, neighbor);
    const std::vector<ModuleMD::SNAP::Neighbor> empty;
    const double isolated_energy = calculator.compute(empty, 0).energy;
    const ModuleMD::SNAP::ForceResult result = calculator.compute_with_analytic_force(neighbors, 0);
    EXPECT_DOUBLE_EQ(result.energy.energy, isolated_energy);
    for (int direction = 0; direction < 3; ++direction)
    {
        EXPECT_DOUBLE_EQ(result.center_force[static_cast<std::size_t>(direction)], 0.0);
        EXPECT_DOUBLE_EQ(result.neighbor_forces[0][static_cast<std::size_t>(direction)], 0.0);
    }
    neighbor.weight = 1.0;
    const std::vector<ModuleMD::SNAP::Neighbor> active(1, neighbor);
    EXPECT_THROW(calculator.compute(active, 0), std::domain_error);
    EXPECT_THROW(calculator.compute_with_analytic_force(active, 0), std::domain_error);
}

TEST(SnapCalculator, AngularQuadraticForcesMatchLammpsAndStrain)
{
    const ModuleMD::SNAP::Calculator calculator(make_angular_quadratic_model());
    const std::array<std::array<double, 3>, 3> positions = {{{0.0, 0.0, 0.0},
                                                              {1.058354, 0.0, 0.0},
                                                              {0.2645885, 1.3229425, 0.0}}};
    const SmallCellResult result = evaluate_small_cell(calculator, positions);
    const double lammps_forces[3][2] = {{67.2928038408958, 57.5476436778783},
                                         {-92.1399449575303, 60.5944497536836},
                                         {24.8471411166345, -118.142093431562}};
    EXPECT_NEAR(result.energy, 54.909561, 1.0e-5);
    for (std::size_t atom = 0; atom < positions.size(); ++atom)
    {
        for (int direction = 0; direction < 2; ++direction)
        {
            EXPECT_NEAR(result.forces[atom][static_cast<std::size_t>(direction)],
                        lammps_forces[atom][direction], 2.0e-4);
        }
    }

    const double strain = 1.0e-6;
    std::array<std::array<double, 3>, 3> expanded = positions;
    std::array<std::array<double, 3>, 3> contracted = positions;
    for (std::size_t atom = 0; atom < positions.size(); ++atom)
    {
        expanded[atom][0] *= 1.0 + strain;
        contracted[atom][0] *= 1.0 - strain;
    }
    const double expanded_energy = evaluate_small_cell(calculator, expanded).energy;
    const double contracted_energy = evaluate_small_cell(calculator, contracted).energy;
    const double strain_derivative = (expanded_energy - contracted_energy) / (2.0 * strain);
    EXPECT_NEAR(result.virial_xx, -strain_derivative, 1.0e-3);
}

} // namespace
