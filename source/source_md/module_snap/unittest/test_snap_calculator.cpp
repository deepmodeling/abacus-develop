#include "../snap_calculator.h"

#include <gtest/gtest.h>

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

TEST(SnapCalculator, RejectsMismatchedCoefficientCount)
{
    ModuleMD::SNAP::ModelData model = make_model();
    model.coefficients[0].push_back(4.0);
    EXPECT_THROW(ModuleMD::SNAP::Calculator calculator(model), std::invalid_argument);
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

} // namespace
