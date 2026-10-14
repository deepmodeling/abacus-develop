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
