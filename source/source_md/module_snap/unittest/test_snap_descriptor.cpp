#include "../snap_descriptor.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

namespace
{

ModuleMD::SNAP::ModelParameters make_parameters()
{
    ModuleMD::SNAP::ModelParameters parameters;
    parameters.rcut = 4.5;
    parameters.rcutfac = 1.0;
    parameters.rmin = 0.0;
    parameters.rfac0 = 0.99363;
    parameters.twojmax = 2;
    parameters.bzero_flag = false;
    parameters.elements.push_back("Si");
    parameters.element_weights.push_back(1.0);
    parameters.element_radii.push_back(2.25);
    return parameters;
}

TEST(SnapDescriptor, AllocatesExpectedArraySizes)
{
    const ModuleMD::SNAP::DescriptorEvaluator evaluator(make_parameters());
    EXPECT_EQ(evaluator.u_size(), static_cast<std::size_t>(1 * 1 + 2 * 2 + 3 * 3));
    EXPECT_EQ(evaluator.bispectrum_size(), static_cast<std::size_t>(5));
}

TEST(SnapDescriptor, SelfContributionIsDiagonal)
{
    const ModuleMD::SNAP::DescriptorEvaluator evaluator(make_parameters());
    const std::vector<ModuleMD::SNAP::Neighbor> neighbors;
    const std::vector<std::complex<double> > u = evaluator.accumulate_u(neighbors, 0);

    EXPECT_DOUBLE_EQ(u[0].real(), 1.0);
    EXPECT_DOUBLE_EQ(u[0].imag(), 0.0);
    EXPECT_DOUBLE_EQ(u[4].real(), 1.0);
    EXPECT_DOUBLE_EQ(u[4].imag(), 0.0);
    EXPECT_DOUBLE_EQ(u[13].real(), 1.0);
    EXPECT_DOUBLE_EQ(u[13].imag(), 0.0);
}

TEST(SnapDescriptor, OutsideCutoffDoesNotChangeEnvironment)
{
    const ModuleMD::SNAP::DescriptorEvaluator evaluator(make_parameters());
    ModuleMD::SNAP::Neighbor neighbor;
    neighbor.x = 5.0;
    const std::vector<ModuleMD::SNAP::Neighbor> neighbors(1, neighbor);
    const std::vector<std::complex<double> > u = evaluator.accumulate_u(neighbors, 0);
    const std::vector<ModuleMD::SNAP::Neighbor> no_neighbors;
    const std::vector<std::complex<double> > self = evaluator.accumulate_u(no_neighbors, 0);

    ASSERT_EQ(u.size(), self.size());
    for (std::size_t i = 0; i < u.size(); ++i)
    {
        EXPECT_NEAR(u[i].real(), self[i].real(), 1.0e-14);
        EXPECT_NEAR(u[i].imag(), self[i].imag(), 1.0e-14);
    }
}

TEST(SnapDescriptor, BispectrumIsFinite)
{
    const ModuleMD::SNAP::DescriptorEvaluator evaluator(make_parameters());
    ModuleMD::SNAP::Neighbor neighbor;
    neighbor.x = 1.2;
    neighbor.y = 0.4;
    neighbor.z = -0.8;
    const std::vector<ModuleMD::SNAP::Neighbor> neighbors(1, neighbor);
    const std::vector<std::complex<double> > u = evaluator.accumulate_u(neighbors, 0);
    const std::vector<double> bispectrum = evaluator.compute_bispectrum(u);

    for (std::size_t i = 0; i < bispectrum.size(); ++i)
    {
        EXPECT_TRUE(std::isfinite(bispectrum[i]));
    }
}

TEST(SnapDescriptor, SmallNonzeroDistanceIsNotDiscarded)
{
    const ModuleMD::SNAP::DescriptorEvaluator evaluator(make_parameters());
    // This is smaller than the former 1e-14 Angstrom skip threshold.
    const double small_distance = 1.0e-15;
    ModuleMD::SNAP::Neighbor neighbor;
    neighbor.x = small_distance;
    const std::vector<std::complex<double> > u = evaluator.compute_u(neighbor);
    EXPECT_DOUBLE_EQ(u[0].real(), 1.0);

    const ModuleMD::SNAP::UWithDerivative result = evaluator.compute_u_with_derivative(neighbor, 0);
    double largest_derivative = 0.0;
    for (std::size_t i = 0; i < result.derivatives.size(); ++i)
    {
        largest_derivative = std::max(largest_derivative, std::abs(result.derivatives[i]));
    }
    EXPECT_GT(largest_derivative, 0.1);
}

TEST(SnapDescriptor, MatchesLammpsSingleElementReference)
{
    const ModuleMD::SNAP::DescriptorEvaluator evaluator(make_parameters());
    const double positions[3][3] = {{1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {1.0, 3.0, 1.0}};
    std::vector<ModuleMD::SNAP::Neighbor> neighbors;
    for (int atom = 1; atom < 3; ++atom)
    {
        ModuleMD::SNAP::Neighbor neighbor;
        neighbor.x = positions[atom][0] - positions[0][0];
        neighbor.y = positions[atom][1] - positions[0][1];
        neighbor.z = positions[atom][2] - positions[0][2];
        neighbor.cutoff_radius = 4.0;
        neighbors.push_back(neighbor);
    }
    const std::vector<std::complex<double> > u = evaluator.accumulate_u(neighbors, 0);
    const std::vector<double> bispectrum = evaluator.compute_bispectrum(u);
    const double reference[5] = {
        13.0368348609, 15.0992648922, 11.7631725456, 13.7069739819, 6.07477849248};
    ASSERT_EQ(bispectrum.size(), static_cast<std::size_t>(5));
    for (int index = 0; index < 5; ++index)
    {
        EXPECT_NEAR(bispectrum[static_cast<std::size_t>(index)], reference[index], 1.0e-9);
    }
}

} // namespace
