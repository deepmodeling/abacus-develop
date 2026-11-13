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

TEST(SnapDescriptor, OrderZeroBispectrumHasAnalyticValue)
{
    ModuleMD::SNAP::ModelParameters parameters = make_parameters();
    parameters.twojmax = 0;
    parameters.switch_flag = false;
    const ModuleMD::SNAP::DescriptorEvaluator evaluator(parameters);
    ModuleMD::SNAP::Neighbor first;
    first.x = 1.0;
    first.weight = 0.7;
    ModuleMD::SNAP::Neighbor second;
    second.y = 2.0;
    second.weight = 1.3;
    const std::vector<ModuleMD::SNAP::Neighbor> neighbors = {first, second};
    const std::vector<std::complex<double> > u = evaluator.accumulate_u(neighbors, 0);
    const std::vector<double> bispectrum = evaluator.compute_bispectrum(u);
    ASSERT_EQ(bispectrum.size(), static_cast<std::size_t>(1));
    EXPECT_NEAR(bispectrum[0], 27.0, 1.0e-12);
}

TEST(SnapDescriptor, BispectrumIsInvariantUnderRotationAndNeighborOrder)
{
    const ModuleMD::SNAP::DescriptorEvaluator evaluator(make_parameters());
    ModuleMD::SNAP::Neighbor first;
    first.x = 1.0;
    first.y = 0.3;
    first.z = -0.2;
    ModuleMD::SNAP::Neighbor second;
    second.x = -0.4;
    second.y = 1.7;
    second.z = 0.5;
    const std::vector<ModuleMD::SNAP::Neighbor> neighbors = {first, second};
    const std::vector<std::complex<double> > u = evaluator.accumulate_u(neighbors, 0);
    const std::vector<double> bispectrum = evaluator.compute_bispectrum(u);

    const std::vector<ModuleMD::SNAP::Neighbor> reordered = {second, first};
    const std::vector<std::complex<double> > reordered_u = evaluator.accumulate_u(reordered, 0);
    const std::vector<double> reordered_bispectrum = evaluator.compute_bispectrum(reordered_u);
    ASSERT_EQ(bispectrum.size(), reordered_bispectrum.size());

    std::vector<ModuleMD::SNAP::Neighbor> rotated = neighbors;
    for (ModuleMD::SNAP::Neighbor& neighbor : rotated)
    {
        const double old_x = neighbor.x;
        neighbor.x = -neighbor.y;
        neighbor.y = old_x;
    }
    const std::vector<std::complex<double> > rotated_u = evaluator.accumulate_u(rotated, 0);
    const std::vector<double> rotated_bispectrum = evaluator.compute_bispectrum(rotated_u);
    ASSERT_EQ(bispectrum.size(), rotated_bispectrum.size());
    for (std::size_t index = 0; index < bispectrum.size(); ++index)
    {
        EXPECT_NEAR(bispectrum[index], reordered_bispectrum[index], 1.0e-10);
        EXPECT_NEAR(bispectrum[index], rotated_bispectrum[index], 1.0e-10);
    }
}

} // namespace
