#include "../xc_ncgga_radial.h"

#include "gtest/gtest.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace
{

double dot(const std::array<double, 3>& left, const std::array<double, 3>& right)
{
    double result = 0.0;
    for (int component = 0; component < 3; ++component)
    {
        result += left[component] * right[component];
    }
    return result;
}

std::array<double, 3> shifted(const std::array<double, 3>& point,
                              const std::array<double, 3>& direction,
                              const double step)
{
    std::array<double, 3> result = point;
    for (int component = 0; component < 3; ++component)
    {
        result[component] += step * direction[component];
    }
    return result;
}

double directional_hessian(const ModuleXC::NcggaRadialPoint& point, const std::array<double, 3>& direction)
{
    double result = 0.0;
    for (int row = 0; row < 3; ++row)
    {
        for (int column = 0; column < 3; ++column)
        {
            result += direction[row] * point.jacobian(row, column) * direction[column];
        }
    }
    return result;
}

ModuleXC::NcggaSpinMapPoint make_spin_map(const std::array<double, 4>& state, const double eta)
{
    const std::array<double, 3> magnetization = {{state[1], state[2], state[3]}};
    return ModuleXC::make_ncgga_spin_map_point(state[0], ModuleXC::make_ncgga_radial_point(magnetization, eta));
}

TEST(NcggaRadial, RejectsNonPositiveEta)
{
    const std::array<double, 3> magnetization = {{0.1, 0.0, 0.0}};
    EXPECT_THROW(ModuleXC::make_ncgga_radial_point(magnetization, 0.0), std::invalid_argument);
    EXPECT_THROW(ModuleXC::make_ncgga_radial_point(magnetization, -0.5), std::invalid_argument);
}

TEST(NcggaRadial, MatchesZeroInteriorJoinAndExteriorValues)
{
    const double eta = 0.5;
    const ModuleXC::NcggaRadialPoint zero = ModuleXC::make_ncgga_radial_point({{0.0, 0.0, 0.0}}, eta);
    EXPECT_DOUBLE_EQ(zero.value, 0.0);
    for (int row = 0; row < 3; ++row)
    {
        EXPECT_DOUBLE_EQ(zero.direction[row], 0.0);
        EXPECT_DOUBLE_EQ(zero.gradient[row], 0.0);
        for (int column = 0; column < 3; ++column)
        {
            EXPECT_DOUBLE_EQ(zero.jacobian(row, column), 0.0);
        }
    }

    const ModuleXC::NcggaRadialPoint interior = ModuleXC::make_ncgga_radial_point({{0.5 * eta, 0.0, 0.0}}, eta);
    EXPECT_DOUBLE_EQ(interior.value, 11.0 * eta / 32.0);
    EXPECT_DOUBLE_EQ(interior.gradient[0], 23.0 / 16.0);
    EXPECT_DOUBLE_EQ(interior.jacobian(0, 0), 3.0 / (2.0 * eta));
    EXPECT_DOUBLE_EQ(interior.jacobian(1, 1), 23.0 / (8.0 * eta));
    EXPECT_DOUBLE_EQ(interior.jacobian(2, 2), 23.0 / (8.0 * eta));

    const ModuleXC::NcggaRadialPoint join = ModuleXC::make_ncgga_radial_point({{eta, 0.0, 0.0}}, eta);
    EXPECT_DOUBLE_EQ(join.value, eta);
    EXPECT_DOUBLE_EQ(join.gradient[0], 1.0);
    EXPECT_DOUBLE_EQ(join.jacobian(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(join.jacobian(1, 1), 1.0 / eta);
    EXPECT_DOUBLE_EQ(join.jacobian(2, 2), 1.0 / eta);

    const ModuleXC::NcggaRadialPoint exterior = ModuleXC::make_ncgga_radial_point({{3.0, 4.0, 0.0}}, eta);
    EXPECT_DOUBLE_EQ(exterior.value, 5.0);
    EXPECT_DOUBLE_EQ(exterior.gradient[0], 0.6);
    EXPECT_DOUBLE_EQ(exterior.gradient[1], 0.8);
    EXPECT_NEAR(exterior.jacobian(0, 0), 0.128, 1.0e-15);
    EXPECT_NEAR(exterior.jacobian(0, 1), -0.096, 1.0e-15);
    EXPECT_NEAR(exterior.jacobian(1, 1), 0.072, 1.0e-15);
    EXPECT_NEAR(exterior.jacobian(2, 2), 0.2, 1.0e-15);
}

TEST(NcggaRadial, CompactBranchConvergesC2AtJoin)
{
    const double eta = 0.5;
    const double coarse_delta = 1.0e-3;
    const double fine_delta = 0.5 * coarse_delta;
    const ModuleXC::NcggaRadialPoint coarse
        = ModuleXC::make_ncgga_radial_point({{eta * (1.0 - coarse_delta), 0.0, 0.0}}, eta);
    const ModuleXC::NcggaRadialPoint fine
        = ModuleXC::make_ncgga_radial_point({{eta * (1.0 - fine_delta), 0.0, 0.0}}, eta);

    EXPECT_LE(std::abs(fine.value - eta), 0.51 * std::abs(coarse.value - eta));
    EXPECT_LE(std::abs(fine.gradient[0] - 1.0), 0.26 * std::abs(coarse.gradient[0] - 1.0));
    EXPECT_LE(std::abs(fine.jacobian(0, 0)), 0.51 * std::abs(coarse.jacobian(0, 0)));
    EXPECT_LE(std::abs(fine.jacobian(1, 1) - 1.0 / eta), 0.51 * std::abs(coarse.jacobian(1, 1) - 1.0 / eta));
}

TEST(NcggaRadial, CompactBranchConvergesC2AtOrigin)
{
    const double eta = 0.5;
    const double coarse_radius = 1.0e-3 * eta;
    const double fine_radius = 0.5 * coarse_radius;
    const ModuleXC::NcggaRadialPoint coarse = ModuleXC::make_ncgga_radial_point({{coarse_radius, 0.0, 0.0}}, eta);
    const ModuleXC::NcggaRadialPoint fine = ModuleXC::make_ncgga_radial_point({{fine_radius, 0.0, 0.0}}, eta);

    EXPECT_LE(fine.value, 0.13 * coarse.value);
    EXPECT_LE(std::abs(fine.gradient[0]), 0.26 * std::abs(coarse.gradient[0]));
    EXPECT_LE(std::abs(fine.jacobian(0, 0)), 0.51 * std::abs(coarse.jacobian(0, 0)));
    EXPECT_LE(std::abs(fine.jacobian(1, 1)), 0.51 * std::abs(coarse.jacobian(1, 1)));
}

TEST(NcggaRadial, GradientMatchesCentralFiniteDifference)
{
    const double eta = 0.5;
    const double step = 2.0e-6 * eta;
    const std::array<double, 3> direction = {{0.31, -0.27, 0.19}};
    const std::array<double, 3> interior = {{0.12, -0.08, 0.05}};
    const std::array<double, 3> exterior = {{0.62, -0.41, 0.23}};
    const std::array<double, 3> join = {{eta, 0.0, 0.0}};
    const std::array<std::array<double, 3>, 3> points = {{interior, exterior, join}};

    for (std::size_t sample = 0; sample < points.size(); ++sample)
    {
        const ModuleXC::NcggaRadialPoint point = ModuleXC::make_ncgga_radial_point(points[sample], eta);
        const double finite_difference
            = (ModuleXC::make_ncgga_radial_point(shifted(points[sample], direction, step), eta).value
               - ModuleXC::make_ncgga_radial_point(shifted(points[sample], direction, -step), eta).value)
              / (2.0 * step);
        EXPECT_NEAR(finite_difference, dot(point.gradient, direction), 2.0e-10);
    }
}

TEST(NcggaRadial, JacobianMatchesGradientFiniteDifference)
{
    const double eta = 0.5;
    const double step = 2.0e-6 * eta;
    const std::array<double, 3> direction = {{0.31, -0.27, 0.19}};
    const std::array<double, 3> interior = {{0.12, -0.08, 0.05}};
    const std::array<double, 3> exterior = {{0.62, -0.41, 0.23}};
    const std::array<std::array<double, 3>, 2> points = {{interior, exterior}};

    for (std::size_t sample = 0; sample < points.size(); ++sample)
    {
        const ModuleXC::NcggaRadialPoint point = ModuleXC::make_ncgga_radial_point(points[sample], eta);
        const ModuleXC::NcggaRadialPoint plus
            = ModuleXC::make_ncgga_radial_point(shifted(points[sample], direction, step), eta);
        const ModuleXC::NcggaRadialPoint minus
            = ModuleXC::make_ncgga_radial_point(shifted(points[sample], direction, -step), eta);
        for (int row = 0; row < 3; ++row)
        {
            double analytic = 0.0;
            for (int column = 0; column < 3; ++column)
            {
                analytic += point.jacobian(row, column) * direction[column];
                EXPECT_NEAR(point.jacobian(row, column), point.jacobian(column, row), 1.0e-15);
            }
            const double finite_difference = (plus.gradient[row] - minus.gradient[row]) / (2.0 * step);
            EXPECT_NEAR(finite_difference, analytic, 2.0e-9);
        }
    }
}

TEST(NcggaRadial, DirectionalSecondDerivativeMatchesValueFiniteDifference)
{
    const double eta = 0.5;
    const double step = 2.0e-4 * eta;
    const std::array<double, 3> direction = {{0.31, -0.27, 0.19}};
    const std::array<double, 3> interior = {{0.12, -0.08, 0.05}};
    const std::array<double, 3> exterior = {{0.62, -0.41, 0.23}};
    const std::array<std::array<double, 3>, 2> points = {{interior, exterior}};

    for (std::size_t sample = 0; sample < points.size(); ++sample)
    {
        const ModuleXC::NcggaRadialPoint point = ModuleXC::make_ncgga_radial_point(points[sample], eta);
        const double plus = ModuleXC::make_ncgga_radial_point(shifted(points[sample], direction, step), eta).value;
        const double minus = ModuleXC::make_ncgga_radial_point(shifted(points[sample], direction, -step), eta).value;
        const double finite_difference = (plus - 2.0 * point.value + minus) / (step * step);
        EXPECT_NEAR(finite_difference, directional_hessian(point, direction), 2.0e-7);
    }
}

TEST(NcggaRadial, SpinDensityMapJacobianMatchesFiniteDifferenceInEachBranch)
{
    const double eta = 0.5;
    const double step = 2.0e-7;
    const std::array<std::array<double, 4>, 5> states = {{{{1.20, 0.30, -0.20, 0.10}},
                                                          {{-1.20, 0.30, -0.20, 0.10}},
                                                          {{0.20, 0.42, -0.31, 0.16}},
                                                          {{-0.20, 0.42, -0.31, 0.16}},
                                                          {{0.90, 0.12, -0.08, 0.05}}}};

    for (std::size_t sample = 0; sample < states.size(); ++sample)
    {
        const ModuleXC::NcggaSpinMapPoint point = make_spin_map(states[sample], eta);
        EXPECT_GE(point.spin_density[0], 0.0);
        EXPECT_GE(point.spin_density[1], 0.0);
        EXPECT_NEAR(point.spin_density[0] + point.spin_density[1], std::abs(states[sample][0]), 1.0e-15);
        EXPECT_NEAR(point.spin_density[0] - point.spin_density[1], point.clipped_magnitude, 1.0e-15);

        double maximum_error = 0.0;
        for (int channel = 0; channel < 4; ++channel)
        {
            std::array<double, 4> plus_state = states[sample];
            std::array<double, 4> minus_state = states[sample];
            plus_state[channel] += step;
            minus_state[channel] -= step;
            const ModuleXC::NcggaSpinMapPoint plus = make_spin_map(plus_state, eta);
            const ModuleXC::NcggaSpinMapPoint minus = make_spin_map(minus_state, eta);
            for (int spin = 0; spin < 2; ++spin)
            {
                const double finite_difference = (plus.spin_density[spin] - minus.spin_density[spin]) / (2.0 * step);
                maximum_error = std::max(maximum_error, std::abs(finite_difference - point.jacobian(spin, channel)));
                EXPECT_NEAR(finite_difference, point.jacobian(spin, channel), 3.0e-9)
                    << "sample=" << sample << " spin=" << spin << " channel=" << channel;
            }
        }
        std::cout << std::setprecision(17) << "NCGGA_SPIN_MAP_FD sample=" << sample
                  << " total_density=" << states[sample][0] << " absolute_density=" << point.absolute_density
                  << " radial_value=" << point.radial.value << " clipped_magnitude=" << point.clipped_magnitude
                  << " saturated=" << point.saturated << " step=" << step << " maximum_error=" << maximum_error << '\n';
    }
}

TEST(NcggaRadial, SpinDensityMapDefinesAbsAndSaturationKinkConventions)
{
    const double eta = 0.1;
    const ModuleXC::NcggaRadialPoint nonzero_radial = ModuleXC::make_ncgga_radial_point({{0.30, 0.0, 0.0}}, eta);
    const ModuleXC::NcggaSpinMapPoint zero_density = ModuleXC::make_ncgga_spin_map_point(0.0, nonzero_radial);
    EXPECT_TRUE(zero_density.saturated);
    EXPECT_DOUBLE_EQ(zero_density.spin_density[0], 0.0);
    EXPECT_DOUBLE_EQ(zero_density.spin_density[1], 0.0);
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int channel = 0; channel < 4; ++channel)
        {
            EXPECT_DOUBLE_EQ(zero_density.jacobian(spin, channel), 0.0);
        }
    }

    const ModuleXC::NcggaRadialPoint equality_radial = ModuleXC::make_ncgga_radial_point({{0.50, 0.0, 0.0}}, eta);
    const ModuleXC::NcggaSpinMapPoint positive = ModuleXC::make_ncgga_spin_map_point(0.50, equality_radial);
    const ModuleXC::NcggaSpinMapPoint negative = ModuleXC::make_ncgga_spin_map_point(-0.50, equality_radial);
    EXPECT_TRUE(positive.saturated);
    EXPECT_TRUE(negative.saturated);
    EXPECT_DOUBLE_EQ(positive.jacobian(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(positive.jacobian(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(negative.jacobian(0, 0), -1.0);
    EXPECT_DOUBLE_EQ(negative.jacobian(1, 0), 0.0);
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int channel = 1; channel < 4; ++channel)
        {
            EXPECT_DOUBLE_EQ(positive.jacobian(spin, channel), 0.0);
            EXPECT_DOUBLE_EQ(negative.jacobian(spin, channel), 0.0);
        }
    }
}

} // namespace
