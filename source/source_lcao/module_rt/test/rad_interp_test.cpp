#include "source_lcao/module_rt/radial_interpolation.h"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace
{

constexpr double interpolation_tolerance = 5.0e-15;

double cubic(const double radius)
{
    return ((0.75 * radius - 1.25) * radius + 0.5) * radius + 2.0;
}

void expect_cubic_interpolation(const std::vector<double>& grid, const std::vector<double>& queries)
{
    std::vector<double> values(grid.size());
    for (size_t i = 0; i < grid.size(); ++i)
    {
        values[i] = cubic(grid[i]);
    }

    module_rt::RadialGridInfo grid_info;
    ASSERT_TRUE(module_rt::analyze_radial_grid(grid.data(),
                                               values.data(),
                                               static_cast<int>(grid.size()),
                                               grid_info));
    double max_error = 0.0;
    for (const double query: queries)
    {
        const double interpolated = module_rt::interpolate_radial(grid.data(),
                                                                  values.data(),
                                                                  static_cast<int>(grid.size()),
                                                                  grid_info,
                                                                  query);
        const double error = std::abs(interpolated - cubic(query));
        max_error = std::max(max_error, error);
        EXPECT_LE(error, interpolation_tolerance) << "query = " << query;
    }
    std::cout << std::scientific << std::setprecision(12) << "radial interpolation mesh=" << grid.size()
              << ", r=[" << grid.front() << ", " << grid.back() << "]: max cubic error = " << max_error
              << std::defaultfloat << std::endl;
}

} // namespace

TEST(RadialInterpolation, UniformGridCoversTailIntervalsAndEndpoint)
{
    const std::vector<double> grid = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25};
    const std::vector<double> queries = {2.1, 2.6, 3.1, 3.25};
    expect_cubic_interpolation(grid, queries);

    std::vector<double> values(grid.size());
    module_rt::RadialGridInfo grid_info;
    EXPECT_TRUE(module_rt::analyze_radial_grid(grid.data(),
                                               values.data(),
                                               static_cast<int>(grid.size()),
                                               grid_info));
    EXPECT_TRUE(grid_info.is_uniform);
}

TEST(RadialInterpolation, NonuniformGridsReproduceCubic)
{
    std::vector<double> logarithmic_grid;
    std::vector<double> shifted_logarithmic_grid;
    for (int i = 0; i < 9; ++i)
    {
        logarithmic_grid.push_back(std::exp(-2.0 + 0.35 * i));
        shifted_logarithmic_grid.push_back(0.4 + std::exp(-3.0 + 0.5 * i));
    }
    const std::vector<double> strongly_nonuniform_grid = {0.2, 0.201, 0.35, 1.4, 1.45, 4.0, 9.0};

    expect_cubic_interpolation(logarithmic_grid,
                               {logarithmic_grid[1] * 1.1,
                                0.5 * (logarithmic_grid[5] + logarithmic_grid[6]),
                                0.5 * (logarithmic_grid[7] + logarithmic_grid[8])});
    expect_cubic_interpolation(shifted_logarithmic_grid,
                               {0.5 * (shifted_logarithmic_grid[0] + shifted_logarithmic_grid[1]),
                                0.5 * (shifted_logarithmic_grid[6] + shifted_logarithmic_grid[7]),
                                shifted_logarithmic_grid.back()});
    expect_cubic_interpolation(strongly_nonuniform_grid, {0.2005, 1.43, 3.2, 8.5});

    std::vector<double> values(strongly_nonuniform_grid.size());
    module_rt::RadialGridInfo grid_info;
    ASSERT_TRUE(module_rt::analyze_radial_grid(strongly_nonuniform_grid.data(),
                                               values.data(),
                                               static_cast<int>(strongly_nonuniform_grid.size()),
                                               grid_info));
    EXPECT_FALSE(grid_info.is_uniform);
}

TEST(RadialInterpolation, PreservesGridValuesAndRejectsExtrapolation)
{
    const std::vector<double> grid = {0.3, 0.5, 1.1, 2.0, 3.7};
    const std::vector<double> values = {9.0, -1.0, 4.5, 8.0, -3.0};
    module_rt::RadialGridInfo grid_info;
    ASSERT_TRUE(module_rt::analyze_radial_grid(grid.data(),
                                               values.data(),
                                               static_cast<int>(grid.size()),
                                               grid_info));

    for (size_t i = 0; i < grid.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(module_rt::interpolate_radial(grid.data(),
                                                      values.data(),
                                                      static_cast<int>(grid.size()),
                                                      grid_info,
                                                      grid[i]),
                         values[i]);
    }

    EXPECT_DOUBLE_EQ(module_rt::interpolate_radial(grid.data(), values.data(), 5, grid_info, 0.2), 0.0);
    EXPECT_DOUBLE_EQ(module_rt::interpolate_radial(grid.data(), values.data(), 5, grid_info, 4.0), 0.0);
    EXPECT_DOUBLE_EQ(module_rt::interpolate_radial(grid.data(),
                                                  values.data(),
                                                  5,
                                                  grid_info,
                                                  std::numeric_limits<double>::quiet_NaN()),
                     0.0);
    EXPECT_DOUBLE_EQ(module_rt::interpolate_radial(grid.data(),
                                                  values.data(),
                                                  5,
                                                  grid_info,
                                                  std::numeric_limits<double>::infinity()),
                     0.0);
}

TEST(RadialInterpolation, HandlesSmallMeshes)
{
    double max_small_mesh_error = 0.0;
    {
        const std::vector<double> grid = {0.4};
        const std::vector<double> values = {3.0};
        module_rt::RadialGridInfo grid_info;
        ASSERT_TRUE(module_rt::analyze_radial_grid(grid.data(), values.data(), 1, grid_info));
        EXPECT_DOUBLE_EQ(module_rt::interpolate_radial(grid.data(), values.data(), 1, grid_info, 0.4), 3.0);
        EXPECT_DOUBLE_EQ(module_rt::interpolate_radial(grid.data(), values.data(), 1, grid_info, 0.5), 0.0);
    }

    {
        const std::vector<double> grid = {0.4, 1.2};
        const std::vector<double> values = {-0.2, 1.4};
        module_rt::RadialGridInfo grid_info;
        ASSERT_TRUE(module_rt::analyze_radial_grid(grid.data(), values.data(), 2, grid_info));
        const double error
            = std::abs(module_rt::interpolate_radial(grid.data(), values.data(), 2, grid_info, 0.8) - 0.6);
        max_small_mesh_error = std::max(max_small_mesh_error, error);
        EXPECT_LE(error, interpolation_tolerance);
    }

    {
        const std::vector<double> grid = {0.4, 0.7, 1.2};
        std::vector<double> values;
        for (const double radius: grid)
        {
            values.push_back(radius * radius - 0.5 * radius + 1.0);
        }
        module_rt::RadialGridInfo grid_info;
        ASSERT_TRUE(module_rt::analyze_radial_grid(grid.data(), values.data(), 3, grid_info));
        for (const double query: {0.4, 0.55, 0.9, 1.2})
        {
            const double expected = query * query - 0.5 * query + 1.0;
            const double error
                = std::abs(module_rt::interpolate_radial(grid.data(), values.data(), 3, grid_info, query) - expected);
            max_small_mesh_error = std::max(max_small_mesh_error, error);
            EXPECT_LE(error, interpolation_tolerance) << "query = " << query;
        }
    }

    expect_cubic_interpolation({0.4, 0.7, 1.2, 2.0}, {0.4, 0.55, 0.9, 1.6, 2.0});
    std::cout << std::scientific << std::setprecision(12)
              << "radial interpolation mesh=1/2/3: max polynomial error = " << max_small_mesh_error
              << std::defaultfloat << std::endl;
}

TEST(RadialInterpolation, RejectsInvalidRadialData)
{
    const double values[] = {1.0, 2.0, 3.0};
    module_rt::RadialGridInfo grid_info;

    EXPECT_FALSE(module_rt::analyze_radial_grid(nullptr, values, 3, grid_info));
    EXPECT_FALSE(module_rt::analyze_radial_grid(values, nullptr, 3, grid_info));
    EXPECT_FALSE(module_rt::analyze_radial_grid(values, values, 0, grid_info));

    const double duplicate_grid[] = {0.0, 0.5, 0.5};
    const double reversed_grid[] = {0.0, 0.5, 0.4};
    const double nonfinite_grid[] = {0.0, std::numeric_limits<double>::infinity(), 1.0};
    const double nonfinite_values[] = {1.0, std::numeric_limits<double>::quiet_NaN(), 3.0};

    EXPECT_FALSE(module_rt::analyze_radial_grid(duplicate_grid, values, 3, grid_info));
    EXPECT_FALSE(module_rt::analyze_radial_grid(reversed_grid, values, 3, grid_info));
    EXPECT_FALSE(module_rt::analyze_radial_grid(nonfinite_grid, values, 3, grid_info));
    EXPECT_FALSE(module_rt::analyze_radial_grid(values, nonfinite_values, 3, grid_info));
}
