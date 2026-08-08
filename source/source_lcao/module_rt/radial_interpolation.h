#ifndef MODULE_RT_RADIAL_INTERPOLATION_H
#define MODULE_RT_RADIAL_INTERPOLATION_H

#include <cmath>

namespace module_rt
{

/**
 * @brief Metadata used to accelerate interpolation on uniform radial grids.
 */
struct RadialGridInfo
{
    double r_min = 0.0;
    double r_max = 0.0;
    double inv_spacing = 0.0;
    bool is_uniform = false;
};

/**
 * @brief Validate and analyze one radial function on the host.
 *
 * A valid grid has at least one point, finite coordinates and values, and
 * strictly increasing coordinates when it contains multiple points.
 */
bool analyze_radial_grid(const double* radial_grid,
                         const double* radial_values,
                         int mesh,
                         RadialGridInfo& grid_info);

/**
 * @brief Validate radial data and terminate with a diagnostic if it is invalid.
 */
RadialGridInfo validate_radial_grid(const double* radial_grid,
                                    const double* radial_values,
                                    int mesh,
                                    const char* caller,
                                    const char* data_name);

#if defined(__CUDACC__)
#define MODULE_RT_HOST_DEVICE __host__ __device__ __forceinline__
#else
#define MODULE_RT_HOST_DEVICE inline
#endif

namespace detail
{

MODULE_RT_HOST_DEVICE bool radial_is_finite(const double value)
{
#if defined(__CUDA_ARCH__)
    return isfinite(value);
#else
    return std::isfinite(value);
#endif
}

MODULE_RT_HOST_DEVICE int clamp_index(const int value, const int lower, const int upper)
{
    return value < lower ? lower : (value > upper ? upper : value);
}

MODULE_RT_HOST_DEVICE int find_radial_interval(const double* radial_grid,
                                               const int mesh,
                                               const RadialGridInfo& grid_info,
                                               const double radius)
{
    if (grid_info.is_uniform)
    {
        int interval = static_cast<int>((radius - grid_info.r_min) * grid_info.inv_spacing);
        interval = clamp_index(interval, 0, mesh - 2);

        // Correct the arithmetic estimate with the actual stored coordinates.
        while (interval > 0 && radius < radial_grid[interval])
        {
            --interval;
        }
        while (interval < mesh - 2 && radius > radial_grid[interval + 1])
        {
            ++interval;
        }
        return interval;
    }

    int lower = 0;
    int upper = mesh - 1;
    while (upper - lower > 1)
    {
        const int middle = lower + (upper - lower) / 2;
        if (radius < radial_grid[middle])
        {
            upper = middle;
        }
        else
        {
            lower = middle;
        }
    }
    return lower;
}

MODULE_RT_HOST_DEVICE double lagrange_interpolate(const double* radial_grid,
                                                  const double* radial_values,
                                                  const int start,
                                                  const int point_count,
                                                  const double radius)
{
    double result = 0.0;
    for (int i = 0; i < point_count; ++i)
    {
        const int ii = start + i;
        double weight = 1.0;
        for (int j = 0; j < point_count; ++j)
        {
            if (i == j)
            {
                continue;
            }
            const int jj = start + j;
            weight *= (radius - radial_grid[jj]) / (radial_grid[ii] - radial_grid[jj]);
        }
        result += weight * radial_values[ii];
    }
    return result;
}

} // namespace detail

/**
 * @brief Interpolate a radial function without extrapolating beyond its grid.
 *
 * Four or more points use a local four-point cubic Lagrange interpolant. The
 * final intervals use the final four grid points instead of returning zero.
 * Grids with one, two, or three points use constant, linear, or quadratic
 * interpolation, respectively.
 */
MODULE_RT_HOST_DEVICE double interpolate_radial(const double* radial_grid,
                                                const double* radial_values,
                                                const int mesh,
                                                const RadialGridInfo& grid_info,
                                                const double radius)
{
    if (radial_grid == nullptr || radial_values == nullptr || mesh <= 0 || !detail::radial_is_finite(radius)
        || radius < grid_info.r_min || radius > grid_info.r_max)
    {
        return 0.0;
    }

    if (radius == radial_grid[0])
    {
        return radial_values[0];
    }
    if (mesh == 1)
    {
        return 0.0;
    }
    if (radius == radial_grid[mesh - 1])
    {
        return radial_values[mesh - 1];
    }

    const int interval = detail::find_radial_interval(radial_grid, mesh, grid_info, radius);
    if (radius == radial_grid[interval])
    {
        return radial_values[interval];
    }
    if (radius == radial_grid[interval + 1])
    {
        return radial_values[interval + 1];
    }

    if (mesh < 4)
    {
        return detail::lagrange_interpolate(radial_grid, radial_values, 0, mesh, radius);
    }

    const int start = detail::clamp_index(interval - 1, 0, mesh - 4);
    return detail::lagrange_interpolate(radial_grid, radial_values, start, 4, radius);
}

#undef MODULE_RT_HOST_DEVICE

} // namespace module_rt

#endif
