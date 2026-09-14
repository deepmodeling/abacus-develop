#include "radial_interpolation.h"

#include "source_base/tool_quit.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace module_rt
{

bool analyze_radial_grid(const double* radial_grid,
                         const double* radial_values,
                         const int mesh,
                         RadialGridInfo& grid_info)
{
    grid_info = RadialGridInfo();
    if (radial_grid == nullptr || radial_values == nullptr || mesh <= 0)
    {
        return false;
    }

    for (int i = 0; i < mesh; ++i)
    {
        if (!std::isfinite(radial_grid[i]) || !std::isfinite(radial_values[i]))
        {
            return false;
        }
        if (i > 0 && radial_grid[i] <= radial_grid[i - 1])
        {
            return false;
        }
    }

    grid_info.r_min = radial_grid[0];
    grid_info.r_max = radial_grid[mesh - 1];
    if (mesh == 1)
    {
        return true;
    }

    const double spacing = (grid_info.r_max - grid_info.r_min) / static_cast<double>(mesh - 1);
    grid_info.inv_spacing = 1.0 / spacing;
    grid_info.is_uniform = true;

    const double scale = std::max(1.0, std::max(std::abs(grid_info.r_min), std::abs(grid_info.r_max)));
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * scale;
    for (int i = 1; i < mesh - 1; ++i)
    {
        const double expected = grid_info.r_min + static_cast<double>(i) * spacing;
        if (std::abs(radial_grid[i] - expected) > tolerance)
        {
            grid_info.is_uniform = false;
            grid_info.inv_spacing = 0.0;
            break;
        }
    }
    return true;
}

RadialGridInfo validate_radial_grid(const double* radial_grid,
                                    const double* radial_values,
                                    const int mesh,
                                    const char* caller,
                                    const char* data_name)
{
    RadialGridInfo grid_info;
    if (!analyze_radial_grid(radial_grid, radial_values, mesh, grid_info))
    {
        ModuleBase::WARNING_QUIT(caller,
                                 std::string("Invalid radial data for ") + data_name
                                     + ": pointers must be non-null, mesh must be positive, values must be finite, "
                                       "and coordinates must be strictly increasing.");
    }
    return grid_info;
}

} // namespace module_rt
