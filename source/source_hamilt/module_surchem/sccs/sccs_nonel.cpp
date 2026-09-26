#include "sccs_nonel.h"

#include <cmath>
#include <stdexcept>

namespace ModuleSccs
{
namespace
{

std::size_t grid_size(const UniformGrid& grid)
{
    if (grid.nx <= 0 || grid.ny <= 0 || grid.nz <= 0)
    {
        throw std::invalid_argument("SCCS non-electrostatic grid dimensions must be positive");
    }
    if (!std::isfinite(grid.spacing_x) || !std::isfinite(grid.spacing_y)
        || !std::isfinite(grid.spacing_z) || grid.spacing_x <= 0.0
        || grid.spacing_y <= 0.0 || grid.spacing_z <= 0.0)
    {
        throw std::invalid_argument("SCCS non-electrostatic grid spacings must be positive and finite");
    }
    return static_cast<std::size_t>(grid.nx) * static_cast<std::size_t>(grid.ny)
           * static_cast<std::size_t>(grid.nz);
}

std::size_t index_3d(const int ix, const int iy, const int iz, const UniformGrid& grid)
{
    return (static_cast<std::size_t>(ix) * static_cast<std::size_t>(grid.ny)
            + static_cast<std::size_t>(iy))
               * static_cast<std::size_t>(grid.nz)
           + static_cast<std::size_t>(iz);
}

void validate_parameters(const NonElectrostaticParameters& parameters)
{
    if (!std::isfinite(parameters.surface_tension) || !std::isfinite(parameters.pressure)
        || !std::isfinite(parameters.surface_regularization))
    {
        throw std::invalid_argument("SCCS non-electrostatic parameters must be finite");
    }
    if (parameters.surface_regularization <= 0.0)
    {
        throw std::invalid_argument("SCCS surface regularization must be positive");
    }
}

} // namespace

NonElectrostaticResult evaluate_non_electrostatic(const UniformGrid& grid,
                                                  const NonElectrostaticParameters& parameters,
                                                  const std::vector<double>& solute,
                                                  const std::vector<double>& dsolute_drho)
{
    const std::size_t size = grid_size(grid);
    validate_parameters(parameters);
    if (solute.size() != size || dsolute_drho.size() != size)
    {
        throw std::invalid_argument("SCCS non-electrostatic arrays do not match the grid size");
    }

    const double volume_element = grid.spacing_x * grid.spacing_y * grid.spacing_z;
    std::vector<double> derivative_solute(size, parameters.pressure * volume_element);
    NonElectrostaticResult result;
    result.density_potential.assign(size, 0.0);

    for (int ix = 0; ix < grid.nx; ++ix)
    {
        const int ix_next = (ix + 1) % grid.nx;
        for (int iy = 0; iy < grid.ny; ++iy)
        {
            const int iy_next = (iy + 1) % grid.ny;
            for (int iz = 0; iz < grid.nz; ++iz)
            {
                const int iz_next = (iz + 1) % grid.nz;
                const std::size_t center = index_3d(ix, iy, iz, grid);
                const std::size_t next_x = index_3d(ix_next, iy, iz, grid);
                const std::size_t next_y = index_3d(ix, iy_next, iz, grid);
                const std::size_t next_z = index_3d(ix, iy, iz_next, grid);
                const double gradient_x = (solute[next_x] - solute[center]) / grid.spacing_x;
                const double gradient_y = (solute[next_y] - solute[center]) / grid.spacing_y;
                const double gradient_z = (solute[next_z] - solute[center]) / grid.spacing_z;
                const double norm = std::sqrt(gradient_x * gradient_x + gradient_y * gradient_y
                                              + gradient_z * gradient_z
                                              + parameters.surface_regularization
                                                    * parameters.surface_regularization);
                const double surface_density = norm - parameters.surface_regularization;
                const double coefficient = parameters.surface_tension * volume_element / norm;

                result.surface += surface_density * volume_element;
                result.volume += solute[center] * volume_element;

                const double derivative_x = coefficient * gradient_x / grid.spacing_x;
                const double derivative_y = coefficient * gradient_y / grid.spacing_y;
                const double derivative_z = coefficient * gradient_z / grid.spacing_z;
                derivative_solute[center] -= derivative_x + derivative_y + derivative_z;
                derivative_solute[next_x] += derivative_x;
                derivative_solute[next_y] += derivative_y;
                derivative_solute[next_z] += derivative_z;
            }
        }
    }

    result.surface_energy = parameters.surface_tension * result.surface;
    result.volume_energy = parameters.pressure * result.volume;
    for (std::size_t index = 0; index < size; ++index)
    {
        result.density_potential[index]
            = derivative_solute[index] * dsolute_drho[index] / volume_element;
    }
    return result;
}

} // namespace ModuleSccs
