#ifndef SCCS_NONEL_H
#define SCCS_NONEL_H

#include <vector>

namespace ModuleSccs
{

struct UniformGrid
{
    int nx = 0;
    int ny = 0;
    int nz = 0;
    double spacing_x = 0.0;
    double spacing_y = 0.0;
    double spacing_z = 0.0;
};

struct NonElectrostaticParameters
{
    double surface_tension = 0.0;
    double pressure = 0.0;
    double surface_regularization = 0.0;
};

struct NonElectrostaticResult
{
    double surface = 0.0;
    double volume = 0.0;
    double surface_energy = 0.0;
    double volume_energy = 0.0;
    std::vector<double> density_potential;
};

NonElectrostaticResult evaluate_non_electrostatic(const UniformGrid& grid,
                                                  const NonElectrostaticParameters& parameters,
                                                  const std::vector<double>& solute,
                                                  const std::vector<double>& dsolute_drho);

} // namespace ModuleSccs

#endif
