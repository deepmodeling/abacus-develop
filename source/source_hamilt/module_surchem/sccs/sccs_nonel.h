#ifndef SCCS_NONEL_H
#define SCCS_NONEL_H

#include <vector>

namespace ModuleSccs
{

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

} // namespace ModuleSccs

#endif
