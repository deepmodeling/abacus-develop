#ifndef SCCS_CAVITY_H
#define SCCS_CAVITY_H

namespace ModuleSccs
{

struct CavityParameters
{
    double density_min = 0.0;
    double density_max = 0.0;
    double epsilon_bulk = 0.0;
};

struct CavityPoint
{
    double solute = 0.0;
    double dsolute_drho = 0.0;
    double epsilon = 0.0;
    double depsilon_drho = 0.0;
};

void validate_cavity_parameters(const CavityParameters& parameters);

CavityPoint evaluate_cavity(double density, const CavityParameters& parameters);

} // namespace ModuleSccs

#endif
