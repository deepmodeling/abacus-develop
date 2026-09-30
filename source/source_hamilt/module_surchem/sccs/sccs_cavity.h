#ifndef SCCS_CAVITY_H
#define SCCS_CAVITY_H

namespace ModuleSccs
{

struct CavityParameters
{
    double density_min = 0.0;
    double density_max = 0.0;
    double epsilon_bulk = 0.0;
    // Environ deriv_lowpass_p1/p2: both positive filter the switching-function
    // derivatives by 0.5 erfc(p1 G^2/Gcut^2 - p2); both non-positive disable it.
    double lowpass_p1 = -1.0;
    double lowpass_p2 = -1.0;
};

struct CavityPoint
{
    double solute = 0.0;
    double dsolute_drho = 0.0;
    double epsilon = 0.0;
    double depsilon_drho = 0.0;
};

void validate_cavity_parameters(const CavityParameters& parameters);

bool uses_switching_lowpass(const CavityParameters& parameters);

// Per grid point; the caller validates parameters once with
// validate_cavity_parameters.
CavityPoint evaluate_cavity(double density, const CavityParameters& parameters);

} // namespace ModuleSccs

#endif
