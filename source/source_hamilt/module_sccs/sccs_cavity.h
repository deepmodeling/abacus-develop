#ifndef SCCS_CAVITY_H
#define SCCS_CAVITY_H

namespace ModuleSccs
{
// Densities are in electrons/Bohr^3. The cavity uses the electronic density.
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
    double d2solute_drho2 = 0.0;
    double epsilon = 0.0;
    double depsilon_drho = 0.0;
};

// Parameters are validated by the INPUT reader. Negative Fourier ringing
// uses the bulk limit.
CavityPoint evaluate_cavity(double density, const CavityParameters& parameters);
} // namespace ModuleSccs

#endif
