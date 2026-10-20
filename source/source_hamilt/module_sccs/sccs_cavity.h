#ifndef SCCS_CAVITY_H
#define SCCS_CAVITY_H

namespace ModuleSccs
{
// Environ solvent-aware interface (Andreussi et al., JCTC 15, 1996 (2019)):
// cavity voids too small for a solvent molecule are filled with solute.
// Lengths are in Bohr; validated by the INPUT reader.
struct SolventAwareParameters
{
    double solvent_radius = 0.0; // zero disables the filling
    double radial_scale = 2.0; // probe radius over solvent_radius
    double radial_spread = 0.5; // erfc spread of the probe sphere
    double filling_threshold = 0.825; // probe solute fraction that fills a point
    double filling_spread = 0.02; // erfc width of the filling step
};

// Densities are in electrons/Bohr^3. The cavity uses the electronic density.
struct CavityParameters
{
    double density_min = 0.0;
    double density_max = 0.0;
    double epsilon_bulk = 0.0;
    // Both positive enable the PCC switching-function derivative filter.
    double lowpass_p1 = -1.0;
    double lowpass_p2 = -1.0;
    SolventAwareParameters solvent_aware;
};

struct CavityPoint
{
    double solute = 0.0;
    double dsolute_drho = 0.0;
    double d2solute_drho2 = 0.0;
    double epsilon = 0.0;
    double depsilon_drho = 0.0;
};

bool uses_lowpass(const CavityParameters& parameters);

bool uses_sa(const SolventAwareParameters& parameters);

// Parameters are validated by the INPUT reader. Negative Fourier ringing
// uses the bulk limit.
CavityPoint evaluate_cavity(double density, const CavityParameters& parameters);
} // namespace ModuleSccs

#endif
