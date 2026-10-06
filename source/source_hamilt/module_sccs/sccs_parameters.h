#ifndef SCCS_PARAMETERS_H
#define SCCS_PARAMETERS_H

#include "sccs_cavity.h"

#include <string>

namespace ModuleSccs
{
enum class Preset { Custom, Vacuum, WaterNeutral, WaterCation, WaterAnion };

struct SccsConfig
{
    CavityParameters cavity;
    double surface_tension = 0.0; // Hartree/Bohr^2
    double pressure = 0.0; // Hartree/Bohr^3
    double surface_regularization = 0.0;
};

struct PolarizationSolverParameters
{
    int max_iterations = 200;
    double tolerance_rms = 1e-10;
    double tolerance_max = 1e-8;
};

// Exact, lower-case names as accepted by the INPUT reader; an unknown name
// stops the run.
Preset parse_preset(const std::string& name);
// Published presets only; custom configurations are supplied by callers.
SccsConfig make_sccs_config(Preset preset);
double dyn_per_cm_to_hartree_per_bohr2(double value);
double gpa_to_hartree_per_bohr3(double value);
} // namespace ModuleSccs

#endif
