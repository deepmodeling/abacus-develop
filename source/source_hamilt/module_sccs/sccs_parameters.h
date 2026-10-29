#ifndef SCCS_PARAMETERS_H
#define SCCS_PARAMETERS_H

#include "sccs_cavity.h"

#include <string>
#include <vector>

namespace ModuleSccs
{
enum class Boundary { Periodic, Pcc0d, Pcc2d };

enum class Preset { Custom, Vacuum, WaterNeutral, WaterCation, WaterAnion };

struct SccsConfig
{
    CavityParameters cavity;
    Boundary boundary = Boundary::Periodic;
    int pcc_2d_axis = 2;
    bool core_electrons = false;
    std::vector<double> core_spreads = {0.5}; // Bohr: singleton skips Z == zv; per-atom list overrides; <= 0 disables
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
SccsConfig preset_config(Preset preset);
double dyn_cm_to_au(double value);
double gpa_to_au(double value);
} // namespace ModuleSccs

#endif
