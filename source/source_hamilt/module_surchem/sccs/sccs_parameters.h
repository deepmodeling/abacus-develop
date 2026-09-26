#ifndef SCCS_PARAMETERS_H
#define SCCS_PARAMETERS_H

#include "sccs_cavity.h"

#include <string>

namespace ModuleSccs
{

enum class Preset
{
    Custom,
    Vacuum,
    WaterNeutral,
    WaterCation,
    WaterAnion
};

enum class Boundary
{
    Periodic,
    Pcc0d,
    Pcc2d
};

struct SccsConfig
{
    CavityParameters cavity;
    double surface_tension = 0.0;
    double pressure = 0.0;
    double surface_regularization = 0.0;
    Boundary boundary = Boundary::Periodic;
    int max_iterations = 0;
    std::string mixing_method = "linear";
    int mixing_history = 8;
    double mixing = 0.0;
    bool adaptive_mixing = false;
    double mixing_min = 0.1;
    double mixing_max = 0.8;
    double tolerance_rms = 0.0;
    double tolerance_max = 0.0;
};

Preset parse_preset(const std::string& value);

Boundary parse_boundary(const std::string& value);

SccsConfig water_preset(Preset preset);

SccsConfig vacuum_preset();

void validate_config(const SccsConfig& config);

double dyn_per_cm_to_hartree_per_bohr2(double value);

double gpa_to_hartree_per_bohr3(double value);

} // namespace ModuleSccs

#endif
