#ifndef SCCS_PARAMETERS_H
#define SCCS_PARAMETERS_H

#include "sccs_cavity.h"

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

// Failure leaves the output unchanged. Custom configurations are supplied by callers.
bool parse_preset(const std::string& name, Preset& preset, std::string& error);
bool make_sccs_config(Preset preset, SccsConfig& config, std::string& error);
bool validate_config(const SccsConfig& config, std::string& error);
double dyn_per_cm_to_hartree_per_bohr2(double value);
double gpa_to_hartree_per_bohr3(double value);
} // namespace ModuleSccs

#endif
