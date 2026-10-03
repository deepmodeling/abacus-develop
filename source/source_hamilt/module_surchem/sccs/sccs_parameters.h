#ifndef SCCS_PARAMETERS_H
#define SCCS_PARAMETERS_H

#include "../pcc/pcc_boundary.h"
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

struct SccsConfig
{
    CavityParameters cavity;
    double surface_tension = 0.0;
    double pressure = 0.0;
    double surface_regularization = 0.0;
    ModulePcc::Boundary boundary = ModulePcc::Boundary::Periodic;
    int max_iterations = 0;
    double tolerance_rms = 0.0;
    double tolerance_max = 0.0;
    bool check_fixed_point = false;
    // ENVIRON solvent_mode 'full': the cavity density adds a Gaussian of the
    // valence charge with core_spread (bohr) on every non-hydrogen atom.
    bool core_electrons = false;
    double core_spread = 0.5;
};

Preset parse_preset(const std::string& value);

SccsConfig water_preset(Preset preset);

SccsConfig vacuum_preset();

void validate_config(const SccsConfig& config);

double dyn_per_cm_to_hartree_per_bohr2(double value);

double gpa_to_hartree_per_bohr3(double value);

} // namespace ModuleSccs

#endif
