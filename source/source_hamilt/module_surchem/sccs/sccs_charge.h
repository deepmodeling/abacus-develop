#ifndef SCCS_CHARGE_H
#define SCCS_CHARGE_H

#include "../common/charge_reduction.h"

#include <vector>

namespace ModuleSccs
{

struct ChargeDensity
{
    std::vector<double> electron;
    std::vector<double> solute;
    double electron_count = 0.0;
    double ionic_charge = 0.0;
    double net_charge = 0.0;
};

std::vector<double> sum_electron_density(const std::vector<std::vector<double>>& spin_density,
                                         int charge_channels);

// Solute charge q = rho_ion - n. The ionic source must integrate to the
// valence charge; the electron count is reported, not enforced, because the
// grid density of a converging SCF need not integrate to it exactly.
ChargeDensity assemble_charge_density(const std::vector<double>& electron_density,
                                      const std::vector<double>& ionic_density,
                                      double volume_element,
                                      double expected_ionic_charge,
                                      double normalization_tolerance,
                                      const ModuleSurchem::ChargeReduction& reduction);

} // namespace ModuleSccs

#endif
