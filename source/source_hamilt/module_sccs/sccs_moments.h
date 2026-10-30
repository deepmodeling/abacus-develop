#ifndef SCCS_MOMENTS_H
#define SCCS_MOMENTS_H

#include "source_estate/charge_moments.h"

#include <vector>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{
// Moments of a local grid charge about the common origin of positions,
// reduced over the PW pool, so every rank holds the same values.
elecstate::ChargeMoments pool_moments(const std::vector<double>& charge,
                                      const std::vector<ModuleBase::Vector3<double>>& positions,
                                      const ModulePW::PW_Basis& basis);
} // namespace ModuleSccs

#endif
