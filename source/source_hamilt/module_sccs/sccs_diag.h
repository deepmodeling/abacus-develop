#ifndef SCCS_DIAG_H
#define SCCS_DIAG_H

#include <vector>

namespace ModulePW { class PW_Basis; }
namespace ModuleSccs
{
struct SccsResponse;

// Environ dielectric_of_potential; diagnostic density, not a new solve source.
void cont_pol_charge(const std::vector<double>& charge,
                     const SccsResponse& response,
                     const ModulePW::PW_Basis& basis,
                     double tpiba,
                     std::vector<double>& density);
} // namespace ModuleSccs
#endif
