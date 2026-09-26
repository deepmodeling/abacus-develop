#ifndef SCCS_PW_NONEL_H
#define SCCS_PW_NONEL_H

#include "sccs_charge.h"
#include "sccs_nonel.h"

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{

NonElectrostaticResult evaluate_pw_non_electrostatic(
    const ModulePW::PW_Basis& basis,
    double tpiba,
    double volume_element,
    const NonElectrostaticParameters& parameters,
    const std::vector<double>& solute,
    const std::vector<double>& dsolute_drho,
    const ChargeReduction& reduction);

} // namespace ModuleSccs

#endif
