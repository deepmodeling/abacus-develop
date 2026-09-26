#ifndef SCCS_PW_FORCE_H
#define SCCS_PW_FORCE_H

#include "source_base/matrix.h"

#include <vector>

class UnitCell;

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{

ModuleBase::matrix smooth_ionic_force_hartree(
    const UnitCell& cell,
    const ModulePW::PW_Basis& basis,
    const ModuleBase::matrix& radial_local_potential_rydberg,
    const std::vector<double>& reaction_potential_hartree);

} // namespace ModuleSccs

#endif
