#ifndef SCCS_IONIC_FORCE_H
#define SCCS_IONIC_FORCE_H

#include <vector>

namespace ModuleBase
{
template <class T> class Vector3;
}
namespace ModulePW
{
class PW_Basis;
}
namespace unitcell
{
struct AtomData;
}
namespace ModuleSccs
{
// Explicit Gaussian-ion derivative of the reaction energy at fixed electron
// density. reaction_potential is in Ha; force is in Ha/Bohr, atom-major.
// Collective over the PW pool.
void gaussian_ionic_force(const std::vector<unitcell::AtomData>& atoms,
                          const std::vector<double>& reaction_potential,
                          const ModulePW::PW_Basis& basis,
                          double tpiba,
                          double spread,
                          std::vector<ModuleBase::Vector3<double>>& forces);
} // namespace ModuleSccs
#endif
