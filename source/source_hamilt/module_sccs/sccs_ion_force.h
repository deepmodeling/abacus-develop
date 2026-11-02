#ifndef SCCS_ION_FORCE_H
#define SCCS_ION_FORCE_H

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
void gauss_ion_force(const std::vector<unitcell::AtomData>& atoms,
                     const std::vector<double>& reaction_potential,
                     const ModulePW::PW_Basis& basis,
                     double tpiba,
                     double spread,
                     std::vector<ModuleBase::Vector3<double>>& forces);

// One positive width per atom.
void gauss_ion_force(const std::vector<unitcell::AtomData>& atoms,
                     const std::vector<double>& reaction_potential,
                     const ModulePW::PW_Basis& basis,
                     double tpiba,
                     const std::vector<double>& spreads,
                     std::vector<ModuleBase::Vector3<double>>& forces);

// Full cavity: a singleton skips known Z == zv atoms; per-atom lists override
// that choice. Nonpositive widths always disable the atom.
void gauss_core_frc(const std::vector<unitcell::AtomData>& atoms,
                    const std::vector<double>& reaction_potential,
                    const ModulePW::PW_Basis& basis,
                    double tpiba,
                    const std::vector<double>& atom_spreads,
                    std::vector<ModuleBase::Vector3<double>>& forces);
} // namespace ModuleSccs
#endif
