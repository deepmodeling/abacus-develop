#ifndef EXPERIMENTAL_GAUSSIAN_H
#define EXPERIMENTAL_GAUSSIAN_H
#include <vector>
class UnitCell;
namespace ModulePW { class PW_Basis; }
namespace ModuleBase { class matrix; }
namespace ModuleSccs {
// Spread (bohr) of the Gaussian ionic source, rho ~ exp(-r^2/spread^2), shared by
// the SCCS reaction energy and its ionic force; ENVIRON's default atomicspread.
const double gaussian_ion_spread = 0.5;
std::vector<double> gaussian_ionic_density(const UnitCell&, const ModulePW::PW_Basis&, double);
ModuleBase::matrix gaussian_ionic_force(const UnitCell&, const ModulePW::PW_Basis&, double,
                                       const std::vector<double>&);
// ENVIRON solvent_mode 'full' core electrons: a Gaussian of the valence charge
// zv and the given spread on every atom except hydrogen (ENVIRON corespread
// 1e-10 for H). It fills the pseudo-valence hole at the nuclei in the cavity
// density. The force is minus the derivative of the integral of the cavity
// potential times this density.
std::vector<double> gaussian_core_density(const UnitCell&, const ModulePW::PW_Basis&, double);
ModuleBase::matrix gaussian_core_force(const UnitCell&, const ModulePW::PW_Basis&, double,
                                      const std::vector<double>&);
}
#endif
