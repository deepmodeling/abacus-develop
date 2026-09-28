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
}
#endif
