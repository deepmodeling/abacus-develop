#ifndef SCCS_IONIC_CHARGE_H
#define SCCS_IONIC_CHARGE_H

#include <string>
#include <vector>

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
// rho ~ exp(-r^2/spread^2), in Bohr; preserves the original SCCS ionic spread.
const double gaussian_ion_spread = 0.5;

// AtomData is extracted by source_cell; positions are Cartesian Bohr.
// Collective over the PW pool. Failure leaves density unchanged.
bool gaussian_ionic_density(const std::vector<unitcell::AtomData>& atoms,
                            const ModulePW::PW_Basis& basis,
                            double tpiba,
                            double spread,
                            std::vector<double>& density,
                            std::string& error);
} // namespace ModuleSccs

#endif
