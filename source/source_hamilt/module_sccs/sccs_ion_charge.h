#ifndef SCCS_ION_CHARGE_H
#define SCCS_ION_CHARGE_H

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
// rho ~ exp(-r^2/spread^2), in Bohr; the Environ default ionic spread.
const double gaussian_ion_spread = 0.5;

// AtomData is extracted by source_cell; positions are Cartesian Bohr.
// Collective over the PW pool.
void gauss_ion_rho(const std::vector<unitcell::AtomData>& atoms,
                   const ModulePW::PW_Basis& basis,
                   double tpiba,
                   double spread,
                   std::vector<double>& density);

// One positive width per atom.
void gauss_ion_rho(const std::vector<unitcell::AtomData>& atoms,
                   const ModulePW::PW_Basis& basis,
                   double tpiba,
                   const std::vector<double>& spreads,
                   std::vector<double>& density);

// Resolve full-cavity widths without changing the physical ionic charge.
// atom_spreads holds one value or one per atom (checked when SCCS is set up).
void prep_core_gauss(const std::vector<unitcell::AtomData>& atoms,
                     const std::vector<double>& atom_spreads,
                     std::vector<unitcell::AtomData>& core_atoms,
                     std::vector<double>& widths);

// Full cavity, as Environ 3.1.1 solvent_mode 'full' (core_electrons in
// ions.f90, used by boundary_electronic.f90): one Gaussian per atom with its
// valence charge and width corespread, added to the electron density that
// defines the cavity.
// A singleton skips known Z == zv atoms; per-atom lists override
// that choice. Nonpositive widths always disable the atom.
void gauss_core_rho(const std::vector<unitcell::AtomData>& atoms,
                    const ModulePW::PW_Basis& basis,
                    double tpiba,
                    const std::vector<double>& atom_spreads,
                    std::vector<double>& density);
} // namespace ModuleSccs

#endif
