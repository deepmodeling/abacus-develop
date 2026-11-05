#include "sccs_ion_charge.h"

#include "source_base/constants.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/cell_tools.h"

#include <cmath>
#include <complex>

namespace ModuleSccs
{
void gauss_ion_rho(const std::vector<unitcell::AtomData>& atoms,
                   const ModulePW::PW_Basis& basis,
                   double tpiba,
                   double spread,
                   std::vector<double>& density)
{
    const std::vector<double> spreads(atoms.size(), spread);
    gauss_ion_rho(atoms, basis, tpiba, spreads, density);
}

void gauss_ion_rho(const std::vector<unitcell::AtomData>& atoms,
                   const ModulePW::PW_Basis& basis,
                   double tpiba,
                   const std::vector<double>& spreads,
                   std::vector<double>& density)
{
    std::vector<std::complex<double>> charge(basis.npw, 0.0);
    const double tpiba2 = tpiba * tpiba;
    for (std::size_t ia = 0; ia < atoms.size(); ++ia)
    {
        const unitcell::AtomData& atom = atoms[ia];
        const double spread = spreads[ia];
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            const double exponent = -0.25 * spread * spread * tpiba2 * basis.gg[ig];
            const double phase = tpiba * (basis.gcar[ig] * atom.position);
            const std::complex<double> phase_factor = ModuleBase::NEG_IMAG_UNIT * phase;
            charge[ig] += atom.valence_charge / basis.omega * std::exp(exponent) * std::exp(phase_factor);
        }
    }
    density.resize(basis.nrxx);
    basis.recip2real(charge.data(), density.data());
}

void prep_core_gauss(const std::vector<unitcell::AtomData>& atoms,
                     const std::vector<double>& atom_spreads,
                     std::vector<unitcell::AtomData>& core_atoms,
                     std::vector<double>& widths)
{
    core_atoms = atoms;
    widths.assign(atoms.size(), gaussian_ion_spread);
    // A valence charge equal to the atomic number within round-off means no core.
    const double charge_tolerance = 1e-8;
    for (std::size_t ia = 0; ia < atoms.size(); ++ia)
    {
        const std::size_t index = atom_spreads.size() == 1 ? 0 : ia;
        const double width = atom_spreads[index];
        const double charge_difference = atoms[ia].atomic_number - atoms[ia].valence_charge;
        const bool no_core = atoms[ia].atomic_number > 0 && std::abs(charge_difference) < charge_tolerance;
        const bool automatic_skip = atom_spreads.size() == 1 && no_core;
        if (width <= 0.0 || automatic_skip)
        {
            core_atoms[ia].valence_charge = 0.0;
        }
        else
        {
            widths[ia] = width;
        }
    }
}

void gauss_core_rho(const std::vector<unitcell::AtomData>& atoms,
                    const ModulePW::PW_Basis& basis,
                    double tpiba,
                    const std::vector<double>& atom_spreads,
                    std::vector<double>& density)
{
    std::vector<unitcell::AtomData> core_atoms;
    std::vector<double> widths;
    prep_core_gauss(atoms, atom_spreads, core_atoms, widths);
    gauss_ion_rho(core_atoms, basis, tpiba, widths, density);
}
} // namespace ModuleSccs
