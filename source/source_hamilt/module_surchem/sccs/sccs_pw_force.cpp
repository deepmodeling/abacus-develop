#include "sccs_pw_force.h"

#include "source_base/constants.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace ModuleSccs
{

ModuleBase::matrix smooth_ionic_force_hartree(
    const UnitCell& cell,
    const ModulePW::PW_Basis& basis,
    const ModuleBase::matrix& radial_local_potential_rydberg,
    const std::vector<double>& reaction_potential_hartree)
{
    if (cell.nat <= 0 || cell.ntype <= 0 || basis.npw <= 0 || basis.nrxx <= 0
        || basis.gg == nullptr || basis.gcar == nullptr || basis.ig2igg == nullptr
        || reaction_potential_hartree.size() != static_cast<std::size_t>(basis.nrxx)
        || radial_local_potential_rydberg.nr != cell.ntype)
    {
        throw std::invalid_argument("SCCS smooth ionic force inputs are inconsistent");
    }
    if (!std::isfinite(cell.omega) || cell.omega <= 0.0
        || !std::isfinite(cell.tpiba) || cell.tpiba <= 0.0)
    {
        throw std::invalid_argument("SCCS smooth ionic force requires a valid cell");
    }

    std::vector<std::complex<double>> reaction_potential_g(basis.npw);
    basis.real2recip(reaction_potential_hartree.data(), reaction_potential_g.data());
    ModuleBase::matrix force(cell.nat, 3);
    const double tpiba2 = cell.tpiba * cell.tpiba;
    int atom_index = 0;
    for (int atom_type = 0; atom_type < cell.ntype; ++atom_type)
    {
        for (int atom = 0; atom < cell.atoms[atom_type].na; ++atom)
        {
            for (int ig = 0; ig < basis.npw; ++ig)
            {
                if (ig == basis.ig_gge0 || basis.gg[ig] == 0.0)
                {
                    continue;
                }
                const int radial_index = basis.ig2igg[ig];
                if (radial_index < 0 || radial_index >= radial_local_potential_rydberg.nc)
                {
                    throw std::out_of_range("SCCS local-potential reciprocal index is invalid");
                }
                const std::complex<double> phase
                    = std::exp(ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI
                               * (basis.gcar[ig] * cell.atoms[atom_type].tau[atom]));
                const double coulomb_rydberg
                    = ModuleBase::e2 * ModuleBase::FOUR_PI / (tpiba2 * basis.gg[ig]);
                const std::complex<double> ionic_charge_g
                    = -radial_local_potential_rydberg(atom_type, radial_index) * phase
                      / coulomb_rydberg;
                const double spectral_derivative
                    = std::imag(std::conj(reaction_potential_g[ig]) * ionic_charge_g);
                force(atom_index, 0)
                    -= cell.omega * cell.tpiba * basis.gcar[ig][0] * spectral_derivative;
                force(atom_index, 1)
                    -= cell.omega * cell.tpiba * basis.gcar[ig][1] * spectral_derivative;
                force(atom_index, 2)
                    -= cell.omega * cell.tpiba * basis.gcar[ig][2] * spectral_derivative;
            }
            ++atom_index;
        }
    }
    return force;
}

} // namespace ModuleSccs
