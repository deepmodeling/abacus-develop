#include "sccs_gaussian_ion.h"

#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>

namespace ModuleSccs
{
namespace
{

// Fourier coefficient at G of the Gaussian of charge zv on one atom.
std::complex<double> atom_charge(const UnitCell& cell,
                                 const ModulePW::PW_Basis& basis,
                                 const int type,
                                 const int atom,
                                 const int ig,
                                 const double spread)
{
    const double exponent = -0.25 * spread * spread * cell.tpiba2 * basis.gg[ig];
    const double phase = ModuleBase::TWO_PI * (basis.gcar[ig] * cell.atoms[type].tau[atom]);
    const std::complex<double> phase_factor = ModuleBase::NEG_IMAG_UNIT * phase;
    return cell.atoms[type].ncpp.zv / cell.omega * std::exp(exponent) * std::exp(phase_factor);
}

// Element label of the pseudopotential, without surrounding blanks.
bool is_hydrogen(const UnitCell& cell, const int type)
{
    const std::string& label = cell.atoms[type].ncpp.psd;
    const std::size_t first = label.find_first_not_of(" \t");
    const std::size_t last = label.find_last_not_of(" \t");
    if (first == std::string::npos)
    {
        return false;
    }
    const std::string element = label.substr(first, last - first + 1);
    return element == "H";
}

std::vector<double> gaussian_density(const UnitCell& cell,
                                     const ModulePW::PW_Basis& basis,
                                     const double spread,
                                     const bool skip_hydrogen)
{
    if (basis.gamma_only)
    {
        throw std::runtime_error("SCCS Gaussian ionic source requires the full G basis");
    }
    std::vector<std::complex<double>> charge(basis.npw, 0.0);
    for (int type = 0; type < cell.ntype; ++type)
    {
        if (skip_hydrogen && is_hydrogen(cell, type))
        {
            continue;
        }
        for (int atom = 0; atom < cell.atoms[type].na; ++atom)
        {
            for (int ig = 0; ig < basis.npw; ++ig)
            {
                charge[ig] += atom_charge(cell, basis, type, atom, ig, spread);
            }
        }
    }
    std::vector<double> density(basis.nrxx);
    basis.recip2real(charge.data(), density.data());
    return density;
}

ModuleBase::matrix gaussian_force(const UnitCell& cell,
                                  const ModulePW::PW_Basis& basis,
                                  const double spread,
                                  const std::vector<double>& potential,
                                  const bool skip_hydrogen)
{
    if (basis.gamma_only || potential.size() != static_cast<std::size_t>(basis.nrxx))
    {
        throw std::runtime_error("SCCS Gaussian ionic force input mismatch");
    }
    std::vector<std::complex<double>> potential_g(basis.npw);
    basis.real2recip(potential.data(), potential_g.data());
    ModuleBase::matrix force(cell.nat, 3);
    int iat = 0;
    for (int type = 0; type < cell.ntype; ++type)
    {
        const bool skip = skip_hydrogen && is_hydrogen(cell, type);
        for (int atom = 0; atom < cell.atoms[type].na; ++atom)
        {
            if (!skip)
            {
                for (int ig = 0; ig < basis.npw; ++ig)
                {
                    const std::complex<double> charge = atom_charge(cell, basis, type, atom, ig, spread);
                    const double derivative = std::imag(std::conj(potential_g[ig]) * charge);
                    for (int d = 0; d < 3; ++d)
                    {
                        force(iat, d) -= cell.omega * cell.tpiba * basis.gcar[ig][d] * derivative;
                    }
                }
            }
            ++iat;
        }
    }
    return force;
}

} // namespace

std::vector<double> gaussian_ionic_density(const UnitCell& cell,
                                           const ModulePW::PW_Basis& basis,
                                           const double spread)
{
    return gaussian_density(cell, basis, spread, false);
}

ModuleBase::matrix gaussian_ionic_force(const UnitCell& cell,
                                        const ModulePW::PW_Basis& basis,
                                        const double spread,
                                        const std::vector<double>& potential)
{
    return gaussian_force(cell, basis, spread, potential, false);
}

std::vector<double> gaussian_core_density(const UnitCell& cell,
                                          const ModulePW::PW_Basis& basis,
                                          const double spread)
{
    return gaussian_density(cell, basis, spread, true);
}

ModuleBase::matrix gaussian_core_force(const UnitCell& cell,
                                       const ModulePW::PW_Basis& basis,
                                       const double spread,
                                       const std::vector<double>& potential)
{
    return gaussian_force(cell, basis, spread, potential, true);
}

} // namespace ModuleSccs
