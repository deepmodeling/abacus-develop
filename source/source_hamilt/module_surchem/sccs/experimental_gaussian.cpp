#include "experimental_gaussian.h"
#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
namespace ModuleSccs {
namespace {
std::complex<double> atom_charge(const UnitCell& cell, const ModulePW::PW_Basis& basis,
                                 int type, int atom, int ig, double sigma)
{
    const double exponent = -0.25 * sigma * sigma * cell.tpiba2 * basis.gg[ig];
    const double phase = ModuleBase::TWO_PI * (basis.gcar[ig] * cell.atoms[type].tau[atom]);
    const std::complex<double> phase_factor = ModuleBase::NEG_IMAG_UNIT * phase;
    return cell.atoms[type].ncpp.zv / cell.omega * std::exp(exponent) * std::exp(phase_factor);
}
// Element label of the pseudopotential, without surrounding blanks.
bool is_hydrogen(const UnitCell& cell, int type)
{
    const std::string& label = cell.atoms[type].ncpp.psd;
    const std::size_t first = label.find_first_not_of(" \t");
    const std::size_t last = label.find_last_not_of(" \t");
    if (first == std::string::npos)
        return false;
    const std::string element = label.substr(first, last - first + 1);
    return element == "H";
}
std::vector<double> gaussian_density(const UnitCell& cell, const ModulePW::PW_Basis& basis,
                                     double sigma, bool skip_hydrogen)
{
    if (basis.gamma_only) throw std::runtime_error("Experimental Gaussian source requires full G basis");
    std::vector<std::complex<double>> charge(basis.npw, 0.0);
    for (int type = 0; type < cell.ntype; ++type)
    {
        if (skip_hydrogen && is_hydrogen(cell, type))
            continue;
        for (int atom = 0; atom < cell.atoms[type].na; ++atom)
            for (int ig = 0; ig < basis.npw; ++ig)
                charge[ig] += atom_charge(cell,basis,type,atom,ig,sigma);
    }
    std::vector<double> density(basis.nrxx);
    basis.recip2real(charge.data(),density.data());
    return density;
}
ModuleBase::matrix gaussian_force(const UnitCell& cell, const ModulePW::PW_Basis& basis,
                                  double sigma, const std::vector<double>& potential,
                                  bool skip_hydrogen)
{
    if (basis.gamma_only || potential.size() != static_cast<std::size_t>(basis.nrxx))
        throw std::runtime_error("Experimental Gaussian force input mismatch");
    std::vector<std::complex<double>> vg(basis.npw);
    basis.real2recip(potential.data(),vg.data());
    ModuleBase::matrix force(cell.nat,3);
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
                    const double derivative = std::imag(std::conj(vg[ig])
                        *atom_charge(cell,basis,type,atom,ig,sigma));
                    for (int d = 0; d < 3; ++d)
                        force(iat,d) -= cell.omega*cell.tpiba*basis.gcar[ig][d]*derivative;
                }
            }
            ++iat;
        }
    }
    return force;
}
}
std::vector<double> gaussian_ionic_density(const UnitCell& cell,
                                          const ModulePW::PW_Basis& basis, double sigma)
{
    return gaussian_density(cell, basis, sigma, false);
}
ModuleBase::matrix gaussian_ionic_force(const UnitCell& cell,
                                        const ModulePW::PW_Basis& basis, double sigma,
                                        const std::vector<double>& potential)
{
    return gaussian_force(cell, basis, sigma, potential, false);
}
std::vector<double> gaussian_core_density(const UnitCell& cell,
                                         const ModulePW::PW_Basis& basis, double sigma)
{
    return gaussian_density(cell, basis, sigma, true);
}
ModuleBase::matrix gaussian_core_force(const UnitCell& cell,
                                       const ModulePW::PW_Basis& basis, double sigma,
                                       const std::vector<double>& potential)
{
    return gaussian_force(cell, basis, sigma, potential, true);
}
}
