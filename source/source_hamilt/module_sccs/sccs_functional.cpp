#include "sccs_functional.h"
#include "sccs_pw_coulomb.h"

#include "source_base/parallel_reduce.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <cmath>
#include <complex>
#include <utility>

namespace ModuleSccs
{
bool evaluate_functional(const std::vector<double>& charge,
                          const SccsResponse& response,
                          const SccsConfig& config,
                          const ModulePW::PW_Basis& basis,
                          double tpiba,
                          FunctionalResult& result,
                          std::string& error)
{
    if (!validate_pw_grid(basis, tpiba, error) || !validate_grid_values(charge, basis, error)
        || !validate_grid_values(response.polarization.potential, basis, error)
        || !validate_grid_values(response.cavity_potential, basis, error)
        || !validate_grid_values(response.solute, basis, error)
        || !validate_grid_values(response.dsolute_drho, basis, error))
    {
        return false;
    }
    double invalid = 0.0;
    if (!validate_config(config, error))
    {
        invalid = 1.0;
    }
    Parallel_Reduce::reduce_max_pool(basis.poolnproc, invalid);
    if (invalid != 0.0)
    {
        error = "SCCS functional requires valid parameters on every pool rank";
        return false;
    }
    const std::size_t size = charge.size();
    const double dv = basis.omega / basis.nxyz;
    PeriodicCoulombOperator coulomb(basis, tpiba);
    std::vector<double> vacuum;
    if (!coulomb.apply_potential(charge, vacuum, error))
    {
        return false;
    }
    FunctionalResult candidate;
    candidate.reaction_potential.resize(size);
    candidate.electron_potential.resize(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        const double reaction = response.polarization.potential[i] - vacuum[i];
        candidate.reaction_potential[i] = reaction;
        candidate.reaction_energy += 0.5 * charge[i] * reaction * dv;
        candidate.electron_potential[i] = -reaction + response.cavity_potential[i];
        candidate.volume += response.solute[i] * dv;
    }
    Parallel_Reduce::reduce_pool(candidate.reaction_energy);
    Parallel_Reduce::reduce_pool(candidate.volume);
    candidate.volume_energy = config.pressure * candidate.volume;

    // Spectral surface derivative: dE/dn = (p - gamma div(grad s/|grad s|)) s'.
    std::vector<std::complex<double>> solute_g(basis.npw);
    basis.real2recip(response.solute.data(), solute_g.data());
    std::vector<ModuleBase::Vector3<double>> unit_gradient(size);
    XC_Functional::grad_rho(solute_g.data(), unit_gradient.data(), &basis, tpiba);
    const double eta = config.surface_regularization;
    for (std::size_t i = 0; i < size; ++i)
    {
        const double norm_squared = unit_gradient[i] * unit_gradient[i] + eta * eta;
        const double norm = std::sqrt(norm_squared);
        candidate.surface += (norm - eta) * dv;
        unit_gradient[i] /= norm;
    }
    Parallel_Reduce::reduce_pool(candidate.surface);
    candidate.surface_energy = config.surface_tension * candidate.surface;
    std::vector<double> divergence(size);
    XC_Functional::grad_dot(unit_gradient.data(), divergence.data(), &basis, tpiba);
    for (std::size_t i = 0; i < size; ++i)
    {
        const double nonel = (config.pressure - config.surface_tension * divergence[i])
                             * response.dsolute_drho[i];
        candidate.electron_potential[i] += nonel;
    }
    if (!validate_grid_values(candidate.reaction_potential, basis, error)
        || !validate_grid_values(candidate.electron_potential, basis, error))
    {
        return false;
    }
    if (!std::isfinite(candidate.reaction_energy) || !std::isfinite(candidate.surface_energy)
        || !std::isfinite(candidate.volume_energy))
    {
        error = "SCCS functional energy is not finite";
        return false;
    }
    result = std::move(candidate);
    return true;
}
} // namespace ModuleSccs
