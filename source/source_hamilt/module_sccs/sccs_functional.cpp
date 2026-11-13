#include "sccs_functional.h"
#include "sccs_parameters.h"
#include "sccs_response.h"
#include "sccs_sa.h"
#include "sccs_thread_sum.h"
#include "sccs_pw_coulomb.h"

#include "source_base/parallel_reduce.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <cmath>
#include <complex>
#include <utility>

namespace ModuleSccs
{
namespace
{
// Spectral surface int(|grad s|_r - r), r = surface_regularization, and
// dE/ds = pressure - gamma div(grad s/|grad s|_r).
void spec_surface(const SccsResponse& response,
                  const SccsConfig& config,
                  const ModulePW::PW_Basis& basis,
                  double tpiba,
                  double& surface,
                  std::vector<double>& boundary_potential)
{
    const std::size_t size = response.solute.size();
    const double dv = basis.omega / basis.nxyz;
    std::vector<std::complex<double>> solute_g(basis.npw);
    basis.real2recip(response.solute.data(), solute_g.data());
    std::vector<ModuleBase::Vector3<double>> unit_gradient(size);
    XC_Functional::grad_rho(solute_g.data(), unit_gradient.data(), &basis, tpiba);
    const double eta = config.surface_regularization;
    const auto add_point = [&](std::size_t i, std::array<double, 1>& sum) {
        const double norm_squared = unit_gradient[i] * unit_gradient[i] + eta * eta;
        const double norm = std::sqrt(norm_squared);
        sum[0] += (norm - eta) * dv;
        unit_gradient[i] /= norm;
    };
    surface = thread_sums<1>(size, add_point)[0];
    Parallel_Reduce::reduce_pool(surface);
    std::vector<double> divergence(size);
    XC_Functional::grad_dot(unit_gradient.data(), divergence.data(), &basis, tpiba);
    boundary_potential.resize(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        boundary_potential[i] = config.pressure - config.surface_tension * divergence[i];
    }
}

// Surface of the filled boundary from the chain gradient of s_sa (Environ
// 'chain'); the spectral gradient of s_sa rings where the filling switches
// within a grid spacing, so every boundary condition uses it. density_potential
// is the exact derivative of the discrete surface and volume energies (see
// sa_nonelec), bounded where the gradient vanishes.
void filled_nonelec(const SccsResponse& response,
                    const SccsConfig& config,
                    const std::vector<double>& probe_kernel,
                    const ModulePW::PW_Basis& basis,
                    double tpiba,
                    double& surface,
                    std::vector<double>& density_potential)
{
    const double eta = config.surface_regularization;
    SolventAwareNonelectrostatic filled = sa_nonelec(
        response.solvent_aware, response.dsolute_drho, probe_kernel, basis, tpiba, eta, config.surface_tension,
        config.pressure);
    const double dv = basis.omega / basis.nxyz;
    const std::vector<ModuleBase::Vector3<double>>& gradient = filled.gradient;
    const auto add_point = [&](std::size_t i, std::array<double, 1>& sum) {
        const double norm_squared = gradient[i].norm2() + eta * eta;
        const double norm = std::sqrt(norm_squared);
        sum[0] += (norm - eta) * dv;
    };
    surface = thread_sums<1>(gradient.size(), add_point)[0];
    Parallel_Reduce::reduce_pool(surface);
    density_potential.swap(filled.density_potential);
}
} // namespace

void eval_functional(const std::vector<double>& charge,
                      const SccsResponse& response,
                      const SccsConfig& config,
                      const ModulePW::PW_Basis& basis,
                      double tpiba,
                      FunctionalResult& result)
{
    PeriodicCoulombOperator coulomb(basis, tpiba);
    const std::vector<double> no_probe;
    eval_functional(charge, response, config, no_probe, basis, tpiba, coulomb, result);
}

void eval_functional(const std::vector<double>& charge,
                      const SccsResponse& response,
                      const SccsConfig& config,
                      const std::vector<double>& probe_kernel,
                      const ModulePW::PW_Basis& basis,
                      double tpiba,
                      CoulombOperator& coulomb,
                      FunctionalResult& result)
{
    const std::size_t size = charge.size();
    const double dv = basis.omega / basis.nxyz;
    std::vector<double> vacuum;
    coulomb.apply_potential(charge, vacuum);
    FunctionalResult candidate;
    candidate.reaction_potential.resize(size);
    candidate.electron_potential.resize(size);
    // The cavity-density derivative is needed only for the forces on the core
    // Gaussians of the full cavity.
    const bool core_electrons = config.core_electrons;
    if (core_electrons)
    {
        candidate.cavity_potential = response.cavity_potential;
    }
    const auto add_point = [&](std::size_t i, std::array<double, 2>& sums) {
        const double reaction = response.polarization.potential[i] - vacuum[i];
        candidate.reaction_potential[i] = reaction;
        candidate.electron_potential[i] = -reaction + response.cavity_potential[i];
        sums[0] += 0.5 * charge[i] * reaction * dv;
        sums[1] += response.solute[i] * dv;
    };
    const std::array<double, 2> sums = thread_sums<2>(size, add_point);
    candidate.reaction_energy = sums[0];
    candidate.volume = sums[1];
    Parallel_Reduce::reduce_pool(candidate.reaction_energy);
    Parallel_Reduce::reduce_pool(candidate.volume);
    candidate.volume_energy = config.pressure * candidate.volume;

    std::vector<double> nonel;
    if (response.solvent_aware.enabled())
    {
        filled_nonelec(response, config, probe_kernel, basis, tpiba, candidate.surface, nonel);
    }
    else
    {
        std::vector<double> boundary_potential;
        spec_surface(response, config, basis, tpiba, candidate.surface, boundary_potential);
        boundary_to_rho(response, probe_kernel, basis, boundary_potential, nonel);
    }
    candidate.surface_energy = config.surface_tension * candidate.surface;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        candidate.electron_potential[i] += nonel[i];
        if (core_electrons)
        {
            candidate.cavity_potential[i] += nonel[i];
        }
    }
    result = std::move(candidate);
}
} // namespace ModuleSccs
