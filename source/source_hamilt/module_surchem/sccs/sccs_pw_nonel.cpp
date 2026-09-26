#include "sccs_pw_nonel.h"

#include "sccs_pw_coulomb.h"

#include "source_base/constants.h"
#include "source_basis/module_pw/pw_basis.h"

#include <cmath>
#include <complex>
#include <stdexcept>

namespace ModuleSccs
{

NonElectrostaticResult evaluate_pw_non_electrostatic(
    const ModulePW::PW_Basis& basis,
    const double tpiba,
    const double volume_element,
    const NonElectrostaticParameters& parameters,
    const std::vector<double>& solute,
    const std::vector<double>& dsolute_drho,
    const ChargeReduction& reduction)
{
    if (solute.size() != static_cast<std::size_t>(basis.nrxx)
        || dsolute_drho.size() != solute.size() || solute.empty())
    {
        throw std::invalid_argument("SCCS PW non-electrostatic arrays must match the local real-space grid");
    }
    if (!std::isfinite(tpiba) || tpiba <= 0.0 || !std::isfinite(volume_element)
        || volume_element <= 0.0 || !std::isfinite(parameters.surface_tension)
        || !std::isfinite(parameters.pressure)
        || !std::isfinite(parameters.surface_regularization)
        || parameters.surface_regularization <= 0.0)
    {
        throw std::invalid_argument("SCCS PW non-electrostatic parameters must be finite and valid");
    }

    const std::vector<ModuleBase::Vector3<double>> gradient
        = periodic_gradient(solute, basis, tpiba);
    std::vector<ModuleBase::Vector3<double>> unit_gradient(solute.size());
    NonElectrostaticResult result;
    result.density_potential.resize(solute.size());
    for (std::size_t index = 0; index < solute.size(); ++index)
    {
        if (!std::isfinite(solute[index]) || !std::isfinite(dsolute_drho[index]))
        {
            throw std::domain_error("SCCS PW non-electrostatic inputs must be finite");
        }
        const double norm = std::sqrt(gradient[index].x * gradient[index].x
                                      + gradient[index].y * gradient[index].y
                                      + gradient[index].z * gradient[index].z
                                      + parameters.surface_regularization
                                            * parameters.surface_regularization);
        unit_gradient[index].x = gradient[index].x / norm;
        unit_gradient[index].y = gradient[index].y / norm;
        unit_gradient[index].z = gradient[index].z / norm;
        result.surface += (norm - parameters.surface_regularization) * volume_element;
        result.volume += solute[index] * volume_element;
    }
    reduction.reduce_sum(result.surface);
    reduction.reduce_sum(result.volume);

    std::vector<std::complex<double>> component_g(basis.npw);
    std::vector<std::complex<double>> divergence_g(basis.npw, std::complex<double>());
    std::vector<double> component_r(solute.size());
    for (int direction = 0; direction < 3; ++direction)
    {
        for (std::size_t index = 0; index < solute.size(); ++index)
        {
            component_r[index] = unit_gradient[index][direction];
        }
        basis.real2recip(component_r.data(), component_g.data());
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            divergence_g[ig] += ModuleBase::IMAG_UNIT * tpiba
                                * basis.gcar[ig][direction] * component_g[ig];
        }
    }
    std::vector<double> divergence(solute.size());
    basis.recip2real(divergence_g.data(), divergence.data());

    result.surface_energy = parameters.surface_tension * result.surface;
    result.volume_energy = parameters.pressure * result.volume;
    for (std::size_t index = 0; index < solute.size(); ++index)
    {
        result.density_potential[index]
            = (parameters.pressure - parameters.surface_tension * divergence[index])
              * dsolute_drho[index];
    }
    return result;
}

} // namespace ModuleSccs
