#include "sccs_charge.h"

#include <cmath>
#include <stdexcept>

namespace ModuleSccs
{

void SerialChargeReduction::reduce_sum(double& value) const
{
    if (!std::isfinite(value))
    {
        throw std::domain_error("SCCS charge integral must be finite before reduction");
    }
}

void SerialChargeReduction::reduce_sum(double* values, const int count) const
{
    if (count < 0 || (count > 0 && values == nullptr))
    {
        throw std::invalid_argument("SCCS charge-array reduction requires valid storage and size");
    }
    for (int index = 0; index < count; ++index)
    {
        if (!std::isfinite(values[index]))
        {
            throw std::domain_error("SCCS charge array must be finite before reduction");
        }
    }
}

std::vector<double> sum_electron_density(const std::vector<std::vector<double>>& spin_density,
                                         const int charge_channels)
{
    if (spin_density.empty() || charge_channels <= 0
        || charge_channels > static_cast<int>(spin_density.size()))
    {
        throw std::invalid_argument("SCCS electron density requires valid charge-bearing spin channels");
    }
    const std::size_t size = spin_density[0].size();
    if (size == 0)
    {
        throw std::invalid_argument("SCCS electron density grid must not be empty");
    }

    std::vector<double> electron_density(size, 0.0);
    for (int channel = 0; channel < charge_channels; ++channel)
    {
        if (spin_density[channel].size() != size)
        {
            throw std::invalid_argument("SCCS spin-density arrays must have the same size");
        }
        for (std::size_t index = 0; index < size; ++index)
        {
            if (!std::isfinite(spin_density[channel][index]))
            {
                throw std::domain_error("SCCS electron density must be finite");
            }
            electron_density[index] += spin_density[channel][index];
        }
    }
    return electron_density;
}

ChargeDensity assemble_charge_density(const std::vector<double>& electron_density,
                                      const std::vector<double>& ionic_density,
                                      const double volume_element,
                                      const double expected_electron_count,
                                      const double expected_ionic_charge,
                                      const double normalization_tolerance,
                                      const ChargeReduction& reduction)
{
    if (electron_density.empty() || electron_density.size() != ionic_density.size())
    {
        throw std::invalid_argument("SCCS electron and ionic density arrays must have the same non-zero size");
    }
    if (!std::isfinite(volume_element) || volume_element <= 0.0
        || !std::isfinite(expected_electron_count) || expected_electron_count < 0.0
        || !std::isfinite(expected_ionic_charge) || expected_ionic_charge < 0.0
        || !std::isfinite(normalization_tolerance) || normalization_tolerance <= 0.0)
    {
        throw std::invalid_argument("SCCS charge normalization inputs must be finite and physically valid");
    }

    ChargeDensity result;
    result.electron = electron_density;
    result.ionic = ionic_density;
    result.solute.resize(electron_density.size());
    for (std::size_t index = 0; index < electron_density.size(); ++index)
    {
        if (!std::isfinite(electron_density[index]) || !std::isfinite(ionic_density[index]))
        {
            throw std::domain_error("SCCS charge densities must be finite");
        }
        result.electron_count += electron_density[index] * volume_element;
        result.ionic_charge += ionic_density[index] * volume_element;
        result.solute[index] = ionic_density[index] - electron_density[index];
    }
    reduction.reduce_sum(result.electron_count);
    reduction.reduce_sum(result.ionic_charge);
    if (!std::isfinite(result.electron_count) || !std::isfinite(result.ionic_charge))
    {
        throw std::domain_error("SCCS reduced charge integrals must be finite");
    }
    if (std::abs(result.electron_count - expected_electron_count) > normalization_tolerance)
    {
        throw std::runtime_error("SCCS electron density normalization does not match the electron count");
    }
    if (std::abs(result.ionic_charge - expected_ionic_charge) > normalization_tolerance)
    {
        throw std::runtime_error("SCCS ionic density normalization does not match the valence charge");
    }
    result.net_charge = result.ionic_charge - result.electron_count;
    return result;
}

MultipoleMoments reduced_density_moments(
    const std::vector<double>& density,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    const double volume_element,
    const ModuleBase::Vector3<double>& origin,
    const ChargeReduction& reduction)
{
    MultipoleMoments moments = density_moments(density, positions, volume_element, origin);
    double values[5] = {moments.charge,
                        moments.dipole.x,
                        moments.dipole.y,
                        moments.dipole.z,
                        moments.quadrupole_trace};
    reduction.reduce_sum(values, 5);
    moments.charge = values[0];
    moments.dipole.x = values[1];
    moments.dipole.y = values[2];
    moments.dipole.z = values[3];
    moments.quadrupole_trace = values[4];
    return moments;
}

} // namespace ModuleSccs
