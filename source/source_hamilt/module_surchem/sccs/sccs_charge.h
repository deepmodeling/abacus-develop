#ifndef SCCS_CHARGE_H
#define SCCS_CHARGE_H

#include "../pcc/sccs_pcc.h"

#include <vector>

namespace ModuleSccs
{

class ChargeReduction
{
  public:
    virtual ~ChargeReduction() = default;

    virtual void reduce_sum(double& value) const = 0;
    virtual void reduce_sum(double* values, int count) const = 0;
};

class SerialChargeReduction : public ChargeReduction
{
  public:
    void reduce_sum(double& value) const override;
    void reduce_sum(double* values, int count) const override;
};

struct ChargeDensity
{
    std::vector<double> electron;
    std::vector<double> ionic;
    std::vector<double> solute;
    double electron_count = 0.0;
    double ionic_charge = 0.0;
    double net_charge = 0.0;
};

std::vector<double> sum_electron_density(const std::vector<std::vector<double>>& spin_density,
                                         int charge_channels);

ChargeDensity assemble_charge_density(const std::vector<double>& electron_density,
                                      const std::vector<double>& ionic_density,
                                      double volume_element,
                                      double expected_electron_count,
                                      double expected_ionic_charge,
                                      double normalization_tolerance,
                                      const ChargeReduction& reduction);

MultipoleMoments reduced_density_moments(
    const std::vector<double>& density,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    double volume_element,
    const ModuleBase::Vector3<double>& origin,
    const ChargeReduction& reduction);

} // namespace ModuleSccs

#endif
