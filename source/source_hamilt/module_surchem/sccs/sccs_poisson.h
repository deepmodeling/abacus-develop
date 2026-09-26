#ifndef SCCS_POISSON_H
#define SCCS_POISSON_H

#include "source_base/vector3.h"

#include <string>
#include <vector>

namespace ModuleSccs
{

struct ElectrostaticField
{
    std::vector<double> potential;
    std::vector<ModuleBase::Vector3<double>> gradient;
};

class CoulombOperator
{
  public:
    virtual ~CoulombOperator() = default;

    virtual void apply(const std::vector<double>& charge, ElectrostaticField& field) const = 0;
};

class PolarizationReduction
{
  public:
    virtual ~PolarizationReduction() = default;

    virtual void reduce_residual(double& square_sum,
                                 double& maximum,
                                 double& point_count) const = 0;
    virtual void reduce_sum(double& value) const;
};

class SerialPolarizationReduction : public PolarizationReduction
{
  public:
    void reduce_residual(double& square_sum,
                         double& maximum,
                         double& point_count) const override;
};

struct PolarizationSolverParameters
{
    int max_iterations = 0;
    std::string mixing_method = "linear";
    int mixing_history = 8;
    double mixing = 0.0;
    bool adaptive_mixing = false;
    double mixing_min = 0.1;
    double mixing_max = 0.8;
    double tolerance_rms = 0.0;
    double tolerance_max = 0.0;
};

enum class PolarizationStatus
{
    Converged,
    MaxIterations,
    NonFinite
};

struct PolarizationResult
{
    PolarizationStatus status = PolarizationStatus::MaxIterations;
    int iterations = 0;
    double residual_rms = 0.0;
    double residual_max = 0.0;
    double final_mixing = 0.0;
    int mixing_restarts = 0;
    std::vector<double> polarization_charge;
    ElectrostaticField field;
};

void validate_polarization_solver_parameters(const PolarizationSolverParameters& parameters);

PolarizationResult solve_polarization(const std::vector<double>& solute_charge,
                                      const std::vector<double>& epsilon,
                                      const std::vector<ModuleBase::Vector3<double>>& grad_log_epsilon,
                                      const std::vector<double>& initial_polarization_charge,
                                      const PolarizationSolverParameters& parameters,
                                      const CoulombOperator& coulomb);

PolarizationResult solve_polarization(const std::vector<double>& solute_charge,
                                      const std::vector<double>& epsilon,
                                      const std::vector<ModuleBase::Vector3<double>>& grad_log_epsilon,
                                      const std::vector<double>& initial_polarization_charge,
                                      const PolarizationSolverParameters& parameters,
                                      const CoulombOperator& coulomb,
                                      const PolarizationReduction& reduction);

} // namespace ModuleSccs

#endif
