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

// Local-rank measurements of the periodic part of the Coulomb operator.
// Transform time includes PW packing and MPI transposes, not just FFT kernels.
struct CoulombTransformProfile
{
    long forward_calls = 0;
    long inverse_calls = 0;
    double forward_seconds = 0.0;
    double inverse_seconds = 0.0;
    double other_seconds = 0.0;
};

class CoulombOperator
{
  public:
    virtual ~CoulombOperator() = default;

    virtual CoulombTransformProfile transform_profile() const
    {
        return CoulombTransformProfile();
    }

    virtual void apply(const std::vector<double>& charge, ElectrostaticField& field) const = 0;

    // Scalar-only clients avoid gradient transforms when the operator supports it.
    virtual void apply_potential(const std::vector<double>& charge,
                                 std::vector<double>& potential) const
    {
        ElectrostaticField field;
        apply(charge, field);
        potential.swap(field.potential);
    }

    // Polarization iterations need only the gradient. Operators without a
    // specialized implementation can use the complete-field fallback.
    virtual void apply_gradient(const std::vector<double>& charge,
                                std::vector<ModuleBase::Vector3<double>>& gradient) const
    {
        ElectrostaticField field;
        apply(charge, field);
        gradient.swap(field.gradient);
    }

    // Adjoint of charge -> electrostatic field gradient under the grid inner product.
    virtual void apply_gradient_adjoint(
        const std::vector<ModuleBase::Vector3<double>>& field,
        std::vector<double>& result) const = 0;

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
    // Periodic sqrt-CG only: verify v = P(q - K v) after convergence.
    bool check_fixed_point = false;
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
    bool warm_started = false;
    bool fixed_point_checked = false;
    double fixed_point_defect_rms = 0.0;
    double fixed_point_defect_max = 0.0;
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
