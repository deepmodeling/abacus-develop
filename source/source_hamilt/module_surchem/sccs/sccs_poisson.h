#ifndef SCCS_POISSON_H
#define SCCS_POISSON_H

#include "source_base/vector3.h"

#include <vector>

namespace ModuleSccs
{

struct ElectrostaticField
{
    std::vector<double> potential;
    std::vector<ModuleBase::Vector3<double>> gradient;
};

// Local-rank FFT counts of the periodic part of the Coulomb operator.
// Transform time is reported by the PW_Basis timers.
struct CoulombTransformCounts
{
    long forward_calls = 0;
    long inverse_calls = 0;
};

class CoulombOperator
{
  public:
    virtual ~CoulombOperator() = default;

    virtual CoulombTransformCounts transform_counts() const
    {
        return CoulombTransformCounts();
    }

    // True when an analytic open-boundary (PCC) term is added. Such potentials
    // keep their physical gauge instead of the periodic zero mean (ENVIRON's
    // has_corrections criterion in generalized_sqrt).
    virtual bool has_boundary_correction() const
    {
        return false;
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

    // Gradient-only clients. Operators without a specialized implementation
    // can use the complete-field fallback.
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

// Controls of the sqrt-preconditioned CG (solve_chain_sccs_response).
struct PolarizationSolverParameters
{
    int max_iterations = 0;
    double tolerance_rms = 0.0;
    double tolerance_max = 0.0;
    // Verify v = P(q - K v) after convergence (one extra Poisson solve).
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
    bool warm_started = false;
    bool fixed_point_checked = false;
    double fixed_point_defect_rms = 0.0;
    double fixed_point_defect_max = 0.0;
    std::vector<double> polarization_charge;
    ElectrostaticField field;
};

} // namespace ModuleSccs

#endif
