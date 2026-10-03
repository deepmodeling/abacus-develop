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

// Coulomb potential of a grid charge, periodic or with an open-boundary term.
// The sqrt-CG preconditioner and the vacuum reference need only the scalar
// potential; field gradients are taken from the solved potential.
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

    virtual void apply_potential(const std::vector<double>& charge,
                                 std::vector<double>& potential) const = 0;
};

// Controls of the sqrt-preconditioned CG (solve_sccs_response).
struct PolarizationSolverParameters
{
    int max_iterations = 0;
    double tolerance_rms = 0.0;
    double tolerance_max = 0.0;
    // Verify v = P(q - K v) after convergence (one extra Poisson solve).
    bool check_fixed_point = false;
};

struct PolarizationResult
{
    int iterations = 0;
    double residual_rms = 0.0;
    double residual_max = 0.0;
    bool warm_started = false;
    bool fixed_point_checked = false;
    double fixed_point_defect_rms = 0.0;
    double fixed_point_defect_max = 0.0;
    // Open boundaries only: the ENVIRON dielectric_of_potential polarization
    // density, for the PCC moment diagnostics. Empty when periodic.
    std::vector<double> polarization_charge;
    ElectrostaticField field;
};

} // namespace ModuleSccs

#endif
