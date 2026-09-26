#ifndef SCCS_ADJOINT_H
#define SCCS_ADJOINT_H

#include <vector>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{
class CoulombOperator;
class PolarizationReduction;
struct ElectrostaticField;
struct ElectrostaticFunctionalResult;
struct PeriodicSccsResult;
struct PolarizationSolverParameters;

struct AdjointResult
{
    std::vector<double> potential;
    int iterations = 0;
    double residual_rms = 0.0;
    double residual_max = 0.0;
};

// Differentiate the discrete polarization fixed point, including a possible
// ionic-shape energy coefficient times the integrated polarization charge.
AdjointResult evaluate_discrete_electrostatic_derivative(
    const std::vector<double>& solute_charge,
    const PeriodicSccsResult& response,
    const ElectrostaticField& vacuum_field,
    const ModulePW::PW_Basis& basis,
    double tpiba,
    double ionic_shape_coefficient,
    const PolarizationSolverParameters& parameters,
    const std::vector<double>& initial_adjoint,
    const CoulombOperator& coulomb,
    const PolarizationReduction& reduction,
    ElectrostaticFunctionalResult& functional);

} // namespace ModuleSccs

#endif
