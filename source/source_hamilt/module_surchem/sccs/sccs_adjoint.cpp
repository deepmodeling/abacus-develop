#include "sccs_adjoint.h"

#include "sccs_functional.h"
#include "sccs_periodic.h"
#include "sccs_pw_coulomb.h"
#include "source_base/constants.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ModuleSccs
{
namespace
{

double reduced_dot(const std::vector<double>& left,
                   const std::vector<double>& right,
                   const PolarizationReduction& reduction)
{
    double value = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i)
    {
        value += left[i] * right[i];
    }
    reduction.reduce_sum(value);
    return value;
}

bool converged_residual(const std::vector<double>& residual,
                        const PolarizationSolverParameters& parameters,
                        const PolarizationReduction& reduction,
                        AdjointResult& result)
{
    double square_sum = 0.0;
    double maximum = 0.0;
    double count = static_cast<double>(residual.size());
    for (const double value : residual)
    {
        square_sum += value * value;
        maximum = std::max(maximum, std::abs(value));
    }
    reduction.reduce_residual(square_sum, maximum, count);
    if (!std::isfinite(square_sum) || !std::isfinite(maximum) || count <= 0.0)
    {
        throw std::runtime_error("SCCS adjoint residual is not finite");
    }
    result.residual_rms = std::sqrt(square_sum / count);
    result.residual_max = maximum;
    return result.residual_rms <= parameters.tolerance_rms
           && result.residual_max <= parameters.tolerance_max;
}

// For total effective charge t: L t = rho/epsilon, with
// L = I - diag(grad(log(epsilon))/(4*pi)) dot grad(C).
// The adjoint is L^T = I - grad(C)^T diag(grad(log(epsilon))/(4*pi)).
class AdjointOperator
{
  public:
    AdjointOperator(const PeriodicSccsResult& response, const CoulombOperator& coulomb)
        : response_(response), coulomb_(coulomb)
    {
    }

    void apply(const std::vector<double>& values,
               std::vector<ModuleBase::Vector3<double>>& weighted,
               std::vector<double>& result) const
    {
        weighted.resize(values.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            for (int d = 0; d < 3; ++d)
            {
                weighted[i][d] = response_.grad_log_epsilon[i][d] * values[i]
                                 / ModuleBase::FOUR_PI;
            }
        }
        coulomb_.apply_gradient_adjoint(weighted, result);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            result[i] = values[i] - result[i];
        }
    }

  private:
    const PeriodicSccsResult& response_;
    const CoulombOperator& coulomb_;
};

AdjointResult solve_adjoint(const std::vector<double>& rhs,
                            const std::vector<double>& initial,
                            const AdjointOperator& op,
                            const PolarizationSolverParameters& parameters,
                            const PolarizationReduction& reduction)
{
    AdjointResult result;
    result.potential = initial;
    const std::size_t size = rhs.size();
    // Workspace belongs to this solve, so repeated applications reuse storage
    // without introducing mutable state into the Coulomb operator.
    std::vector<ModuleBase::Vector3<double>> weighted(size);
    std::vector<double> residual;
    op.apply(result.potential, weighted, residual);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        residual[i] = rhs[i] - residual[i];
    }
    std::vector<double> shadow = residual;
    std::vector<double> direction(size, 0.0);
    std::vector<double> image(size, 0.0);
    std::vector<double> intermediate(size);
    std::vector<double> intermediate_image;
    double previous_rho = 1.0;
    double alpha = 1.0;
    double omega = 1.0;
    // BiCGSTAB solves the transpose system, which is not symmetric. Always
    // verify the true residual before accepting an iterated residual estimate.
    for (int iteration = 0; iteration <= parameters.max_iterations; ++iteration)
    {
        result.iterations = iteration;
        if (converged_residual(residual, parameters, reduction, result))
        {
            op.apply(result.potential, weighted, residual);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for (std::size_t i = 0; i < size; ++i)
            {
                residual[i] = rhs[i] - residual[i];
            }
            if (converged_residual(residual, parameters, reduction, result))
            {
                return result;
            }
            shadow = residual;
            std::fill(direction.begin(), direction.end(), 0.0);
            std::fill(image.begin(), image.end(), 0.0);
            previous_rho = 1.0;
            alpha = 1.0;
            omega = 1.0;
        }
        if (iteration == parameters.max_iterations)
        {
            break;
        }
        const double rho = reduced_dot(shadow, residual, reduction);
        if (!std::isfinite(rho) || rho == 0.0 || omega == 0.0)
        {
            throw std::runtime_error("SCCS adjoint BiCGSTAB breakdown");
        }
        const double beta = (rho / previous_rho) * (alpha / omega);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (std::size_t i = 0; i < size; ++i)
        {
            direction[i] = residual[i] + beta * (direction[i] - omega * image[i]);
        }
        op.apply(direction, weighted, image);
        const double denominator = reduced_dot(shadow, image, reduction);
        if (!std::isfinite(denominator) || denominator == 0.0)
        {
            throw std::runtime_error("SCCS adjoint BiCGSTAB singular direction");
        }
        alpha = rho / denominator;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (std::size_t i = 0; i < size; ++i)
        {
            intermediate[i] = residual[i] - alpha * image[i];
            result.potential[i] += alpha * direction[i];
        }
        if (converged_residual(intermediate, parameters, reduction, result))
        {
            residual = intermediate;
            previous_rho = rho;
            continue;
        }
        op.apply(intermediate, weighted, intermediate_image);
        const double image_norm = reduced_dot(intermediate_image, intermediate_image, reduction);
        if (!std::isfinite(image_norm) || image_norm == 0.0)
        {
            throw std::runtime_error("SCCS adjoint BiCGSTAB singular stabilization");
        }
        const double projection = reduced_dot(intermediate_image, intermediate, reduction);
        omega = projection / image_norm;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (std::size_t i = 0; i < size; ++i)
        {
            result.potential[i] += omega * intermediate[i];
            residual[i] = intermediate[i] - omega * intermediate_image[i];
        }
        previous_rho = rho;
    }
    throw std::runtime_error("SCCS discrete adjoint iteration did not converge");
}

} // namespace

// Differentiate the converged discrete polarization equation, not its continuum
// limit. For L t = q/epsilon, solve L^T lambda = Cq + 2*c, where c is the
// ionic-shape coefficient. Return dE/dq for ionic forces and dE/dn (including
// dq/dn = -1) for the electronic Hamiltonian. See discrete_derivative.md.
AdjointResult evaluate_discrete_electrostatic_derivative(
    const std::vector<double>& cavity_density,
    const CavityParameters& cavity_parameters,
    const std::vector<double>& solute_charge,
    const PeriodicSccsResult& response,
    const ElectrostaticField& vacuum_field,
    const ModulePW::PW_Basis& basis,
    const double tpiba,
    const double ionic_shape_coefficient,
    const PolarizationSolverParameters& parameters,
    const std::vector<double>& initial_adjoint,
    const CoulombOperator& coulomb,
    const PolarizationReduction& reduction,
    ElectrostaticFunctionalResult& functional)
{
    validate_cavity_parameters(cavity_parameters);
    const std::size_t size = solute_charge.size();
    if (size == 0 || cavity_density.size() != size || response.epsilon.size() != size
        || response.depsilon_drho.size() != size
        || response.density_gradient.size() != size
        || response.grad_log_epsilon.size() != size || vacuum_field.potential.size() != size
        || response.polarization.field.potential.size() != size
        || response.polarization.field.gradient.size() != size
        || (!initial_adjoint.empty() && initial_adjoint.size() != size)
        || !std::isfinite(ionic_shape_coefficient))
    {
        throw std::invalid_argument("SCCS adjoint input arrays have inconsistent sizes");
    }
    validate_polarization_solver_parameters(parameters);
    std::vector<double> rhs(size);
    std::vector<double> initial(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        rhs[i] = vacuum_field.potential[i] + 2.0 * ionic_shape_coefficient;
        initial[i] = initial_adjoint.empty()
                         ? response.epsilon[i] * response.polarization.field.potential[i]
                         : initial_adjoint[i];
    }
    const AdjointOperator op(response, coulomb);
    AdjointResult result = solve_adjoint(rhs, initial, op, parameters, reduction);
    // The forward response already evaluated Dn on this grid. Reuse the raw
    // gradient; recovering it from grad(log(epsilon)) would divide by zero
    // outside the cavity transition and lose the original discrete field.
    const std::vector<ModuleBase::Vector3<double>>& density_gradient = response.density_gradient;
    const double density_ratio = cavity_parameters.density_max / cavity_parameters.density_min;
    const double log_width = std::log(density_ratio);
    const double log_bulk = std::log(cavity_parameters.epsilon_bulk);
    std::vector<ModuleBase::Vector3<double>> weighted_gradient(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        for (int d = 0; d < 3; ++d)
        {
            weighted_gradient[i][d]
                = result.potential[i] * response.polarization.field.gradient[i][d]
                  * response.depsilon_drho[i] / response.epsilon[i];
        }
    }
    // For g(n) = f'(n) Dn, delta g = f'(n) D(delta n) + f''(n) Dn delta n.
    // The first term uses D^T (negative divergence); the second is local below.
    const std::vector<double> negative_divergence
        = periodic_negative_divergence(weighted_gradient, basis, tpiba);
    functional.charge_potential.resize(size);
    functional.electron_potential.resize(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        const double epsilon = response.epsilon[i];
        const double charge_derivative
            = 0.5 * (response.polarization.field.potential[i] + result.potential[i] / epsilon)
              - vacuum_field.potential[i] - ionic_shape_coefficient;
        double second_log_derivative = 0.0;
        const double density = cavity_density[i];
        if (density > cavity_parameters.density_min && density < cavity_parameters.density_max)
        {
            const double local_density_ratio = cavity_parameters.density_max / density;
            const double x = std::log(local_density_ratio) / log_width;
            const double angle = ModuleBase::TWO_PI * x;
            second_log_derivative
                = log_bulk * ((1.0 - std::cos(angle)) / log_width
                              + ModuleBase::TWO_PI * std::sin(angle) / (log_width * log_width))
                  / (density * density);
        }
        double gradient_dot = 0.0;
        for (int d = 0; d < 3; ++d)
        {
            gradient_dot += density_gradient[i][d] * response.polarization.field.gradient[i][d];
        }
        const double cavity_derivative
            = -0.5 * result.potential[i] * solute_charge[i] * response.depsilon_drho[i]
                  / (epsilon * epsilon)
              + (negative_divergence[i] + result.potential[i] * second_log_derivative * gradient_dot)
                    / (2.0 * ModuleBase::FOUR_PI);
        functional.charge_potential[i] = charge_derivative;
        functional.electron_potential[i] = -charge_derivative + cavity_derivative;
    }
    return result;
}

} // namespace ModuleSccs
