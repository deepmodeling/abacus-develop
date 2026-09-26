#include "sccs_poisson.h"

#include "source_base/constants.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ModuleSccs
{
namespace
{

bool finite_vector(const ModuleBase::Vector3<double>& value)
{
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

double dot(const ModuleBase::Vector3<double>& left, const ModuleBase::Vector3<double>& right)
{
    return left.x * right.x + left.y * right.y + left.z * right.z;
}

void validate_input(const std::vector<double>& solute_charge,
                    const std::vector<double>& epsilon,
                    const std::vector<ModuleBase::Vector3<double>>& grad_log_epsilon,
                    const std::vector<double>& initial_polarization_charge)
{
    const std::size_t size = solute_charge.size();
    if (size == 0)
    {
        throw std::invalid_argument("SCCS polarization solver requires a non-empty grid");
    }
    if (epsilon.size() != size || grad_log_epsilon.size() != size)
    {
        throw std::invalid_argument("SCCS polarization input arrays must have the same size");
    }
    if (!initial_polarization_charge.empty() && initial_polarization_charge.size() != size)
    {
        throw std::invalid_argument("SCCS initial polarization charge has the wrong size");
    }
    for (std::size_t index = 0; index < size; ++index)
    {
        if (!std::isfinite(solute_charge[index]) || !std::isfinite(epsilon[index])
            || epsilon[index] < 1.0 || !finite_vector(grad_log_epsilon[index]))
        {
            throw std::domain_error("SCCS polarization inputs must be finite and epsilon must be at least one");
        }
        if (!initial_polarization_charge.empty()
            && !std::isfinite(initial_polarization_charge[index]))
        {
            throw std::domain_error("SCCS initial polarization charge must be finite");
        }
    }
}

bool valid_field(const ElectrostaticField& field, const std::size_t size)
{
    if (field.potential.size() != size || field.gradient.size() != size)
    {
        throw std::runtime_error("SCCS Coulomb operator returned arrays with the wrong size");
    }
    for (std::size_t index = 0; index < size; ++index)
    {
        if (!std::isfinite(field.potential[index]) || !finite_vector(field.gradient[index]))
        {
            return false;
        }
    }
    return true;
}

void evaluate_field(const std::vector<double>& solute_charge,
                    const std::vector<double>& polarization_charge,
                    const CoulombOperator& coulomb,
                    ElectrostaticField& field)
{
    std::vector<double> total_charge(solute_charge.size(), 0.0);
    for (std::size_t index = 0; index < solute_charge.size(); ++index)
    {
        total_charge[index] = solute_charge[index] + polarization_charge[index];
    }
    coulomb.apply(total_charge, field);
}

double global_dot(const std::vector<double>& left,
                  const std::vector<double>& right,
                  const PolarizationReduction& reduction)
{
    double value = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index)
    {
        value += left[index] * right[index];
    }
    reduction.reduce_sum(value);
    return value;
}

bool solve_dense_system(std::vector<double>& matrix,
                        std::vector<double>& right_hand_side,
                        const int dimension)
{
    for (int column = 0; column < dimension; ++column)
    {
        int pivot = column;
        double pivot_magnitude = std::abs(matrix[column * dimension + column]);
        for (int row = column + 1; row < dimension; ++row)
        {
            const double magnitude = std::abs(matrix[row * dimension + column]);
            if (magnitude > pivot_magnitude)
            {
                pivot = row;
                pivot_magnitude = magnitude;
            }
        }
        if (!std::isfinite(pivot_magnitude)
            || pivot_magnitude <= 100.0 * std::numeric_limits<double>::epsilon())
        {
            return false;
        }
        if (pivot != column)
        {
            for (int entry = column; entry < dimension; ++entry)
            {
                std::swap(matrix[column * dimension + entry],
                          matrix[pivot * dimension + entry]);
            }
            std::swap(right_hand_side[column], right_hand_side[pivot]);
        }
        const double diagonal = matrix[column * dimension + column];
        for (int row = column + 1; row < dimension; ++row)
        {
            const double factor = matrix[row * dimension + column] / diagonal;
            matrix[row * dimension + column] = 0.0;
            for (int entry = column + 1; entry < dimension; ++entry)
            {
                matrix[row * dimension + entry]
                    -= factor * matrix[column * dimension + entry];
            }
            right_hand_side[row] -= factor * right_hand_side[column];
        }
    }
    for (int row = dimension - 1; row >= 0; --row)
    {
        double value = right_hand_side[row];
        for (int column = row + 1; column < dimension; ++column)
        {
            value -= matrix[row * dimension + column] * right_hand_side[column];
        }
        const double diagonal = matrix[row * dimension + row];
        if (!std::isfinite(diagonal)
            || std::abs(diagonal) <= 100.0 * std::numeric_limits<double>::epsilon())
        {
            return false;
        }
        right_hand_side[row] = value / diagonal;
        if (!std::isfinite(right_hand_side[row]))
        {
            return false;
        }
    }
    return true;
}

void linear_mixing_step(const std::vector<double>& current,
                        const std::vector<double>& residual,
                        const double mixing,
                        std::vector<double>& next)
{
    for (std::size_t index = 0; index < current.size(); ++index)
    {
        next[index] = current[index] + mixing * residual[index];
    }
}

bool safe_coefficients(const std::vector<double>& coefficients, const int count)
{
    const double coefficient_limit = 20.0;
    for (int index = 0; index < count; ++index)
    {
        if (!std::isfinite(coefficients[index])
            || std::abs(coefficients[index]) > coefficient_limit)
        {
            return false;
        }
    }
    return true;
}

bool pulay_mixing_step(const std::vector<std::vector<double>>& values,
                       const std::vector<std::vector<double>>& residuals,
                       const double mixing,
                       const PolarizationReduction& reduction,
                       std::vector<double>& next)
{
    const int history_size = static_cast<int>(values.size());
    if (history_size < 2)
    {
        return false;
    }
    double scale = 0.0;
    std::vector<double> gram(history_size * history_size, 0.0);
    for (int row = 0; row < history_size; ++row)
    {
        for (int column = row; column < history_size; ++column)
        {
            const double value = global_dot(residuals[row], residuals[column], reduction);
            gram[row * history_size + column] = value;
            gram[column * history_size + row] = value;
        }
        scale = std::max(scale, std::abs(gram[row * history_size + row]));
    }
    if (!std::isfinite(scale) || scale <= std::numeric_limits<double>::min())
    {
        return false;
    }
    const int dimension = history_size + 1;
    std::vector<double> matrix(dimension * dimension, 0.0);
    std::vector<double> coefficients(dimension, 0.0);
    const double regularization = 1.0e-10;
    for (int row = 0; row < history_size; ++row)
    {
        for (int column = 0; column < history_size; ++column)
        {
            matrix[row * dimension + column]
                = gram[row * history_size + column] / scale;
        }
        matrix[row * dimension + row] += regularization;
        matrix[row * dimension + history_size] = 1.0;
        matrix[history_size * dimension + row] = 1.0;
    }
    coefficients[history_size] = 1.0;
    if (!solve_dense_system(matrix, coefficients, dimension)
        || !safe_coefficients(coefficients, history_size))
    {
        return false;
    }
    std::fill(next.begin(), next.end(), 0.0);
    for (int history = 0; history < history_size; ++history)
    {
        for (std::size_t index = 0; index < next.size(); ++index)
        {
            next[index] += coefficients[history]
                           * (values[history][index] + mixing * residuals[history][index]);
        }
    }
    return true;
}

bool anderson_mixing_step(const std::vector<std::vector<double>>& values,
                          const std::vector<std::vector<double>>& residuals,
                          const double mixing,
                          const PolarizationReduction& reduction,
                          std::vector<double>& next)
{
    const int history_size = static_cast<int>(values.size());
    if (history_size < 2)
    {
        return false;
    }
    const int difference_count = history_size - 1;
    std::vector<std::vector<double>> value_differences(
        difference_count, std::vector<double>(next.size(), 0.0));
    std::vector<std::vector<double>> residual_differences(
        difference_count, std::vector<double>(next.size(), 0.0));
    for (int history = 0; history < difference_count; ++history)
    {
        for (std::size_t index = 0; index < next.size(); ++index)
        {
            value_differences[history][index]
                = values[history + 1][index] - values[history][index];
            residual_differences[history][index]
                = residuals[history + 1][index] - residuals[history][index];
        }
    }
    double scale = 0.0;
    std::vector<double> matrix(difference_count * difference_count, 0.0);
    for (int row = 0; row < difference_count; ++row)
    {
        for (int column = row; column < difference_count; ++column)
        {
            const double value = global_dot(residual_differences[row],
                                            residual_differences[column],
                                            reduction);
            matrix[row * difference_count + column] = value;
            matrix[column * difference_count + row] = value;
        }
        scale = std::max(scale, std::abs(matrix[row * difference_count + row]));
    }
    if (!std::isfinite(scale) || scale <= std::numeric_limits<double>::min())
    {
        return false;
    }
    std::vector<double> coefficients(difference_count, 0.0);
    const std::vector<double>& current_residual = residuals.back();
    const double regularization = 1.0e-10;
    for (int row = 0; row < difference_count; ++row)
    {
        coefficients[row]
            = global_dot(residual_differences[row], current_residual, reduction) / scale;
        for (int column = 0; column < difference_count; ++column)
        {
            matrix[row * difference_count + column] /= scale;
        }
        matrix[row * difference_count + row] += regularization;
    }
    if (!solve_dense_system(matrix, coefficients, difference_count)
        || !safe_coefficients(coefficients, difference_count))
    {
        return false;
    }
    linear_mixing_step(values.back(), current_residual, mixing, next);
    for (int history = 0; history < difference_count; ++history)
    {
        for (std::size_t index = 0; index < next.size(); ++index)
        {
            next[index]
                -= coefficients[history]
                   * (value_differences[history][index]
                      + mixing * residual_differences[history][index]);
        }
    }
    return true;
}

bool finite_vector(const std::vector<double>& values)
{
    for (std::size_t index = 0; index < values.size(); ++index)
    {
        if (!std::isfinite(values[index]))
        {
            return false;
        }
    }
    return true;
}

struct AdaptiveMixingState
{
    double previous_residual_rms = std::numeric_limits<double>::infinity();
    double current_mixing = 0.0;
    int rapid_decrease_count = 0;
    int residual_rise_count = 0;
};

bool update_adaptive_mixing(
    const PolarizationSolverParameters& parameters,
    const double residual_rms,
    AdaptiveMixingState& state,
    std::vector<std::vector<double>>& value_history,
    std::vector<std::vector<double>>& residual_history,
    int& restart_count)
{
    if (!parameters.adaptive_mixing
        || !std::isfinite(state.previous_residual_rms))
    {
        return false;
    }
    const double residual_ratio = residual_rms / state.previous_residual_rms;
    if (residual_ratio < 0.7)
    {
        ++state.rapid_decrease_count;
        state.residual_rise_count = 0;
        if (state.rapid_decrease_count >= 3)
        {
            state.current_mixing = std::min(parameters.mixing_max,
                                            1.1 * state.current_mixing);
            state.rapid_decrease_count = 0;
        }
        return false;
    }
    state.rapid_decrease_count = 0;
    if (residual_ratio <= 1.1)
    {
        state.residual_rise_count = 0;
        return false;
    }
    ++state.residual_rise_count;
    if (residual_ratio <= 2.0 && state.residual_rise_count < 2)
    {
        return false;
    }
    state.current_mixing = std::max(parameters.mixing_min,
                                    0.5 * state.current_mixing);
    value_history.clear();
    residual_history.clear();
    ++restart_count;
    state.residual_rise_count = 0;
    return true;
}

bool valid_mixing_controls(const PolarizationSolverParameters& parameters)
{
    const bool mixing_valid = std::isfinite(parameters.mixing)
                              && parameters.mixing > 0.0
                              && parameters.mixing <= 1.0;
    const bool bounds_valid = std::isfinite(parameters.mixing_min)
                              && std::isfinite(parameters.mixing_max)
                              && parameters.mixing_min > 0.0
                              && parameters.mixing_min <= parameters.mixing_max
                              && parameters.mixing_max <= 1.0;
    const bool initial_value_valid = !parameters.adaptive_mixing
                                     || (parameters.mixing >= parameters.mixing_min
                                         && parameters.mixing <= parameters.mixing_max);
    return mixing_valid && bounds_valid && initial_value_valid;
}

} // namespace

void PolarizationReduction::reduce_sum(double& value) const
{
    if (!std::isfinite(value))
    {
        throw std::domain_error("SCCS reduction input must be finite");
    }
}

void validate_polarization_solver_parameters(const PolarizationSolverParameters& parameters)
{
    if (parameters.max_iterations <= 0)
    {
        throw std::invalid_argument("SCCS polarization maximum iteration count must be positive");
    }
    if (parameters.mixing_method != "linear" && parameters.mixing_method != "pulay"
        && parameters.mixing_method != "anderson")
    {
        throw std::invalid_argument("unknown SCCS polarization mixing method: "
                                    + parameters.mixing_method);
    }
    if (parameters.mixing_history < 2)
    {
        throw std::invalid_argument("SCCS accelerated-mixing history must be at least two");
    }
    if (!valid_mixing_controls(parameters))
    {
        throw std::invalid_argument(
            "SCCS mixing must satisfy 0 < min <= initial <= max <= 1 when adaptive mixing is enabled");
    }
    if (!std::isfinite(parameters.tolerance_rms) || parameters.tolerance_rms <= 0.0
        || !std::isfinite(parameters.tolerance_max) || parameters.tolerance_max <= 0.0)
    {
        throw std::invalid_argument("SCCS polarization residual tolerances must be positive and finite");
    }
}

void SerialPolarizationReduction::reduce_residual(double& square_sum,
                                                  double& maximum,
                                                  double& point_count) const
{
    if (!std::isfinite(square_sum) || !std::isfinite(maximum) || !std::isfinite(point_count))
    {
        throw std::domain_error("SCCS residual values must be finite before reduction");
    }
}

PolarizationResult solve_polarization(
    const std::vector<double>& solute_charge,
    const std::vector<double>& epsilon,
    const std::vector<ModuleBase::Vector3<double>>& grad_log_epsilon,
    const std::vector<double>& initial_polarization_charge,
    const PolarizationSolverParameters& parameters,
    const CoulombOperator& coulomb)
{
    const SerialPolarizationReduction reduction;
    return solve_polarization(solute_charge,
                              epsilon,
                              grad_log_epsilon,
                              initial_polarization_charge,
                              parameters,
                              coulomb,
                              reduction);
}

PolarizationResult solve_polarization(
    const std::vector<double>& solute_charge,
    const std::vector<double>& epsilon,
    const std::vector<ModuleBase::Vector3<double>>& grad_log_epsilon,
    const std::vector<double>& initial_polarization_charge,
    const PolarizationSolverParameters& parameters,
    const CoulombOperator& coulomb,
    const PolarizationReduction& reduction)
{
    validate_polarization_solver_parameters(parameters);
    validate_input(solute_charge, epsilon, grad_log_epsilon, initial_polarization_charge);

    const std::size_t size = solute_charge.size();
    PolarizationResult result;
    result.polarization_charge.assign(size, 0.0);
    result.final_mixing = parameters.mixing;
    if (!initial_polarization_charge.empty())
    {
        result.polarization_charge = initial_polarization_charge;
    }

    std::vector<double> trial(size, 0.0);
    std::vector<double> residual(size, 0.0);
    std::vector<double> next(size, 0.0);
    std::vector<std::vector<double>> value_history;
    std::vector<std::vector<double>> residual_history;
    AdaptiveMixingState mixing_state;
    mixing_state.current_mixing = parameters.mixing;
    for (int iteration = 1; iteration <= parameters.max_iterations; ++iteration)
    {
        evaluate_field(solute_charge, result.polarization_charge, coulomb, result.field);
        result.iterations = iteration;
        if (!valid_field(result.field, size))
        {
            result.status = PolarizationStatus::NonFinite;
            return result;
        }

        double residual_square_sum = 0.0;
        result.residual_max = 0.0;
        for (std::size_t index = 0; index < size; ++index)
        {
            const double dielectric_source
                = dot(grad_log_epsilon[index], result.field.gradient[index])
                  / ModuleBase::FOUR_PI;
            const double screening_source
                = -(epsilon[index] - 1.0) * solute_charge[index] / epsilon[index];
            trial[index] = dielectric_source + screening_source;
            residual[index] = trial[index] - result.polarization_charge[index];
            if (!std::isfinite(trial[index]) || !std::isfinite(residual[index]))
            {
                result.status = PolarizationStatus::NonFinite;
                return result;
            }
            residual_square_sum += residual[index] * residual[index];
            result.residual_max = std::max(result.residual_max, std::abs(residual[index]));
        }
        double point_count = static_cast<double>(size);
        reduction.reduce_residual(residual_square_sum, result.residual_max, point_count);
        if (!std::isfinite(residual_square_sum) || !std::isfinite(result.residual_max)
            || !std::isfinite(point_count) || point_count <= 0.0)
        {
            result.status = PolarizationStatus::NonFinite;
            return result;
        }
        result.residual_rms = std::sqrt(residual_square_sum / point_count);

        if (result.residual_rms <= parameters.tolerance_rms
            && result.residual_max <= parameters.tolerance_max)
        {
            result.polarization_charge = trial;
            result.status = PolarizationStatus::Converged;
            evaluate_field(solute_charge, result.polarization_charge, coulomb, result.field);
            if (!valid_field(result.field, size))
            {
                result.status = PolarizationStatus::NonFinite;
            }
            return result;
        }

        const bool accelerated = parameters.mixing_method != "linear";
        const bool force_linear_step
            = update_adaptive_mixing(parameters,
                                     result.residual_rms,
                                     mixing_state,
                                     value_history,
                                     residual_history,
                                     result.mixing_restarts);
        if (!parameters.adaptive_mixing && accelerated
            && result.residual_rms > 4.0 * mixing_state.previous_residual_rms)
        {
            value_history.clear();
            residual_history.clear();
        }
        value_history.push_back(result.polarization_charge);
        residual_history.push_back(residual);
        if (static_cast<int>(value_history.size()) > parameters.mixing_history)
        {
            value_history.erase(value_history.begin());
            residual_history.erase(residual_history.begin());
        }

        bool mixed = false;
        if (!force_linear_step && parameters.mixing_method == "pulay")
        {
            mixed = pulay_mixing_step(value_history,
                                      residual_history,
                                      mixing_state.current_mixing,
                                      reduction,
                                      next);
        }
        else if (!force_linear_step && parameters.mixing_method == "anderson")
        {
            mixed = anderson_mixing_step(value_history,
                                         residual_history,
                                         mixing_state.current_mixing,
                                         reduction,
                                         next);
        }
        if (!mixed || !finite_vector(next))
        {
            linear_mixing_step(result.polarization_charge,
                               residual,
                               mixing_state.current_mixing,
                               next);
            if (accelerated && !finite_vector(next))
            {
                result.status = PolarizationStatus::NonFinite;
                return result;
            }
        }
        result.polarization_charge.swap(next);
        result.final_mixing = mixing_state.current_mixing;
        mixing_state.previous_residual_rms = result.residual_rms;
    }

    result.status = PolarizationStatus::MaxIterations;
    evaluate_field(solute_charge, result.polarization_charge, coulomb, result.field);
    if (!valid_field(result.field, size))
    {
        result.status = PolarizationStatus::NonFinite;
    }
    return result;
}

} // namespace ModuleSccs
