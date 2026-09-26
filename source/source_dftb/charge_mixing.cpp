#include "source_dftb/charge_mixing.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace ModuleDFTB
{
namespace
{

double vector_dot(const std::vector<double>& lhs, const std::vector<double>& rhs)
{
    double result = 0.0;
    for (std::size_t i = 0; i < lhs.size(); ++i) result += lhs[i] * rhs[i];
    return result;
}

bool solve_dense_system(std::vector<std::vector<double>> matrix,
                        std::vector<double> rhs,
                        std::vector<double>* solution)
{
    if (solution == nullptr || matrix.empty() || matrix.size() != rhs.size()) return false;
    const std::size_t n = matrix.size();
    double matrix_scale = 0.0;
    for (const auto& row : matrix)
    {
        if (row.size() != n) return false;
        for (const double value : row)
        {
            if (!std::isfinite(value)) return false;
            matrix_scale = std::max(matrix_scale, std::abs(value));
        }
    }
    if (!(matrix_scale > 0.0) || !std::isfinite(matrix_scale)) return false;
    for (const double value : rhs)
        if (!std::isfinite(value)) return false;

    for (std::size_t column = 0; column < n; ++column)
    {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < n; ++row)
            if (std::abs(matrix[row][column]) > std::abs(matrix[pivot][column])) pivot = row;
        if (std::abs(matrix[pivot][column]) < 1.0e-14 * matrix_scale) return false;
        std::swap(matrix[pivot], matrix[column]);
        std::swap(rhs[pivot], rhs[column]);
        const double pivot_value = matrix[column][column];
        for (std::size_t j = column; j < n; ++j) matrix[column][j] /= pivot_value;
        rhs[column] /= pivot_value;
        for (std::size_t row = 0; row < n; ++row)
        {
            if (row == column) continue;
            const double factor = matrix[row][column];
            for (std::size_t j = column; j < n; ++j) matrix[row][j] -= factor * matrix[column][j];
            rhs[row] -= factor * rhs[column];
        }
    }
    for (const double value : rhs)
        if (!std::isfinite(value) || std::abs(value) > 1.0e8) return false;
    *solution = std::move(rhs);
    return true;
}

bool pulay_coefficients(const std::vector<std::vector<double>>& residuals,
                        std::vector<double>* coefficients)
{
    const std::size_t count = residuals.size();
    if (count < 2 || coefficients == nullptr || residuals.front().empty()) return false;
    const std::size_t n_atoms = residuals.front().size();
    std::vector<std::vector<double>> matrix(count + 1, std::vector<double>(count + 1, 0.0));
    double diagonal_scale = 0.0;
    for (std::size_t i = 0; i < count; ++i)
    {
        if (residuals[i].size() != n_atoms) return false;
        for (std::size_t j = 0; j < count; ++j)
            matrix[i][j] = vector_dot(residuals[i], residuals[j]) / static_cast<double>(n_atoms);
        diagonal_scale = std::max(diagonal_scale, matrix[i][i]);
    }
    if (!(diagonal_scale > 0.0) || !std::isfinite(diagonal_scale)) return false;
    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t j = 0; j < count; ++j) matrix[i][j] /= diagonal_scale;
        matrix[i][i] += 1.0e-10;
        matrix[i][count] = 1.0;
        matrix[count][i] = 1.0;
    }
    std::vector<double> rhs(count + 1, 0.0);
    rhs[count] = 1.0;
    std::vector<double> solution;
    if (!solve_dense_system(std::move(matrix), std::move(rhs), &solution)) return false;
    coefficients->assign(solution.begin(), solution.begin() + count);
    double coefficient_sum = 0.0;
    for (const double coefficient : *coefficients)
    {
        if (std::abs(coefficient) > 20.0) return false;
        coefficient_sum += coefficient;
    }
    return std::isfinite(coefficient_sum) && std::abs(coefficient_sum - 1.0) <= 1.0e-8;
}

void restore_total_charge(std::vector<double>* charges, const double target)
{
    double current = 0.0;
    for (const double charge : *charges) current += charge;
    const double correction = (target - current) / static_cast<double>(charges->size());
    for (double& charge : *charges) charge += correction;
}

std::vector<double> linear_mix(const std::vector<double>& charges,
                               const std::vector<double>& residual,
                               const double parameter)
{
    std::vector<double> mixed(charges.size(), 0.0);
    for (std::size_t i = 0; i < charges.size(); ++i) mixed[i] = charges[i] + parameter * residual[i];
    return mixed;
}

} // namespace

DftbChargeMixer::DftbChargeMixer(const DftbChargeMixerParameters& parameters) : parameters_(parameters)
{
    if ((parameters_.method != "linear" && parameters_.method != "pulay" && parameters_.method != "broyden")
        || !(parameters_.mixing_parameter > 0.0 && parameters_.mixing_parameter <= 1.0)
        || parameters_.history < 2 || parameters_.history > 20
        || !(parameters_.inverse_jacobi_weight > 0.0) || !std::isfinite(parameters_.inverse_jacobi_weight)
        || !(parameters_.minimal_weight > 0.0) || !std::isfinite(parameters_.minimal_weight)
        || !(parameters_.maximal_weight >= parameters_.minimal_weight) || !std::isfinite(parameters_.maximal_weight)
        || !(parameters_.weight_factor > 0.0) || !std::isfinite(parameters_.weight_factor))
        throw std::invalid_argument("Invalid native DFTB charge-mixer controls");
}

std::vector<double> DftbChargeMixer::mix(const std::vector<double>& charges,
                                         const std::vector<double>& residual,
                                         const double target_total_charge)
{
    if (charges.empty() || charges.size() != residual.size() || !std::isfinite(target_total_charge))
        throw std::invalid_argument("Native DFTB mixer needs finite total charge and equally sized, non-empty vectors");
    for (std::size_t i = 0; i < charges.size(); ++i)
        if (!std::isfinite(charges[i]) || !std::isfinite(residual[i]))
            throw std::invalid_argument("Native DFTB mixer received a non-finite charge or residual");

    std::vector<double> next;
    bool accelerated = false;
    bool acceleration_history_ready = false;
    if (parameters_.method == "pulay")
    {
        charge_history_.push_back(charges);
        residual_history_.push_back(residual);
        while (charge_history_.size() > static_cast<std::size_t>(parameters_.history))
        {
            charge_history_.erase(charge_history_.begin());
            residual_history_.erase(residual_history_.begin());
        }
        acceleration_history_ready = residual_history_.size() >= 2;
        std::vector<double> coefficients;
        if (pulay_coefficients(residual_history_, &coefficients))
        {
            next.assign(charges.size(), 0.0);
            for (std::size_t history = 0; history < coefficients.size(); ++history)
            {
                const std::vector<double> mixed = linear_mix(charge_history_[history], residual_history_[history],
                                                              parameters_.mixing_parameter);
                for (std::size_t atom = 0; atom < charges.size(); ++atom)
                    next[atom] += coefficients[history] * mixed[atom];
            }
            accelerated = std::all_of(next.begin(), next.end(), [](double value) { return std::isfinite(value); });
        }
    }
    else if (parameters_.method == "broyden")
    {
        if (has_previous_)
        {
            std::vector<double> difference(residual.size(), 0.0);
            std::vector<double> charge_step(charges.size(), 0.0);
            double difference_norm_squared = 0.0;
            double residual_norm_squared = 0.0;
            for (std::size_t i = 0; i < residual.size(); ++i)
            {
                difference[i] = residual[i] - previous_residual_[i];
                charge_step[i] = charges[i] - previous_charges_[i];
                difference_norm_squared += difference[i] * difference[i];
                residual_norm_squared += residual[i] * residual[i];
            }
            const double difference_norm = std::sqrt(difference_norm_squared);
            const double residual_norm = std::sqrt(residual_norm_squared);
            if (difference_norm > std::numeric_limits<double>::epsilon())
            {
                std::vector<double> normalized_difference(difference.size(), 0.0);
                std::vector<double> update_difference(difference.size(), 0.0);
                for (std::size_t i = 0; i < difference.size(); ++i)
                {
                    normalized_difference[i] = difference[i] / difference_norm;
                    update_difference[i] = (charge_step[i] + parameters_.mixing_parameter * difference[i])
                                           / difference_norm;
                }
                normalized_residual_differences_.push_back(std::move(normalized_difference));
                update_differences_.push_back(std::move(update_difference));
                const double weight = residual_norm > parameters_.weight_factor / parameters_.maximal_weight
                                          ? parameters_.weight_factor / residual_norm
                                          : parameters_.maximal_weight;
                broyden_weights_.push_back(std::max(parameters_.minimal_weight,
                                                    std::min(parameters_.maximal_weight, weight)));
                while (broyden_weights_.size() > static_cast<std::size_t>(parameters_.history))
                {
                    broyden_weights_.erase(broyden_weights_.begin());
                    normalized_residual_differences_.erase(normalized_residual_differences_.begin());
                    update_differences_.erase(update_differences_.begin());
                }
            }

            const std::size_t count = broyden_weights_.size();
            acceleration_history_ready = count > 0;
            if (count > 0)
            {
                std::vector<std::vector<double>> matrix(count, std::vector<double>(count, 0.0));
                std::vector<double> rhs(count, 0.0);
                for (std::size_t i = 0; i < count; ++i)
                {
                    rhs[i] = broyden_weights_[i]
                             * vector_dot(normalized_residual_differences_[i], residual);
                    for (std::size_t j = 0; j < count; ++j)
                        matrix[i][j] = broyden_weights_[i] * broyden_weights_[j]
                                       * vector_dot(normalized_residual_differences_[i],
                                                    normalized_residual_differences_[j]);
                    matrix[i][i] += parameters_.inverse_jacobi_weight * parameters_.inverse_jacobi_weight;
                }
                std::vector<double> coefficients;
                if (solve_dense_system(std::move(matrix), std::move(rhs), &coefficients))
                {
                    next = linear_mix(charges, residual, parameters_.mixing_parameter);
                    for (std::size_t history = 0; history < count; ++history)
                    {
                        const double scale = broyden_weights_[history] * coefficients[history];
                        for (std::size_t atom = 0; atom < charges.size(); ++atom)
                            next[atom] -= scale * update_differences_[history][atom];
                    }
                    accelerated = std::all_of(next.begin(), next.end(), [](double value) { return std::isfinite(value); });
                }
            }
        }
        previous_charges_ = charges;
        previous_residual_ = residual;
        has_previous_ = true;
    }

    if (!accelerated)
    {
        next = linear_mix(charges, residual, parameters_.mixing_parameter);
        last_step_ = parameters_.method == "linear" ? "linear"
                     : (acceleration_history_ready ? "linear fallback" : "linear startup");
    }
    else
    {
        last_step_ = parameters_.method == "pulay" ? "Pulay" : "Broyden";
    }
    for (const double value : next)
        if (!std::isfinite(value)) throw std::runtime_error("Native DFTB charge mixer generated a non-finite value");
    restore_total_charge(&next, target_total_charge);
    return next;
}

} // namespace ModuleDFTB
