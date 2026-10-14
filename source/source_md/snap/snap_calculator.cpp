#include "snap_calculator.h"

#include <stdexcept>

namespace ModuleMD
{
namespace SNAP
{

Calculator::Calculator(const ModelData& model)
    : model_(model), descriptor_evaluator_(model.parameters), descriptor_size_(descriptor_evaluator_.bispectrum_size())
{
    if (model_.coefficients.size() != model_.parameters.elements.size())
    {
        throw std::invalid_argument("SNAP coefficient blocks do not match the element list.");
    }
    const std::size_t expected_count = descriptor_size_ + 1
                                       + (model_.parameters.quadratic_flag
                                              ? descriptor_size_ * (descriptor_size_ + 1) / 2
                                              : 0);
    for (std::size_t element = 0; element < model_.coefficients.size(); ++element)
    {
        if (model_.coefficients[element].size() != expected_count)
        {
            throw std::invalid_argument("SNAP coefficient count does not match the descriptor size.");
        }
    }
}

EnergyResult Calculator::compute(const std::vector<Neighbor>& neighbors, const std::size_t central_element) const
{
    if (central_element >= model_.coefficients.size())
    {
        throw std::out_of_range("SNAP central element index is outside the model.");
    }

    EnergyResult result;
    result.descriptor = descriptor_evaluator_.compute_bispectrum(descriptor_evaluator_.accumulate_u(neighbors));
    if (model_.parameters.bzero_flag)
    {
        const std::vector<double> bzero =
            descriptor_evaluator_.compute_bispectrum(descriptor_evaluator_.accumulate_u(std::vector<Neighbor>()));
        for (std::size_t index = 0; index < result.descriptor.size(); ++index)
        {
            result.descriptor[index] -= bzero[index];
        }
    }
    const std::vector<double>& coefficients = model_.coefficients[central_element];
    result.energy = coefficients[0];
    for (std::size_t index = 0; index < result.descriptor.size(); ++index)
    {
        result.energy += coefficients[index + 1] * result.descriptor[index];
    }
    if (model_.parameters.quadratic_flag)
    {
        std::size_t coefficient = result.descriptor.size() + 1;
        for (std::size_t i = 0; i < result.descriptor.size(); ++i)
        {
            for (std::size_t j = i; j < result.descriptor.size(); ++j)
            {
                result.energy += coefficients[coefficient++] * result.descriptor[i] * result.descriptor[j]
                                 * (i == j ? 0.5 : 1.0);
            }
        }
    }
    return result;
}

ForceResult Calculator::compute_with_analytic_force(const std::vector<Neighbor>& neighbors,
                                                    const std::size_t central_element) const
{
    if (central_element >= model_.coefficients.size())
    {
        throw std::out_of_range("SNAP central element index is outside the model.");
    }
    ForceResult result;
    result.energy = compute(neighbors, central_element);
    result.center_force.assign(3, 0.0);
    result.neighbor_forces.assign(neighbors.size(), std::vector<double>(3, 0.0));
    const std::vector<double>& coefficients = model_.coefficients[central_element];
    for (std::size_t neighbor_index = 0; neighbor_index < neighbors.size(); ++neighbor_index)
    {
        for (int direction = 0; direction < 3; ++direction)
        {
            const UWithDerivative u_derivative =
                descriptor_evaluator_.compute_u_with_derivative(neighbors[neighbor_index], direction);
            std::vector<std::complex<double> > u = descriptor_evaluator_.accumulate_u(neighbors);
            std::vector<std::complex<double> > du(u.size(), std::complex<double>(0.0, 0.0));
            for (std::size_t index = 0; index < du.size(); ++index)
            {
                du[index] = u_derivative.derivatives[index];
            }
            const std::vector<double> db = descriptor_evaluator_.compute_bispectrum_derivative(u, du);
            double derivative = 0.0;
            for (std::size_t index = 0; index < db.size(); ++index)
            {
                derivative += coefficients[index + 1] * db[index];
            }
            if (model_.parameters.quadratic_flag)
            {
                std::size_t coefficient = db.size() + 1;
                for (std::size_t i = 0; i < db.size(); ++i)
                {
                    for (std::size_t j = i; j < db.size(); ++j)
                    {
                        const double factor = i == j ? 0.5 : 1.0;
                        derivative += coefficients[coefficient++] * factor
                                      * (db[i] * result.energy.descriptor[j]
                                         + result.energy.descriptor[i] * db[j]);
                    }
                }
            }
            result.center_force[static_cast<std::size_t>(direction)] += derivative;
            result.neighbor_forces[neighbor_index][static_cast<std::size_t>(direction)] = -derivative;
        }
    }
    return result;
}

ForceResult Calculator::compute_with_finite_difference_force(const std::vector<Neighbor>& neighbors,
                                                             const std::size_t central_element,
                                                             const double displacement) const
{
    if (displacement <= 0.0)
    {
        throw std::invalid_argument("SNAP finite-difference displacement must be positive.");
    }
    ForceResult result;
    result.energy = compute(neighbors, central_element);
    result.center_force.assign(3, 0.0);
    result.neighbor_forces.assign(neighbors.size(), std::vector<double>(3, 0.0));

    for (std::size_t neighbor_index = 0; neighbor_index < neighbors.size(); ++neighbor_index)
    {
        for (int direction = 0; direction < 3; ++direction)
        {
            std::vector<Neighbor> plus = neighbors;
            std::vector<Neighbor> minus = neighbors;
            if (direction == 0)
            {
                plus[neighbor_index].x += displacement;
                minus[neighbor_index].x -= displacement;
            }
            else if (direction == 1)
            {
                plus[neighbor_index].y += displacement;
                minus[neighbor_index].y -= displacement;
            }
            else
            {
                plus[neighbor_index].z += displacement;
                minus[neighbor_index].z -= displacement;
            }
            const double derivative =
                (compute(plus, central_element).energy - compute(minus, central_element).energy) / (2.0 * displacement);
            result.center_force[static_cast<std::size_t>(direction)] += derivative;
            result.neighbor_forces[neighbor_index][static_cast<std::size_t>(direction)] -= derivative;
        }
    }
    return result;
}

} // namespace SNAP
} // namespace ModuleMD
