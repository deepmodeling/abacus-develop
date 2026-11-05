#include "snap_calculator.h"

#include <stdexcept>

namespace ModuleMD
{
namespace SNAP
{

Calculator::Calculator(const ModelData& model)
    : model_(model),
      descriptor_evaluator_(model.parameters),
      descriptor_size_(0)
{
    const std::size_t channels = descriptor_evaluator_.channel_count();
    descriptor_size_ = descriptor_evaluator_.bispectrum_size() * channels * channels * channels;
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
            const std::string message = "SNAP coefficient count for element " + model_.parameters.elements[element]
                                        + " does not match the descriptor size: expected "
                                        + std::to_string(expected_count) + ", found "
                                        + std::to_string(model_.coefficients[element].size()) + ".";
            throw std::invalid_argument(message);
        }
    }
    if (model_.parameters.bzero_flag)
    {
        const std::vector<Neighbor> empty_neighbors;
        self_spectra_.reserve(model_.coefficients.size());
        for (std::size_t element = 0; element < model_.coefficients.size(); ++element)
        {
            const std::vector<std::complex<double> > self_u
                = descriptor_evaluator_.accumulate_u(empty_neighbors, element);
            const std::vector<double> self_bispectrum = descriptor_evaluator_.compute_bispectrum(self_u);
            self_spectra_.push_back(self_bispectrum);
        }
    }
}

EnergyResult Calculator::compute(const std::vector<Neighbor>& neighbors, const std::size_t central_element) const
{
    if (central_element >= model_.coefficients.size())
    {
        throw std::out_of_range("SNAP central element index is outside the model.");
    }

    const std::vector<std::complex<double> > u = descriptor_evaluator_.accumulate_u(neighbors, central_element);
    return compute_from_u(u, central_element);
}

EnergyResult Calculator::compute_from_u(const std::vector<std::complex<double> >& u,
                                        const std::size_t central_element) const
{
    EnergyResult result;
    result.descriptor = descriptor_evaluator_.compute_bispectrum(u);
    if (model_.parameters.bzero_flag)
    {
        for (std::size_t index = 0; index < result.descriptor.size(); ++index)
        {
            result.descriptor[index] -= self_spectra_[central_element][index];
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
    const std::vector<std::complex<double> > u = descriptor_evaluator_.accumulate_u(neighbors, central_element);
    result.energy = compute_from_u(u, central_element);
    result.center_force.assign(3, 0.0);
    result.neighbor_forces.assign(neighbors.size(), std::vector<double>(3, 0.0));
    const std::vector<double>& coefficients = model_.coefficients[central_element];
    std::vector<double> force_coefficients(result.energy.descriptor.size(), 0.0);
    for (std::size_t i = 0; i < force_coefficients.size(); ++i)
    {
        force_coefficients[i] = coefficients[i + 1];
    }
    if (model_.parameters.quadratic_flag)
    {
        std::size_t coefficient = force_coefficients.size() + 1;
        for (std::size_t i = 0; i < force_coefficients.size(); ++i)
        {
            for (std::size_t j = i; j < force_coefficients.size(); ++j)
            {
                const double value = coefficients[coefficient++];
                if (i == j)
                {
                    force_coefficients[i] += value * result.energy.descriptor[i];
                }
                else
                {
                    force_coefficients[i] += value * result.energy.descriptor[j];
                    force_coefficients[j] += value * result.energy.descriptor[i];
                }
            }
        }
    }
    for (std::size_t neighbor_index = 0; neighbor_index < neighbors.size(); ++neighbor_index)
    {
        for (int direction = 0; direction < 3; ++direction)
        {
            const UWithDerivative u_derivative =
                descriptor_evaluator_.compute_u_with_derivative(neighbors[neighbor_index], direction);
            std::vector<std::complex<double> > du(u.size(), std::complex<double>(0.0, 0.0));
            const std::size_t derivative_offset = model_.parameters.chem_flag
                                                       ? neighbors[neighbor_index].element_index
                                                             * descriptor_evaluator_.u_size()
                                                       : 0;
            const std::vector<std::complex<double> >& one_neighbor_derivative = u_derivative.derivatives;
            for (std::size_t index = 0; index < one_neighbor_derivative.size(); ++index)
            {
                du[derivative_offset + index] = one_neighbor_derivative[index];
            }
            const std::vector<double> db = descriptor_evaluator_.compute_bispectrum_derivative(u, du);
            double derivative = 0.0;
            for (std::size_t index = 0; index < db.size(); ++index)
            {
                derivative += force_coefficients[index] * db[index];
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

CellAccumulator::CellAccumulator(const std::size_t atom_count) : forces_(atom_count)
{
    for (std::size_t atom = 0; atom < forces_.size(); ++atom)
    {
        forces_[atom].fill(0.0);
    }
    virial_.fill(0.0);
}

void CellAccumulator::add(const std::size_t center_index,
                          const std::vector<std::size_t>& neighbor_indices,
                          const std::vector<Neighbor>& neighbors,
                          const ForceResult& result)
{
    if (center_index >= forces_.size() || neighbor_indices.size() != neighbors.size()
        || result.neighbor_forces.size() != neighbors.size() || result.center_force.size() != 3)
    {
        throw std::invalid_argument("SNAP cell contribution has inconsistent atom or force indices.");
    }
    for (std::size_t index = 0; index < neighbor_indices.size(); ++index)
    {
        if (neighbor_indices[index] >= forces_.size() || result.neighbor_forces[index].size() != 3)
        {
            throw std::invalid_argument("SNAP neighbor force index is outside the cell.");
        }
    }

    energy_ += result.energy.energy;
    for (int direction = 0; direction < 3; ++direction)
    {
        forces_[center_index][static_cast<std::size_t>(direction)] +=
            result.center_force[static_cast<std::size_t>(direction)];
    }
    for (std::size_t index = 0; index < neighbors.size(); ++index)
    {
        const Neighbor& neighbor = neighbors[index];
        const double displacement[3] = {neighbor.x, neighbor.y, neighbor.z};
        const std::size_t atom = neighbor_indices[index];
        for (int direction = 0; direction < 3; ++direction)
        {
            forces_[atom][static_cast<std::size_t>(direction)] +=
                result.neighbor_forces[index][static_cast<std::size_t>(direction)];
        }
        for (int row = 0; row < 3; ++row)
        {
            for (int column = 0; column < 3; ++column)
            {
                virial_[static_cast<std::size_t>(3 * row + column)] +=
                    displacement[row] * result.neighbor_forces[index][static_cast<std::size_t>(column)];
            }
        }
    }
}

double CellAccumulator::energy() const
{
    return energy_;
}

const std::vector<std::array<double, 3> >& CellAccumulator::forces() const
{
    return forces_;
}

const std::array<double, 9>& CellAccumulator::virial() const
{
    return virial_;
}

} // namespace SNAP
} // namespace ModuleMD
