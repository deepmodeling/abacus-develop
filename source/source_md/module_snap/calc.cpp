#include "calc.h"

#include <stdexcept>

namespace ModuleMD
{
namespace SNAP
{

Calculator::Calculator(const ModelData& model)
    : model_(model),
      desc_eval_(model.parameters),
      desc_size_(0)
{
    const std::size_t channels = desc_eval_.channel_count();
    desc_size_ = desc_eval_.bispectrum_size() * channels * channels * channels;
    if (model_.coefficients.size() != model_.parameters.elements.size())
    {
        throw std::invalid_argument("SNAP coefficient blocks do not match the element list.");
    }
    const std::size_t expected_count = desc_size_ + 1
                                       + (model_.parameters.quadratic_flag
                                              ? desc_size_ * (desc_size_ + 1) / 2
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
    result.descriptor = desc_eval_.calc_bispec(
        desc_eval_.accumulate_u(neighbors, central_element));
    if (model_.parameters.bzero_flag)
    {
        const std::vector<double> bzero =
            desc_eval_.calc_bispec(
                desc_eval_.accumulate_u(std::vector<Neighbor>(), central_element));
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

ForceResult Calculator::eval_force(const std::vector<Neighbor>& neighbors,
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
    const std::vector<std::complex<double> > u = desc_eval_.accumulate_u(neighbors, central_element);
    for (std::size_t neighbor_index = 0; neighbor_index < neighbors.size(); ++neighbor_index)
    {
        for (int direction = 0; direction < 3; ++direction)
        {
            const UWithDerivative u_derivative =
                desc_eval_.calc_u_deriv(neighbors[neighbor_index], direction);
            std::vector<std::complex<double> > du(u.size(), std::complex<double>(0.0, 0.0));
            const std::size_t deriv_offset = model_.parameters.chem_flag
                                                       ? neighbors[neighbor_index].element_index
                                                             * desc_eval_.u_size()
                                                       : 0;
            const std::vector<std::complex<double> >& one_neigh_deriv = u_derivative.derivatives;
            for (std::size_t index = 0; index < one_neigh_deriv.size(); ++index)
            {
                du[deriv_offset + index] = one_neigh_deriv[index];
            }
            const std::vector<double> db = desc_eval_.calc_b_deriv(u, du);
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

ForceResult Calculator::eval_force_fd(const std::vector<Neighbor>& neighbors,
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
                          const std::vector<std::size_t>& neigh_indices,
                          const std::vector<Neighbor>& neighbors,
                          const ForceResult& result)
{
    if (center_index >= forces_.size() || neigh_indices.size() != neighbors.size()
        || result.neighbor_forces.size() != neighbors.size() || result.center_force.size() != 3)
    {
        throw std::invalid_argument("SNAP cell contribution has inconsistent atom or force indices.");
    }
    for (std::size_t index = 0; index < neigh_indices.size(); ++index)
    {
        if (neigh_indices[index] >= forces_.size() || result.neighbor_forces[index].size() != 3)
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
        const std::size_t atom = neigh_indices[index];
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
