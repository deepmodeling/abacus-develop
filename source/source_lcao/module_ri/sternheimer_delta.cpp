#include "source_lcao/module_ri/sternheimer_delta.h"

#include "source_lcao/module_ri/sternheimer_fd_solver.h"
#include "source_lcao/module_ri/sternheimer_rpa.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace ModuleRI
{
namespace
{

using Complex = SternheimerFDHamiltonian::Complex;
using Vector = SternheimerFDHamiltonian::Vector;

void check_vector_size(const Vector& vector, const std::size_t expected_size, const std::string& context)
{
    if (vector.size() != expected_size)
    {
        throw std::invalid_argument(context + " vector size does not match the Sternheimer delta grid.");
    }
}

void axpy(const Complex alpha, const Vector& x, Vector& y)
{
    check_vector_size(x, y.size(), "Sternheimer delta axpy");
    for (std::size_t ir = 0; ir != y.size(); ++ir)
    {
        y[ir] += alpha * x[ir];
    }
}

double difference_norm(const Vector& lhs, const Vector& rhs, const double volume_element)
{
    check_vector_size(lhs, rhs.size(), "Sternheimer delta difference_norm");
    Vector difference(lhs.size(), Complex(0.0, 0.0));
    for (std::size_t ir = 0; ir != lhs.size(); ++ir)
    {
        difference[ir] = lhs[ir] - rhs[ir];
    }
    return sternheimer_fd_grid_norm(difference, volume_element);
}

std::vector<Vector> collect_virtual_orbitals(const std::vector<SternheimerDeltaVirtualState>& virtual_states)
{
    std::vector<Vector> orbitals;
    orbitals.reserve(virtual_states.size());
    for (const SternheimerDeltaVirtualState& state: virtual_states)
    {
        orbitals.push_back(state.orbital);
    }
    return orbitals;
}

void validate_postprocess_input(const Vector& standard_delta_wavefunction,
                                const SternheimerDeltaPostprocessInput& input)
{
    if (standard_delta_wavefunction.empty())
    {
        throw std::invalid_argument("Sternheimer delta postprocess requires a non-empty standard solution.");
    }
    if (input.volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta postprocess requires a positive grid volume element.");
    }
    if (input.virtual_states.size() != input.perturbation_matrix_elements.size())
    {
        throw std::invalid_argument(
            "Sternheimer delta postprocess requires one perturbation matrix element per virtual state.");
    }

    const std::size_t grid_size = standard_delta_wavefunction.size();
    for (const Vector& occupied: input.occupied_wavefunctions)
    {
        check_vector_size(occupied, grid_size, "Sternheimer delta occupied state");
    }
    for (const SternheimerDeltaVirtualState& state: input.virtual_states)
    {
        check_vector_size(state.orbital, grid_size, "Sternheimer delta virtual orbital");
        check_vector_size(state.residual, grid_size, "Sternheimer delta virtual residual");
    }
}

} // namespace

SternheimerDeltaPostprocessResult postprocess_delta_sternheimer_solution(
    const SternheimerFDHamiltonian::Vector& standard_delta_wavefunction,
    const SternheimerDeltaPostprocessInput& input)
{
    validate_postprocess_input(standard_delta_wavefunction, input);

    auto dot = [&input](const Vector& lhs, const Vector& rhs) {
        return sternheimer_fd_grid_dot(lhs, rhs, input.volume_element);
    };

    SternheimerDeltaPostprocessResult result;
    result.out_wavefunction = standard_delta_wavefunction;
    SternheimerRPA::project_out_subspace(input.occupied_wavefunctions, dot, result.out_wavefunction);
    const std::vector<Vector> virtual_orbitals = collect_virtual_orbitals(input.virtual_states);
    SternheimerRPA::project_out_subspace(virtual_orbitals, dot, result.out_wavefunction);

    result.coefficients.assign(input.virtual_states.size(), Complex(0.0, 0.0));
    for (std::size_t ia = 0; ia != input.virtual_states.size(); ++ia)
    {
        const SternheimerDeltaVirtualState& state = input.virtual_states[ia];
        const Complex denominator(input.occupied_eigenvalue - state.eigenvalue, -input.omega);
        if (std::abs(denominator) < 1.0e-30)
        {
            throw std::runtime_error("Sternheimer delta postprocess found a singular virtual-state denominator.");
        }
        const Complex residual_overlap = sternheimer_fd_grid_dot(state.residual,
                                                                 result.out_wavefunction,
                                                                 input.volume_element);
        result.coefficients[ia] = (input.perturbation_matrix_elements[ia] + residual_overlap) / denominator;
    }

    result.reconstructed_wavefunction = result.out_wavefunction;
    for (std::size_t ia = 0; ia != input.virtual_states.size(); ++ia)
    {
        axpy(result.coefficients[ia], input.virtual_states[ia].orbital, result.reconstructed_wavefunction);
    }
    result.out_norm = sternheimer_fd_grid_norm(result.out_wavefunction, input.volume_element);
    result.reconstruction_error
        = difference_norm(result.reconstructed_wavefunction, standard_delta_wavefunction, input.volume_element);
    return result;
}

SternheimerFDHamiltonian::Complex accumulate_delta_sternheimer_response(
    const SternheimerFDHamiltonian::Vector& probe_wavefunction,
    const SternheimerDeltaPostprocessResult& response,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const double volume_element)
{
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta response requires a positive grid volume element.");
    }
    if (response.coefficients.size() != virtual_states.size())
    {
        throw std::invalid_argument("Sternheimer delta response coefficient count does not match virtual states.");
    }
    check_vector_size(response.out_wavefunction, probe_wavefunction.size(), "Sternheimer delta response out state");

    Complex value = sternheimer_fd_grid_dot(probe_wavefunction, response.out_wavefunction, volume_element);
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        check_vector_size(virtual_states[ia].orbital, probe_wavefunction.size(), "Sternheimer delta response virtual");
        value += response.coefficients[ia]
            * sternheimer_fd_grid_dot(probe_wavefunction, virtual_states[ia].orbital, volume_element);
    }
    return value;
}

} // namespace ModuleRI
