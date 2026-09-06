#include "source_lcao/module_ri/sternheimer_periodic_solver.h"

#include <cmath>
#include <memory>
#include <limits>
#include <stdexcept>
#include <utility>

namespace ModuleRI
{

namespace
{
std::vector<const SternheimerFDHamiltonian::Vector*> fixed_subspace_views(
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states)
{
    std::vector<const SternheimerFDHamiltonian::Vector*> views;
    views.reserve(occupied.size() + virtual_states.size());
    // Match the existing fixed-subspace projection order.
    for (const auto& state: occupied)
    {
        views.push_back(&state);
    }
    for (const auto& state: virtual_states)
    {
        views.push_back(&state.orbital);
    }
    return views;
}
} // namespace

SternheimerPeriodicResponseProjectors::SternheimerPeriodicResponseProjectors(
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const double volume_element,
    const bool pack_batch)
    : occupied(occupied_wavefunctions, volume_element, pack_batch),
      fixed(fixed_subspace_views(occupied_wavefunctions, virtual_states), volume_element, pack_batch)
{
}

std::size_t SternheimerPeriodicResponseProjectors::estimated_storage_bytes(
    const std::size_t grid_size, const std::size_t occupied_count, const std::size_t virtual_count,
    const bool pack_batch)
{
    const auto maximum = std::numeric_limits<std::size_t>::max();
    const auto value_bytes = sizeof(SternheimerFDHamiltonian::Complex);
    const auto metadata_bytes = value_bytes + sizeof(const SternheimerFDHamiltonian::Vector*);
    const auto packed_grid_size = pack_batch ? grid_size : 0;
    if (occupied_count > (maximum - virtual_count) / 2
        || packed_grid_size > (maximum - metadata_bytes) / value_bytes)
    {
        throw std::overflow_error("Periodic Sternheimer shared-projector storage size overflow.");
    }
    const auto count = 2 * occupied_count + virtual_count;
    const auto bytes_per_state = packed_grid_size * value_bytes + metadata_bytes;
    if (count > maximum / bytes_per_state)
    {
        throw std::overflow_error("Periodic Sternheimer shared-projector storage size overflow.");
    }
    return count * bytes_per_state;
}

SternheimerPeriodicLinearResponse solve_sternheimer_periodic_linear_response(
    const bool use_delta_sternheimer,
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    const double reference_eigenvalue,
    const SternheimerFDHamiltonian::Vector& rhs,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<SternheimerFDHamiltonian::Complex>& perturbation_matrix_elements,
    const double omega,
    const double volume_element,
    const SternheimerRPA::SolverOptions& options,
    const SternheimerPeriodicResponseProjectors* shared_projectors)
{
    SternheimerPeriodicLinearResponse result;
    result.projected_rhs = rhs;
    const auto dot = [volume_element](const SternheimerFDHamiltonian::Vector& lhs,
                                      const SternheimerFDHamiltonian::Vector& rhs_vector) {
        return sternheimer_fd_grid_dot(lhs, rhs_vector, volume_element);
    };
    SternheimerRPA::project_out_subspace(occupied_wavefunctions, dot, result.projected_rhs);
    if (use_delta_sternheimer)
    {
        std::unique_ptr<SternheimerPeriodicResponseProjectors> local_projectors;
        if (shared_projectors == nullptr)
        {
            local_projectors = std::make_unique<SternheimerPeriodicResponseProjectors>(
                occupied_wavefunctions, virtual_states, volume_element, false);
            shared_projectors = local_projectors.get();
        }
        SternheimerDeltaLinearResponse response
            = solve_delta_sternheimer_linear_response(hamiltonian,
                                                       shared_projectors->fixed,
                                                       reference_eigenvalue,
                                                       rhs,
                                                       virtual_states,
                                                       perturbation_matrix_elements,
                                                       omega,
                                                       volume_element,
                                                       options);
        result.wavefunction = response.response.reconstructed_wavefunction;
        result.solver = response.solver;
        result.residual_norm = response.residual_norm;
        result.full_grid_equation_residual_norm
            = sternheimer_fd_linear_response_residual_norm(hamiltonian,
                                                            occupied_wavefunctions,
                                                            reference_eigenvalue,
                                                            rhs,
                                                            result.wavefunction,
                                                            omega,
                                                            volume_element);
        result.has_delta_components = true;
        result.delta_components = std::move(response.response);
        return result;
    }

    const SternheimerFDLinearResponse response
        = solve_sternheimer_fd_linear_response(hamiltonian,
                                                occupied_wavefunctions,
                                                reference_eigenvalue,
                                                rhs,
                                                omega,
                                                volume_element,
                                                options);
    result.wavefunction = response.delta_wavefunction;
    result.solver = response.solver;
    result.residual_norm = response.residual_norm;
    result.full_grid_equation_residual_norm = response.residual_norm;
    return result;
}

std::vector<SternheimerPeriodicLinearResponse> solve_sternheimer_periodic_linear_response_batch(
    const bool use_delta_sternheimer,
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    const double reference_eigenvalue,
    const SternheimerFDHamiltonian::Matrix& rhs,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<std::vector<SternheimerFDHamiltonian::Complex>>& perturbation_matrix_elements,
    const double omega,
    const double volume_element,
    const SternheimerRPA::SolverOptions& options,
    const SternheimerPeriodicResponseProjectors* shared_projectors)
{
    if (rhs.size() != perturbation_matrix_elements.size())
    {
        throw std::invalid_argument(
            "Periodic Sternheimer batch response requires one perturbation-element vector per rhs.");
    }
    std::vector<SternheimerPeriodicLinearResponse> results(rhs.size());
    if (!use_delta_sternheimer)
    {
        for (std::size_t column = 0; column != rhs.size(); ++column)
        {
            results[column] = solve_sternheimer_periodic_linear_response(false,
                                                                          hamiltonian,
                                                                          occupied_wavefunctions,
                                                                          reference_eigenvalue,
                                                                          rhs[column],
                                                                          virtual_states,
                                                                          perturbation_matrix_elements[column],
                                                                          omega,
                                                                          volume_element,
                                                                          options);
        }
        return results;
    }

    std::unique_ptr<SternheimerPeriodicResponseProjectors> local_projectors;
    if (shared_projectors == nullptr)
    {
        local_projectors = std::make_unique<SternheimerPeriodicResponseProjectors>(
            occupied_wavefunctions, virtual_states, volume_element);
        shared_projectors = local_projectors.get();
    }
    const auto& occupied_projector = shared_projectors->occupied;
    SternheimerFDHamiltonian::Matrix projected_rhs = rhs;
    occupied_projector.project_batch(projected_rhs);
    for (std::size_t column = 0; column != rhs.size(); ++column)
    {
        results[column].projected_rhs = projected_rhs[column];
    }

    std::vector<SternheimerDeltaLinearResponse> responses
        = solve_delta_sternheimer_linear_response_batch(hamiltonian,
                                                         shared_projectors->fixed,
                                                         reference_eigenvalue,
                                                         rhs,
                                                         virtual_states,
                                                         perturbation_matrix_elements,
                                                         omega,
                                                         volume_element,
                                                         options);
    SternheimerFDHamiltonian::Matrix wavefunctions(rhs.size());
    for (std::size_t column = 0; column != rhs.size(); ++column)
    {
        results[column].wavefunction = responses[column].response.reconstructed_wavefunction;
        wavefunctions[column] = results[column].wavefunction;
        results[column].solver = responses[column].solver;
        results[column].residual_norm = responses[column].residual_norm;
        results[column].has_delta_components = true;
        results[column].delta_components = std::move(responses[column].response);
    }

    SternheimerFDHamiltonian::Matrix projected_wavefunctions = wavefunctions;
    occupied_projector.project_batch(projected_wavefunctions);
    SternheimerFDHamiltonian::Matrix residuals;
    hamiltonian.apply_batch(projected_wavefunctions, residuals);
    const SternheimerFDHamiltonian::Complex shift(-reference_eigenvalue, omega);
#pragma omp parallel for schedule(static)
    for (std::size_t column = 0; column != residuals.size(); ++column)
    {
        for (std::size_t ir = 0; ir != residuals[column].size(); ++ir)
        {
            residuals[column][ir] += shift * projected_wavefunctions[column][ir];
        }
    }
    occupied_projector.project_batch(residuals);
#pragma omp parallel for schedule(static)
    for (std::size_t column = 0; column != residuals.size(); ++column)
    {
        double norm_squared = 0.0;
        for (std::size_t ir = 0; ir != residuals[column].size(); ++ir)
        {
            residuals[column][ir] -= projected_rhs[column][ir];
            norm_squared += std::norm(residuals[column][ir]);
        }
        results[column].full_grid_equation_residual_norm
            = std::sqrt(volume_element * norm_squared);
    }
    return results;
}

} // namespace ModuleRI
