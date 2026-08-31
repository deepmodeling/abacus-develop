#include "source_lcao/module_ri/sternheimer_periodic_solver.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace ModuleRI
{

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
    const SternheimerRPA::SolverOptions& options)
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
        SternheimerDeltaLinearResponse response
            = solve_delta_sternheimer_linear_response(hamiltonian,
                                                       occupied_wavefunctions,
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
    const SternheimerRPA::SolverOptions& options)
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

    const SternheimerSubspaceProjector occupied_projector(occupied_wavefunctions, volume_element, true);
    SternheimerFDHamiltonian::Matrix projected_rhs = rhs;
    occupied_projector.project_batch(projected_rhs);
    for (std::size_t column = 0; column != rhs.size(); ++column)
    {
        results[column].projected_rhs = projected_rhs[column];
    }

    const SternheimerDeltaFixedSubspace fixed_subspace
        = build_delta_sternheimer_fixed_subspace(occupied_wavefunctions, virtual_states);
    std::vector<SternheimerDeltaLinearResponse> responses
        = solve_delta_sternheimer_linear_response_batch(hamiltonian,
                                                         fixed_subspace,
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
