#ifndef STERNHEIMER_PERIODIC_SOLVER_H
#define STERNHEIMER_PERIODIC_SOLVER_H

#include "source_lcao/module_ri/sternheimer_delta.h"
#include "source_lcao/module_ri/sternheimer_fd_solver.h"

#include <stdexcept>
#include <string>

namespace ModuleRI
{

class SternheimerExcitationError : public std::runtime_error
{
  public:
    using std::runtime_error::runtime_error;
};

// Validate one periodic k/k+q pair using all source occupied energies in Ry.
// Return the positive minimum excitation, or throw without modifying inputs.
// The caller supplies grid/k-pair context and handles fatal process termination;
// workers without an owned pair must not call this with empty placeholder data.
double require_positive_periodic_delta_excitations(
    const std::vector<SternheimerDeltaVirtualState>& target_states,
    const std::vector<double>& source_occupied_eigenvalues_ry,
    const std::string& diagnostic_context);

struct SternheimerPeriodicLinearResponse
{
    SternheimerFDHamiltonian::Vector projected_rhs;
    SternheimerFDHamiltonian::Vector wavefunction;
    SternheimerRPA::SolverResult solver;
    double residual_norm = 0.0;
    double full_grid_equation_residual_norm = 0.0;
    bool has_delta_components = false;
    SternheimerDeltaPostprocessResult delta_components;
};

// Build once per k/k+q pair, before worker admission. Basis vectors must stay
// unchanged and alive until every channel worker has finished.
struct SternheimerPeriodicResponseProjectors
{
    SternheimerPeriodicResponseProjectors(
        const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
        const std::vector<SternheimerDeltaVirtualState>& virtual_states,
        double volume_element,
        bool pack_batch = true);

    static std::size_t estimated_storage_bytes(std::size_t grid_size,
                                               std::size_t occupied_count,
                                               std::size_t virtual_count,
                                               bool pack_batch = true);

    SternheimerSubspaceProjector occupied;
    SternheimerSubspaceProjector fixed;
};

SternheimerPeriodicLinearResponse solve_sternheimer_periodic_linear_response(
    bool use_delta_sternheimer,
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    double reference_eigenvalue,
    const SternheimerFDHamiltonian::Vector& rhs,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<SternheimerFDHamiltonian::Complex>& perturbation_matrix_elements,
    double omega,
    double volume_element,
    const SternheimerRPA::SolverOptions& options = SternheimerRPA::SolverOptions(),
    const SternheimerPeriodicResponseProjectors* shared_projectors = nullptr);

std::vector<SternheimerPeriodicLinearResponse> solve_sternheimer_periodic_linear_response_batch(
    bool use_delta_sternheimer,
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    double reference_eigenvalue,
    const SternheimerFDHamiltonian::Matrix& rhs,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<std::vector<SternheimerFDHamiltonian::Complex>>& perturbation_matrix_elements,
    double omega,
    double volume_element,
    const SternheimerRPA::SolverOptions& options = SternheimerRPA::SolverOptions(),
    const SternheimerPeriodicResponseProjectors* shared_projectors = nullptr);

} // namespace ModuleRI

#endif
