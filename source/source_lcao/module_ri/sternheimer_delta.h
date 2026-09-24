#ifndef STERNHEIMER_DELTA_H
#define STERNHEIMER_DELTA_H

#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"
#include "source_lcao/module_ri/sternheimer_kq.h"
#include "source_lcao/module_ri/sternheimer_rpa.h"

#include <array>
#include <cstddef>
#include <complex>
#include <string>
#include <vector>

namespace ModuleRI
{

struct SternheimerDeltaVirtualState
{
    SternheimerFDHamiltonian::Vector orbital;
    SternheimerFDHamiltonian::Vector residual;
    double eigenvalue = 0.0;
};

struct SternheimerDeltaPostprocessInput
{
    std::vector<SternheimerFDHamiltonian::Vector> occupied_wavefunctions;
    std::vector<SternheimerDeltaVirtualState> virtual_states;
    std::vector<SternheimerFDHamiltonian::Complex> perturbation_matrix_elements;
    double occupied_eigenvalue = 0.0;
    double omega = 0.0;
    double volume_element = 0.0;
};

struct SternheimerDeltaPostprocessResult
{
    SternheimerFDHamiltonian::Vector out_wavefunction;
    std::vector<SternheimerFDHamiltonian::Complex> sos_coefficients;
    std::vector<SternheimerFDHamiltonian::Complex> pulay_coefficients;
    std::vector<SternheimerFDHamiltonian::Complex> coefficients;
    SternheimerFDHamiltonian::Vector in_sos_wavefunction;
    SternheimerFDHamiltonian::Vector in_pulay_wavefunction;
    SternheimerFDHamiltonian::Vector reconstructed_wavefunction;
    double out_norm = 0.0;
    double reconstruction_error = 0.0;
};

struct SternheimerDeltaCoefficientComponents
{
    std::vector<SternheimerFDHamiltonian::Complex> sos;
    std::vector<SternheimerFDHamiltonian::Complex> pulay;
    std::vector<SternheimerFDHamiltonian::Complex> total;
};

struct SternheimerDeltaPulayOperatorComponents
{
    SternheimerFDHamiltonian::Vector kinetic;
    SternheimerFDHamiltonian::Vector fixed_local;
    SternheimerFDHamiltonian::Vector hxc_local;
    SternheimerFDHamiltonian::Vector nonlocal;
    SternheimerFDHamiltonian::Vector eigenvalue;
    SternheimerFDHamiltonian::Vector total;
};

struct SternheimerDeltaSubspaceOptions
{
    int max_virtual_states = 0;
    double norm_tolerance = 1.0e-10;
    bool retain_grid_functions = true;
    bool evaluate_full_grid_difference = true;
    // Complete candidate spaces can use a generalized Hermitian solve and
    // blocked grid contractions instead of pivoted vector-by-vector MGS.
    bool use_block_generalized_eigensolver = true;
};

enum class SternheimerDeltaABlockMode
{
    ReferenceValueGradient,
    FullGrid
};

SternheimerDeltaABlockMode parse_sternheimer_delta_a_block_mode(const std::string& name);
const char* sternheimer_delta_a_block_mode_name(SternheimerDeltaABlockMode mode);

int sternheimer_delta_virtual_state_limit(int requested_states,
                                          int candidate_states,
                                          int occupied_states);

struct SternheimerDeltaGridFunction
{
    SternheimerFDHamiltonian::Vector values;
    std::array<SternheimerFDHamiltonian::Vector, 3> gradients;
};

struct SternheimerDeltaGridMatrices
{
    // LAPACK column-major storage: element (row, column) is row + row_count * column.
    std::vector<SternheimerFDHamiltonian::Complex> overlap;
    std::vector<SternheimerFDHamiltonian::Complex> kinetic;
    std::vector<SternheimerFDHamiltonian::Complex> local_potential;
    std::vector<SternheimerFDHamiltonian::Complex> nonlocal;
    std::vector<SternheimerFDHamiltonian::Complex> hamiltonian;
    std::size_t row_count = 0;
    std::size_t column_count = 0;
};

// Add one interleaved AO image sample to grid functions. Function values and
// all analytic gradients receive the same ABACUS Bloch phase.
void accumulate_delta_sternheimer_bloch_samples(
    const std::vector<double>& sampled_values,
    const std::array<std::vector<double>, 3>& sampled_gradients,
    int sample_count,
    int orbital_count,
    std::size_t grid_begin,
    std::size_t function_begin,
    const SternheimerReducedKPoint& kpoint,
    const std::array<int, 3>& lattice_translation,
    std::vector<SternheimerDeltaGridFunction>& functions);

// Contract sampled AO chunks directly into selected LCAO states. This avoids
// retaining one full value-and-gradient grid function per AO.
void accumulate_delta_sternheimer_lcao_state_samples(
    const std::vector<double>& sampled_values,
    const std::array<std::vector<double>, 3>& sampled_gradients,
    int sample_count,
    int orbital_count,
    std::size_t grid_begin,
    std::size_t ao_begin,
    const std::vector<std::vector<SternheimerFDHamiltonian::Complex>>& coefficients,
    const SternheimerReducedKPoint& kpoint,
    const std::array<int, 3>& lattice_translation,
    std::vector<SternheimerDeltaGridFunction>& functions);

SternheimerDeltaGridFunction make_delta_sternheimer_grid_function_with_fd_gradients(
    const SternheimerFDHamiltonian::Vector& values,
    const SternheimerFDHamiltonian::Grid& grid);

// Form one LCAO KS state on the uniform grid. The same AO coefficients are
// applied to function values and every analytic-gradient component.
SternheimerDeltaGridFunction linear_combination_delta_sternheimer_grid_functions(
    const std::vector<SternheimerDeltaGridFunction>& basis_functions,
    const std::vector<SternheimerFDHamiltonian::Complex>& coefficients);

// Build an orthonormal grid-metric basis spanning the same functions. Values
// and analytic gradients are transformed with identical MGS coefficients.
std::vector<SternheimerDeltaGridFunction> orthonormalize_delta_sternheimer_grid_functions(
    const std::vector<SternheimerDeltaGridFunction>& functions,
    double volume_element,
    double norm_tolerance);

// Enumerate Gamma-point periodic copies whose cutoff spheres intersect the
// primary orthogonal cell. Nonperiodic grids return only the zero image.
std::vector<std::array<int, 3>> enumerate_delta_sternheimer_periodic_images(
    const SternheimerFDHamiltonian::Grid& grid,
    const std::array<double, 3>& atom_position,
    double cutoff_radius);

struct SternheimerDeltaSubspace
{
    std::vector<SternheimerDeltaVirtualState> virtual_states;
    std::vector<SternheimerDeltaGridFunction> grid_functions;
    int accepted_candidates = 0;
    int discarded_candidates = 0;
    bool used_block_generalized_eigensolver = false;
    double full_grid_hamiltonian_relative_difference = 0.0;
    double full_grid_hamiltonian_max_abs_difference = 0.0;
};

struct SternheimerDeltaLinearResponse
{
    SternheimerDeltaPostprocessResult response;
    SternheimerRPA::SolverResult solver;
    double residual_norm = 0.0;
};

struct SternheimerDeltaFixedSubspace
{
    std::vector<SternheimerFDHamiltonian::Vector> functions;
};

SternheimerDeltaPostprocessResult postprocess_delta_sternheimer_solution(
    const SternheimerFDHamiltonian::Vector& standard_delta_wavefunction,
    const SternheimerDeltaPostprocessInput& input);

// Split the Pulay first-order wavefunction into the operator terms that form
// Q_out (H_grid - epsilon_a) eta_a. This is a diagnostic decomposition; the
// sum of the five operator terms must reproduce total.
SternheimerDeltaPulayOperatorComponents decompose_delta_sternheimer_pulay_operator_terms(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<double>& fixed_local_potential,
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const SternheimerFDHamiltonian::Vector& out_wavefunction,
    double occupied_eigenvalue,
    double omega,
    double volume_element);

// Matrices use LAPACK column-major storage. This solves
// (H_delta - shift S_delta) c = rhs_delta - h_out by two shared-LU right-hand sides.
SternheimerDeltaCoefficientComponents solve_delta_sternheimer_subspace_coefficients(
    const std::vector<SternheimerFDHamiltonian::Complex>& hamiltonian_matrix,
    const std::vector<SternheimerFDHamiltonian::Complex>& overlap_matrix,
    const std::vector<SternheimerFDHamiltonian::Complex>& rhs,
    const std::vector<SternheimerFDHamiltonian::Complex>& hamiltonian_out_coupling,
    SternheimerFDHamiltonian::Complex shift);

// Reject ill-conditioned sampled Gram matrices and states outside the AO span.
std::vector<std::vector<SternheimerFDHamiltonian::Complex>> recover_delta_grid_ao_coefficients(
    const std::vector<SternheimerDeltaGridFunction>& basis,
    const std::vector<const SternheimerFDHamiltonian::Vector*>& states,
    double volume_element,
    double reconstruction_tolerance,
    double& maximum_relative_error,
    double* gram_reciprocal_condition = nullptr);

SternheimerDeltaSubspace build_delta_sternheimer_subspace(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    std::vector<SternheimerFDHamiltonian::Vector> candidate_orbitals,
    double volume_element,
    const SternheimerDeltaSubspaceOptions& options = SternheimerDeltaSubspaceOptions());

// Assemble the reference-code Delta matrices from function values and analytic
// gradients. The finite-difference kinetic operator is deliberately not used.
SternheimerDeltaGridMatrices assemble_delta_sternheimer_grid_matrices(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerDeltaGridFunction>& basis_functions,
    double volume_element);

// Same reference integrals, using the production subspace's blocked BLAS path.
SternheimerDeltaGridMatrices assemble_delta_sternheimer_grid_matrices_fast(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerDeltaGridFunction>& basis_functions,
    double volume_element);

// Rectangular weak forms <left_i|O|right_j> on the same fine grid/Hamiltonian.
// Gradients are caller-supplied analytic Cartesian derivatives of the sampled
// functions, including their Bloch phase, not fractional or periodic-part derivatives.
// T uses hamiltonian.kinetic_prefactor() (1.0 for Ry); no FD kinetic is applied.
// No projection, normalization, or Hermitian cleanup is performed, even for
// equal-sized bases. Either basis may be empty; dimensions are retained.
// Invalid dimensions/nonfinite inputs throw invalid_argument; unrepresentable
// sizes or nonfinite arithmetic results throw overflow_error.
SternheimerDeltaGridMatrices assemble_delta_sternheimer_cross_matrices_fast(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerDeltaGridFunction>& left_functions,
    const std::vector<SternheimerDeltaGridFunction>& right_functions,
    double volume_element);

// Reference-code path: project values and gradients with identical coefficients,
// apply two-pass modified Gram-Schmidt, then diagonalize the independently
// assembled Hermitian Delta Hamiltonian.
SternheimerDeltaSubspace build_reference_delta_sternheimer_subspace(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerDeltaGridFunction>& occupied_functions,
    std::vector<SternheimerDeltaGridFunction> candidate_functions,
    double volume_element,
    const SternheimerDeltaSubspaceOptions& options = SternheimerDeltaSubspaceOptions());

SternheimerDeltaSubspace build_delta_sternheimer_subspace_by_mode(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerDeltaGridFunction>& occupied_functions,
    std::vector<SternheimerDeltaGridFunction> candidate_functions,
    double volume_element,
    const SternheimerDeltaSubspaceOptions& options,
    SternheimerDeltaABlockMode mode);

std::vector<SternheimerFDHamiltonian::Complex> delta_sternheimer_perturbation_matrix_elements(
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<double>& perturbation_potential,
    const SternheimerFDHamiltonian::Vector& occupied_wavefunction,
    double volume_element);

std::vector<SternheimerFDHamiltonian::Complex> delta_sternheimer_perturbation_matrix_elements(
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const SternheimerFDHamiltonian::Vector& perturbation_potential,
    const SternheimerFDHamiltonian::Vector& occupied_wavefunction,
    double volume_element);

SternheimerDeltaFixedSubspace build_delta_sternheimer_fixed_subspace(
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states);

// Assemble the positive-frequency SOS branch from explicit orthonormal
// virtual states. Energies and omega must use the same units.
SternheimerFDHamiltonian::Vector build_delta_sternheimer_sos_wavefunction(
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<SternheimerFDHamiltonian::Complex>& perturbation_matrix_elements,
    double occupied_eigenvalue,
    double omega);

// All overloads use the supplied perturbation elements for the AO-block rhs (-g).
// Only Q*rhs is taken from the grid. The reconstructed residual checks this hybrid equation.
SternheimerDeltaLinearResponse solve_delta_sternheimer_linear_response(
    const SternheimerFDHamiltonian& hamiltonian,
    const SternheimerDeltaFixedSubspace& fixed_subspace,
    double reference_eigenvalue,
    const SternheimerFDHamiltonian::Vector& rhs,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<SternheimerFDHamiltonian::Complex>& perturbation_matrix_elements,
    double omega,
    double volume_element,
    const SternheimerRPA::SolverOptions& options = SternheimerRPA::SolverOptions());

std::vector<SternheimerDeltaLinearResponse> solve_delta_sternheimer_linear_response_batch(
    const SternheimerFDHamiltonian& hamiltonian,
    const SternheimerDeltaFixedSubspace& fixed_subspace,
    double reference_eigenvalue,
    const SternheimerFDHamiltonian::Matrix& rhs,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<std::vector<SternheimerFDHamiltonian::Complex>>& perturbation_matrix_elements,
    double omega,
    double volume_element,
    const SternheimerRPA::SolverOptions& options = SternheimerRPA::SolverOptions());

// Reuse an immutable projector across channel workers. Its basis must outlive the solve.
SternheimerDeltaLinearResponse solve_delta_sternheimer_linear_response(
    const SternheimerFDHamiltonian& hamiltonian,
    const SternheimerSubspaceProjector& fixed_projector,
    double reference_eigenvalue,
    const SternheimerFDHamiltonian::Vector& rhs,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<SternheimerFDHamiltonian::Complex>& perturbation_matrix_elements,
    double omega,
    double volume_element,
    const SternheimerRPA::SolverOptions& options = SternheimerRPA::SolverOptions());

std::vector<SternheimerDeltaLinearResponse> solve_delta_sternheimer_linear_response_batch(
    const SternheimerFDHamiltonian& hamiltonian,
    const SternheimerSubspaceProjector& fixed_projector,
    double reference_eigenvalue,
    const SternheimerFDHamiltonian::Matrix& rhs,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<std::vector<SternheimerFDHamiltonian::Complex>>& perturbation_matrix_elements,
    double omega,
    double volume_element,
    const SternheimerRPA::SolverOptions& options = SternheimerRPA::SolverOptions());

SternheimerDeltaLinearResponse solve_delta_sternheimer_linear_response(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    double reference_eigenvalue,
    const SternheimerFDHamiltonian::Vector& rhs,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<SternheimerFDHamiltonian::Complex>& perturbation_matrix_elements,
    double omega,
    double volume_element,
    const SternheimerRPA::SolverOptions& options = SternheimerRPA::SolverOptions());

SternheimerFDHamiltonian::Complex accumulate_delta_sternheimer_response(
    const SternheimerFDHamiltonian::Vector& probe_wavefunction,
    const SternheimerDeltaPostprocessResult& response,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    double volume_element);

} // namespace ModuleRI

#endif
