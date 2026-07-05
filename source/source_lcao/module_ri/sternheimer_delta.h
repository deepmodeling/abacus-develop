#ifndef STERNHEIMER_DELTA_H
#define STERNHEIMER_DELTA_H

#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"
#include "source_lcao/module_ri/sternheimer_rpa.h"

#include <complex>
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
    std::vector<SternheimerFDHamiltonian::Complex> coefficients;
    SternheimerFDHamiltonian::Vector reconstructed_wavefunction;
    double out_norm = 0.0;
    double reconstruction_error = 0.0;
};

struct SternheimerDeltaSubspaceOptions
{
    int max_virtual_states = 0;
    double norm_tolerance = 1.0e-10;
};

struct SternheimerDeltaSubspace
{
    std::vector<SternheimerDeltaVirtualState> virtual_states;
    int accepted_candidates = 0;
    int discarded_candidates = 0;
};

struct SternheimerDeltaLinearResponse
{
    SternheimerDeltaPostprocessResult response;
    SternheimerRPA::SolverResult solver;
    double residual_norm = 0.0;
};

SternheimerDeltaPostprocessResult postprocess_delta_sternheimer_solution(
    const SternheimerFDHamiltonian::Vector& standard_delta_wavefunction,
    const SternheimerDeltaPostprocessInput& input);

SternheimerDeltaSubspace build_delta_sternheimer_subspace(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    const std::vector<SternheimerFDHamiltonian::Vector>& candidate_orbitals,
    double volume_element,
    const SternheimerDeltaSubspaceOptions& options = SternheimerDeltaSubspaceOptions());

std::vector<SternheimerFDHamiltonian::Complex> delta_sternheimer_perturbation_matrix_elements(
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<double>& perturbation_potential,
    const SternheimerFDHamiltonian::Vector& occupied_wavefunction,
    double volume_element);

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
