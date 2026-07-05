#ifndef STERNHEIMER_DELTA_H
#define STERNHEIMER_DELTA_H

#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"

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

SternheimerDeltaPostprocessResult postprocess_delta_sternheimer_solution(
    const SternheimerFDHamiltonian::Vector& standard_delta_wavefunction,
    const SternheimerDeltaPostprocessInput& input);

SternheimerFDHamiltonian::Complex accumulate_delta_sternheimer_response(
    const SternheimerFDHamiltonian::Vector& probe_wavefunction,
    const SternheimerDeltaPostprocessResult& response,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    double volume_element);

} // namespace ModuleRI

#endif
