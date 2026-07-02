#ifndef STERNHEIMER_FD_SOLVER_H
#define STERNHEIMER_FD_SOLVER_H

#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"
#include "source_lcao/module_ri/sternheimer_rpa.h"

#include <complex>
#include <vector>

namespace ModuleRI
{

struct SternheimerFDZeroOrderStates
{
    std::vector<double> eigenvalues;
    std::vector<SternheimerFDHamiltonian::Vector> wavefunctions;
    std::vector<double> residual_norms;
};

struct SternheimerFDLanczosOptions
{
    int max_subspace_size = 120;
    double residual_tolerance = 1.0e-8;
    unsigned int initial_seed = 1;
};

struct SternheimerFDLinearResponse
{
    SternheimerFDHamiltonian::Vector delta_wavefunction;
    SternheimerRPA::SolverResult solver;
    double residual_norm = 0.0;
};

SternheimerFDHamiltonian::Complex sternheimer_fd_grid_dot(const SternheimerFDHamiltonian::Vector& lhs,
                                                          const SternheimerFDHamiltonian::Vector& rhs,
                                                          double volume_element);

double sternheimer_fd_grid_norm(const SternheimerFDHamiltonian::Vector& wavefunction, double volume_element);

SternheimerFDZeroOrderStates solve_sternheimer_fd_zero_order_dense(const SternheimerFDHamiltonian& hamiltonian,
                                                                   int num_states,
                                                                   double volume_element,
                                                                   int max_size = 4096);

SternheimerFDZeroOrderStates solve_sternheimer_fd_zero_order_lanczos(
    const SternheimerFDHamiltonian& hamiltonian,
    int num_states,
    double volume_element,
    const SternheimerFDLanczosOptions& options = SternheimerFDLanczosOptions());

SternheimerFDLinearResponse solve_sternheimer_fd_linear_response(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    double reference_eigenvalue,
    const SternheimerFDHamiltonian::Vector& rhs,
    double omega,
    double volume_element,
    const SternheimerRPA::SolverOptions& options = SternheimerRPA::SolverOptions());

} // namespace ModuleRI

#endif
