#ifndef STERNHEIMER_PERIODIC_SOLVER_H
#define STERNHEIMER_PERIODIC_SOLVER_H

#include "source_lcao/module_ri/sternheimer_delta.h"
#include "source_lcao/module_ri/sternheimer_fd_solver.h"

namespace ModuleRI
{

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
    const SternheimerRPA::SolverOptions& options = SternheimerRPA::SolverOptions());

} // namespace ModuleRI

#endif
