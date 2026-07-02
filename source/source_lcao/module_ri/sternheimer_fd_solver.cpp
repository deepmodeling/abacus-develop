#include "source_lcao/module_ri/sternheimer_fd_solver.h"

#include <cmath>
#include <stdexcept>

namespace ModuleRI
{

SternheimerFDHamiltonian::Complex sternheimer_fd_grid_dot(const SternheimerFDHamiltonian::Vector& lhs,
                                                          const SternheimerFDHamiltonian::Vector& rhs,
                                                          const double volume_element)
{
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer FD grid dot requires a positive volume element.");
    }
    if (lhs.size() != rhs.size())
    {
        throw std::invalid_argument("Sternheimer FD grid dot vector sizes do not match.");
    }

    SternheimerFDHamiltonian::Complex value(0.0, 0.0);
    for (std::size_t ir = 0; ir != lhs.size(); ++ir)
    {
        value += std::conj(lhs[ir]) * rhs[ir];
    }
    return volume_element * value;
}

double sternheimer_fd_grid_norm(const SternheimerFDHamiltonian::Vector& wavefunction, const double volume_element)
{
    return std::sqrt(std::real(sternheimer_fd_grid_dot(wavefunction, wavefunction, volume_element)));
}

SternheimerFDZeroOrderStates solve_sternheimer_fd_zero_order_dense(const SternheimerFDHamiltonian& hamiltonian,
                                                                   const int num_states,
                                                                   const double volume_element,
                                                                   const int max_size)
{
    const int grid_size = hamiltonian.grid().size();
    if (num_states <= 0)
    {
        throw std::invalid_argument("Sternheimer FD zero-order solver requires a positive number of states.");
    }
    if (num_states > grid_size)
    {
        throw std::invalid_argument("Sternheimer FD zero-order solver requested more states than grid points.");
    }
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer FD zero-order solver requires a positive volume element.");
    }

    const auto eigenpairs = hamiltonian.diagonalize_dense(max_size);
    const double inv_sqrt_volume = 1.0 / std::sqrt(volume_element);

    SternheimerFDZeroOrderStates states;
    states.eigenvalues.assign(eigenpairs.eigenvalues.begin(), eigenpairs.eigenvalues.begin() + num_states);
    states.wavefunctions.assign(num_states, SternheimerFDHamiltonian::Vector(grid_size));
    states.residual_norms.assign(num_states, 0.0);

    for (int ib = 0; ib != num_states; ++ib)
    {
        for (int ir = 0; ir != grid_size; ++ir)
        {
            states.wavefunctions[ib][ir] = inv_sqrt_volume * eigenpairs.eigenvectors[ib][ir];
        }

        SternheimerFDHamiltonian::Vector hpsi;
        hamiltonian.apply(states.wavefunctions[ib], hpsi);
        SternheimerFDHamiltonian::Vector residual(grid_size);
        for (int ir = 0; ir != grid_size; ++ir)
        {
            residual[ir] = hpsi[ir] - states.eigenvalues[ib] * states.wavefunctions[ib][ir];
        }
        states.residual_norms[ib] = sternheimer_fd_grid_norm(residual, volume_element);
    }

    return states;
}

SternheimerFDLinearResponse solve_sternheimer_fd_linear_response(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    const double reference_eigenvalue,
    const SternheimerFDHamiltonian::Vector& rhs,
    const double omega,
    const double volume_element,
    const SternheimerRPA::SolverOptions& options)
{
    const int grid_size = hamiltonian.grid().size();
    if (static_cast<int>(rhs.size()) != grid_size)
    {
        throw std::invalid_argument("Sternheimer FD linear response rhs size does not match the grid.");
    }
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer FD linear response requires a positive volume element.");
    }
    for (const auto& occupied: occupied_wavefunctions)
    {
        if (static_cast<int>(occupied.size()) != grid_size)
        {
            throw std::invalid_argument("Sternheimer FD linear response occupied-state size does not match the grid.");
        }
    }

    auto dot = [volume_element](const SternheimerFDHamiltonian::Vector& lhs,
                                const SternheimerFDHamiltonian::Vector& rhs_vec) {
        return sternheimer_fd_grid_dot(lhs, rhs_vec, volume_element);
    };

    SternheimerFDHamiltonian::Vector projected_rhs = rhs;
    SternheimerRPA::project_out_subspace(occupied_wavefunctions, dot, projected_rhs);

    SternheimerRPA::LinearProblem problem;
    problem.dot = dot;
    problem.apply = [&hamiltonian, &occupied_wavefunctions, reference_eigenvalue, omega, dot](
                        const SternheimerFDHamiltonian::Vector& input,
                        SternheimerFDHamiltonian::Vector& output) {
        SternheimerFDHamiltonian::Vector pc_input = input;
        SternheimerRPA::project_out_subspace(occupied_wavefunctions, dot, pc_input);

        hamiltonian.apply(pc_input, output);
        const SternheimerFDHamiltonian::Complex shift(-reference_eigenvalue, omega);
        for (std::size_t ir = 0; ir != output.size(); ++ir)
        {
            output[ir] += shift * pc_input[ir];
        }
        SternheimerRPA::project_out_subspace(occupied_wavefunctions, dot, output);
    };

    SternheimerFDLinearResponse response;
    response.delta_wavefunction.assign(grid_size, SternheimerFDHamiltonian::Complex(0.0, 0.0));
    response.solver = SternheimerRPA::solve_bicgstab(problem, projected_rhs, response.delta_wavefunction, options);
    SternheimerRPA::project_out_subspace(occupied_wavefunctions, dot, response.delta_wavefunction);

    SternheimerFDHamiltonian::Vector applied;
    problem.apply(response.delta_wavefunction, applied);
    for (int ir = 0; ir != grid_size; ++ir)
    {
        applied[ir] -= projected_rhs[ir];
    }
    response.residual_norm = sternheimer_fd_grid_norm(applied, volume_element);
    return response;
}

} // namespace ModuleRI
