#include "source_lcao/module_ri/sternheimer_fd_solver.h"

#include "source_base/module_external/lapack_connector.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

namespace ModuleRI
{
namespace
{

using Complex = SternheimerFDHamiltonian::Complex;
using Vector = SternheimerFDHamiltonian::Vector;

void scale_vector(Vector& vector, const Complex factor)
{
    for (Complex& value: vector)
    {
        value *= factor;
    }
}

void axpy(const Complex alpha, const Vector& x, Vector& y)
{
    if (x.size() != y.size())
    {
        throw std::invalid_argument("Sternheimer FD vector axpy sizes do not match.");
    }
    for (std::size_t ir = 0; ir != x.size(); ++ir)
    {
        y[ir] += alpha * x[ir];
    }
}

void normalize_vector(Vector& vector, const double volume_element)
{
    const double norm = sternheimer_fd_grid_norm(vector, volume_element);
    if (norm <= 0.0)
    {
        throw std::runtime_error("Sternheimer FD Lanczos generated a zero vector.");
    }
    scale_vector(vector, Complex(1.0 / norm, 0.0));
}

std::vector<double> diagonalize_real_symmetric(std::vector<double>& matrix, const int size)
{
    char jobz = 'V';
    char uplo = 'U';
    const int lda = size;
    int info = 0;
    std::vector<double> eigenvalues(size, 0.0);
    double work_query = 0.0;
    const int minus_one = -1;

    dsyev_(&jobz, &uplo, &size, matrix.data(), &lda, eigenvalues.data(), &work_query, &minus_one, &info);
    if (info != 0)
    {
        throw std::runtime_error("Sternheimer FD Lanczos Ritz diagonalization workspace query failed.");
    }

    const int lwork = std::max(1, static_cast<int>(std::ceil(work_query)));
    std::vector<double> work(lwork, 0.0);
    dsyev_(&jobz, &uplo, &size, matrix.data(), &lda, eigenvalues.data(), work.data(), &lwork, &info);
    if (info != 0)
    {
        throw std::runtime_error("Sternheimer FD Lanczos Ritz diagonalization failed.");
    }
    return eigenvalues;
}

SternheimerFDZeroOrderStates build_zero_order_states_from_wavefunctions(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<double>& eigenvalues,
    std::vector<Vector> wavefunctions,
    const double volume_element)
{
    SternheimerFDZeroOrderStates states;
    states.eigenvalues = eigenvalues;
    states.wavefunctions = std::move(wavefunctions);
    states.residual_norms.assign(states.eigenvalues.size(), 0.0);

    for (std::size_t ib = 0; ib != states.wavefunctions.size(); ++ib)
    {
        normalize_vector(states.wavefunctions[ib], volume_element);

        Vector hpsi;
        hamiltonian.apply(states.wavefunctions[ib], hpsi);
        Vector residual(hpsi.size());
        for (std::size_t ir = 0; ir != hpsi.size(); ++ir)
        {
            residual[ir] = hpsi[ir] - states.eigenvalues[ib] * states.wavefunctions[ib][ir];
        }
        states.residual_norms[ib] = sternheimer_fd_grid_norm(residual, volume_element);
    }
    return states;
}

} // namespace

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

double sternheimer_fd_linear_response_residual_norm(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    const double reference_eigenvalue,
    const SternheimerFDHamiltonian::Vector& rhs,
    const SternheimerFDHamiltonian::Vector& delta_wavefunction,
    const double omega,
    const double volume_element)
{
    const std::size_t grid_size = static_cast<std::size_t>(hamiltonian.grid().size());
    if (rhs.size() != grid_size || delta_wavefunction.size() != grid_size)
    {
        throw std::invalid_argument("Sternheimer FD residual vectors do not match the grid.");
    }
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer FD residual requires a positive volume element.");
    }
    for (const Vector& occupied: occupied_wavefunctions)
    {
        if (occupied.size() != grid_size)
        {
            throw std::invalid_argument("Sternheimer FD residual occupied state does not match the grid.");
        }
    }

    auto dot = [volume_element](const Vector& lhs, const Vector& rhs_vec) {
        return sternheimer_fd_grid_dot(lhs, rhs_vec, volume_element);
    };
    Vector projected_wavefunction = delta_wavefunction;
    SternheimerRPA::project_out_subspace(occupied_wavefunctions, dot, projected_wavefunction);

    Vector residual;
    hamiltonian.apply(projected_wavefunction, residual);
    const Complex shift(-reference_eigenvalue, omega);
    for (std::size_t ir = 0; ir != residual.size(); ++ir)
    {
        residual[ir] += shift * projected_wavefunction[ir];
    }
    SternheimerRPA::project_out_subspace(occupied_wavefunctions, dot, residual);

    Vector projected_rhs = rhs;
    SternheimerRPA::project_out_subspace(occupied_wavefunctions, dot, projected_rhs);
    for (std::size_t ir = 0; ir != residual.size(); ++ir)
    {
        residual[ir] -= projected_rhs[ir];
    }
    return sternheimer_fd_grid_norm(residual, volume_element);
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

SternheimerFDZeroOrderStates solve_sternheimer_fd_zero_order_lanczos(
    const SternheimerFDHamiltonian& hamiltonian,
    const int num_states,
    const double volume_element,
    const SternheimerFDLanczosOptions& options)
{
    const int grid_size = hamiltonian.grid().size();
    if (num_states <= 0)
    {
        throw std::invalid_argument("Sternheimer FD Lanczos zero-order solver requires a positive number of states.");
    }
    if (num_states > grid_size)
    {
        throw std::invalid_argument("Sternheimer FD Lanczos zero-order solver requested more states than grid points.");
    }
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer FD Lanczos zero-order solver requires a positive volume element.");
    }
    if (options.max_subspace_size < num_states)
    {
        throw std::invalid_argument(
            "Sternheimer FD Lanczos zero-order solver requires max_subspace_size >= num_states.");
    }
    if (options.residual_tolerance < 0.0)
    {
        throw std::invalid_argument(
            "Sternheimer FD Lanczos zero-order solver requires a non-negative residual tolerance.");
    }

    const int subspace_limit = std::min(options.max_subspace_size, grid_size);
    std::mt19937 generator(options.initial_seed);
    std::uniform_real_distribution<double> distribution(-0.5, 0.5);

    Vector q(grid_size, Complex(0.0, 0.0));
    for (Complex& value: q)
    {
        value = Complex(distribution(generator), distribution(generator));
    }
    normalize_vector(q, volume_element);

    std::vector<Vector> basis;
    basis.reserve(subspace_limit);
    std::vector<double> alpha;
    alpha.reserve(subspace_limit);
    std::vector<double> beta;
    beta.reserve(std::max(0, subspace_limit - 1));

    Vector previous(grid_size, Complex(0.0, 0.0));
    double beta_previous = 0.0;
    constexpr double breakdown_tolerance = 1.0e-14;

    for (int iter = 0; iter != subspace_limit; ++iter)
    {
        basis.push_back(q);

        Vector work;
        hamiltonian.apply(q, work);
        if (iter > 0)
        {
            axpy(Complex(-beta_previous, 0.0), previous, work);
        }

        const double alpha_value = std::real(sternheimer_fd_grid_dot(q, work, volume_element));
        alpha.push_back(alpha_value);
        axpy(Complex(-alpha_value, 0.0), q, work);

        for (int pass = 0; pass != 2; ++pass)
        {
            for (const Vector& basis_vector: basis)
            {
                const Complex overlap = sternheimer_fd_grid_dot(basis_vector, work, volume_element);
                axpy(-overlap, basis_vector, work);
            }
        }

        const double beta_value = sternheimer_fd_grid_norm(work, volume_element);
        if (iter + 1 == subspace_limit || beta_value < breakdown_tolerance)
        {
            break;
        }

        beta.push_back(beta_value);
        previous = q;
        q = std::move(work);
        scale_vector(q, Complex(1.0 / beta_value, 0.0));
        beta_previous = beta_value;
    }

    const int subspace_size = static_cast<int>(basis.size());
    if (subspace_size < num_states)
    {
        throw std::runtime_error("Sternheimer FD Lanczos subspace is smaller than the requested state count.");
    }

    std::vector<double> ritz_matrix(subspace_size * subspace_size, 0.0);
    for (int i = 0; i != subspace_size; ++i)
    {
        ritz_matrix[i + subspace_size * i] = alpha[i];
    }
    for (int i = 0; i + 1 != subspace_size; ++i)
    {
        ritz_matrix[i + (i + 1) * subspace_size] = beta[i];
        ritz_matrix[i + 1 + i * subspace_size] = beta[i];
    }

    const std::vector<double> ritz_values = diagonalize_real_symmetric(ritz_matrix, subspace_size);

    std::vector<double> eigenvalues(num_states, 0.0);
    std::vector<Vector> wavefunctions(num_states, Vector(grid_size, Complex(0.0, 0.0)));
    for (int ib = 0; ib != num_states; ++ib)
    {
        eigenvalues[ib] = ritz_values[ib];
        for (int j = 0; j != subspace_size; ++j)
        {
            const double coefficient = ritz_matrix[j + ib * subspace_size];
            axpy(Complex(coefficient, 0.0), basis[j], wavefunctions[ib]);
        }
    }

    return build_zero_order_states_from_wavefunctions(hamiltonian, eigenvalues, std::move(wavefunctions), volume_element);
}

SternheimerFDHamiltonian::Vector build_sternheimer_fd_complete_sos_response(
    const SternheimerFDZeroOrderStates& complete_states,
    const int occupied_state_count,
    const int occupied_state_index,
    const SternheimerFDHamiltonian::Vector& rhs,
    const double omega,
    const double volume_element)
{
    const std::size_t grid_size = rhs.size();
    if (grid_size == 0 || complete_states.wavefunctions.size() != grid_size
        || complete_states.eigenvalues.size() != grid_size)
    {
        throw std::invalid_argument(
            "Sternheimer complete SOS reference requires one eigenpair per finite-grid degree of freedom.");
    }
    if (occupied_state_count <= 0 || occupied_state_count >= static_cast<int>(grid_size)
        || occupied_state_index < 0 || occupied_state_index >= occupied_state_count)
    {
        throw std::invalid_argument("Sternheimer complete SOS reference has invalid occupied-state indices.");
    }
    if (omega < 0.0 || volume_element <= 0.0)
    {
        throw std::invalid_argument(
            "Sternheimer complete SOS reference requires non-negative omega and a positive volume element.");
    }
    for (const auto& wavefunction: complete_states.wavefunctions)
    {
        if (wavefunction.size() != grid_size)
        {
            throw std::invalid_argument("Sternheimer complete SOS eigenvector size does not match the rhs.");
        }
    }

    SternheimerFDHamiltonian::Vector response(grid_size, SternheimerFDHamiltonian::Complex(0.0, 0.0));
    const double occupied_eigenvalue
        = complete_states.eigenvalues[static_cast<std::size_t>(occupied_state_index)];
    for (int state_index = occupied_state_count; state_index != static_cast<int>(grid_size); ++state_index)
    {
        const auto& virtual_state = complete_states.wavefunctions[static_cast<std::size_t>(state_index)];
        const SternheimerFDHamiltonian::Complex denominator(
            complete_states.eigenvalues[static_cast<std::size_t>(state_index)] - occupied_eigenvalue,
            omega);
        if (std::abs(denominator) < 1.0e-30)
        {
            throw std::runtime_error("Sternheimer complete SOS reference found a singular denominator.");
        }
        const SternheimerFDHamiltonian::Complex coefficient
            = sternheimer_fd_grid_dot(virtual_state, rhs, volume_element) / denominator;
        axpy(coefficient, virtual_state, response);
    }
    return response;
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
    response.solver = SternheimerRPA::solve_gmres(problem, projected_rhs, response.delta_wavefunction, options);
    SternheimerRPA::project_out_subspace(occupied_wavefunctions, dot, response.delta_wavefunction);

    response.residual_norm = sternheimer_fd_linear_response_residual_norm(hamiltonian,
                                                                           occupied_wavefunctions,
                                                                           reference_eigenvalue,
                                                                           rhs,
                                                                           response.delta_wavefunction,
                                                                           omega,
                                                                           volume_element);
    return response;
}

} // namespace ModuleRI
