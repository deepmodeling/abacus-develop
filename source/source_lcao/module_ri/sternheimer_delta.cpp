#include "source_lcao/module_ri/sternheimer_delta.h"

#include "source_base/module_external/lapack_connector.h"
#include "source_lcao/module_ri/sternheimer_fd_solver.h"
#include "source_lcao/module_ri/sternheimer_rpa.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace ModuleRI
{
namespace
{

using Complex = SternheimerFDHamiltonian::Complex;
using Vector = SternheimerFDHamiltonian::Vector;

void check_vector_size(const Vector& vector, const std::size_t expected_size, const std::string& context)
{
    if (vector.size() != expected_size)
    {
        throw std::invalid_argument(context + " vector size does not match the Sternheimer delta grid.");
    }
}

void axpy(const Complex alpha, const Vector& x, Vector& y)
{
    check_vector_size(x, y.size(), "Sternheimer delta axpy");
    for (std::size_t ir = 0; ir != y.size(); ++ir)
    {
        y[ir] += alpha * x[ir];
    }
}

double difference_norm(const Vector& lhs, const Vector& rhs, const double volume_element)
{
    check_vector_size(lhs, rhs.size(), "Sternheimer delta difference_norm");
    Vector difference(lhs.size(), Complex(0.0, 0.0));
    for (std::size_t ir = 0; ir != lhs.size(); ++ir)
    {
        difference[ir] = lhs[ir] - rhs[ir];
    }
    return sternheimer_fd_grid_norm(difference, volume_element);
}

void scale_vector(Vector& vector, const Complex factor)
{
    for (Complex& value: vector)
    {
        value *= factor;
    }
}

std::vector<Vector> collect_virtual_orbitals(const std::vector<SternheimerDeltaVirtualState>& virtual_states)
{
    std::vector<Vector> orbitals;
    orbitals.reserve(virtual_states.size());
    for (const SternheimerDeltaVirtualState& state: virtual_states)
    {
        orbitals.push_back(state.orbital);
    }
    return orbitals;
}

std::vector<Vector> collect_fixed_subspace(const std::vector<Vector>& occupied_wavefunctions,
                                           const std::vector<SternheimerDeltaVirtualState>& virtual_states)
{
    std::vector<Vector> subspace = occupied_wavefunctions;
    subspace.reserve(occupied_wavefunctions.size() + virtual_states.size());
    for (const SternheimerDeltaVirtualState& state: virtual_states)
    {
        subspace.push_back(state.orbital);
    }
    return subspace;
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
        throw std::runtime_error("Sternheimer delta virtual-subspace diagonalization workspace query failed.");
    }
    const int lwork = std::max(1, static_cast<int>(std::ceil(work_query)));
    std::vector<double> work(lwork, 0.0);
    dsyev_(&jobz, &uplo, &size, matrix.data(), &lda, eigenvalues.data(), work.data(), &lwork, &info);
    if (info != 0)
    {
        throw std::runtime_error("Sternheimer delta virtual-subspace diagonalization failed.");
    }
    return eigenvalues;
}

void validate_postprocess_input(const Vector& standard_delta_wavefunction,
                                const SternheimerDeltaPostprocessInput& input)
{
    if (standard_delta_wavefunction.empty())
    {
        throw std::invalid_argument("Sternheimer delta postprocess requires a non-empty standard solution.");
    }
    if (input.volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta postprocess requires a positive grid volume element.");
    }
    if (input.virtual_states.size() != input.perturbation_matrix_elements.size())
    {
        throw std::invalid_argument(
            "Sternheimer delta postprocess requires one perturbation matrix element per virtual state.");
    }

    const std::size_t grid_size = standard_delta_wavefunction.size();
    for (const Vector& occupied: input.occupied_wavefunctions)
    {
        check_vector_size(occupied, grid_size, "Sternheimer delta occupied state");
    }
    for (const SternheimerDeltaVirtualState& state: input.virtual_states)
    {
        check_vector_size(state.orbital, grid_size, "Sternheimer delta virtual orbital");
        check_vector_size(state.residual, grid_size, "Sternheimer delta virtual residual");
    }
}

void validate_virtual_states(const std::vector<SternheimerDeltaVirtualState>& virtual_states,
                             const std::size_t grid_size)
{
    for (const SternheimerDeltaVirtualState& state: virtual_states)
    {
        check_vector_size(state.orbital, grid_size, "Sternheimer delta virtual orbital");
        check_vector_size(state.residual, grid_size, "Sternheimer delta virtual residual");
    }
}

} // namespace

SternheimerDeltaPostprocessResult postprocess_delta_sternheimer_solution(
    const SternheimerFDHamiltonian::Vector& standard_delta_wavefunction,
    const SternheimerDeltaPostprocessInput& input)
{
    validate_postprocess_input(standard_delta_wavefunction, input);

    auto dot = [&input](const Vector& lhs, const Vector& rhs) {
        return sternheimer_fd_grid_dot(lhs, rhs, input.volume_element);
    };

    SternheimerDeltaPostprocessResult result;
    result.out_wavefunction = standard_delta_wavefunction;
    SternheimerRPA::project_out_subspace(input.occupied_wavefunctions, dot, result.out_wavefunction);
    const std::vector<Vector> virtual_orbitals = collect_virtual_orbitals(input.virtual_states);
    SternheimerRPA::project_out_subspace(virtual_orbitals, dot, result.out_wavefunction);

    result.coefficients.assign(input.virtual_states.size(), Complex(0.0, 0.0));
    for (std::size_t ia = 0; ia != input.virtual_states.size(); ++ia)
    {
        const SternheimerDeltaVirtualState& state = input.virtual_states[ia];
        const Complex denominator(input.occupied_eigenvalue - state.eigenvalue, -input.omega);
        if (std::abs(denominator) < 1.0e-30)
        {
            throw std::runtime_error("Sternheimer delta postprocess found a singular virtual-state denominator.");
        }
        const Complex residual_overlap = sternheimer_fd_grid_dot(state.residual,
                                                                 result.out_wavefunction,
                                                                 input.volume_element);
        result.coefficients[ia] = (input.perturbation_matrix_elements[ia] + residual_overlap) / denominator;
    }

    result.reconstructed_wavefunction = result.out_wavefunction;
    for (std::size_t ia = 0; ia != input.virtual_states.size(); ++ia)
    {
        axpy(result.coefficients[ia], input.virtual_states[ia].orbital, result.reconstructed_wavefunction);
    }
    result.out_norm = sternheimer_fd_grid_norm(result.out_wavefunction, input.volume_element);
    result.reconstruction_error
        = difference_norm(result.reconstructed_wavefunction, standard_delta_wavefunction, input.volume_element);
    return result;
}

SternheimerDeltaSubspace build_delta_sternheimer_subspace(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    const std::vector<SternheimerFDHamiltonian::Vector>& candidate_orbitals,
    const double volume_element,
    const SternheimerDeltaSubspaceOptions& options)
{
    const int grid_size = hamiltonian.grid().size();
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta subspace requires a positive grid volume element.");
    }
    if (options.max_virtual_states < 0)
    {
        throw std::invalid_argument("Sternheimer delta subspace max_virtual_states must be non-negative.");
    }
    if (options.norm_tolerance < 0.0)
    {
        throw std::invalid_argument("Sternheimer delta subspace norm_tolerance must be non-negative.");
    }
    for (const Vector& occupied: occupied_wavefunctions)
    {
        check_vector_size(occupied, static_cast<std::size_t>(grid_size), "Sternheimer delta occupied state");
    }

    auto dot = [volume_element](const Vector& lhs, const Vector& rhs) {
        return sternheimer_fd_grid_dot(lhs, rhs, volume_element);
    };

    std::vector<Vector> orthonormal_candidates;
    std::vector<Vector> projection_subspace = occupied_wavefunctions;
    for (Vector candidate: candidate_orbitals)
    {
        check_vector_size(candidate, static_cast<std::size_t>(grid_size), "Sternheimer delta candidate orbital");
        SternheimerRPA::project_out_subspace(projection_subspace, dot, candidate);
        const double norm = sternheimer_fd_grid_norm(candidate, volume_element);
        if (norm <= options.norm_tolerance)
        {
            continue;
        }
        scale_vector(candidate, Complex(1.0 / norm, 0.0));
        projection_subspace.push_back(candidate);
        orthonormal_candidates.push_back(std::move(candidate));
        if (options.max_virtual_states > 0
            && static_cast<int>(orthonormal_candidates.size()) >= options.max_virtual_states)
        {
            break;
        }
    }

    SternheimerDeltaSubspace subspace;
    subspace.accepted_candidates = static_cast<int>(orthonormal_candidates.size());
    subspace.discarded_candidates
        = static_cast<int>(candidate_orbitals.size()) - subspace.accepted_candidates;
    if (orthonormal_candidates.empty())
    {
        return subspace;
    }

    const int nvirtual = static_cast<int>(orthonormal_candidates.size());
    std::vector<Vector> h_candidates(nvirtual);
    for (int j = 0; j != nvirtual; ++j)
    {
        hamiltonian.apply(orthonormal_candidates[j], h_candidates[j]);
    }

    std::vector<double> h_matrix(static_cast<std::size_t>(nvirtual) * static_cast<std::size_t>(nvirtual), 0.0);
    for (int j = 0; j != nvirtual; ++j)
    {
        for (int i = 0; i != nvirtual; ++i)
        {
            h_matrix[static_cast<std::size_t>(i) + static_cast<std::size_t>(nvirtual) * static_cast<std::size_t>(j)]
                = std::real(dot(orthonormal_candidates[i], h_candidates[j]));
        }
    }
    const std::vector<double> eigenvalues = diagonalize_real_symmetric(h_matrix, nvirtual);

    subspace.virtual_states.reserve(static_cast<std::size_t>(nvirtual));
    for (int ia = 0; ia != nvirtual; ++ia)
    {
        Vector eta(static_cast<std::size_t>(grid_size), Complex(0.0, 0.0));
        for (int j = 0; j != nvirtual; ++j)
        {
            const double coefficient
                = h_matrix[static_cast<std::size_t>(j) + static_cast<std::size_t>(nvirtual) * static_cast<std::size_t>(ia)];
            axpy(Complex(coefficient, 0.0), orthonormal_candidates[j], eta);
        }

        SternheimerRPA::project_out_subspace(occupied_wavefunctions, dot, eta);
        const double eta_norm = sternheimer_fd_grid_norm(eta, volume_element);
        if (eta_norm <= options.norm_tolerance)
        {
            ++subspace.discarded_candidates;
            continue;
        }
        scale_vector(eta, Complex(1.0 / eta_norm, 0.0));

        Vector residual;
        hamiltonian.apply(eta, residual);
        for (std::size_t ir = 0; ir != residual.size(); ++ir)
        {
            residual[ir] -= eigenvalues[static_cast<std::size_t>(ia)] * eta[ir];
        }

        SternheimerDeltaVirtualState state;
        state.orbital = std::move(eta);
        state.eigenvalue = eigenvalues[static_cast<std::size_t>(ia)];
        subspace.virtual_states.push_back(std::move(state));
    }

    const std::vector<Vector> virtual_orbitals = collect_virtual_orbitals(subspace.virtual_states);
    std::vector<Vector> residual_projector = occupied_wavefunctions;
    residual_projector.insert(residual_projector.end(), virtual_orbitals.begin(), virtual_orbitals.end());
    for (SternheimerDeltaVirtualState& state: subspace.virtual_states)
    {
        hamiltonian.apply(state.orbital, state.residual);
        for (std::size_t ir = 0; ir != state.residual.size(); ++ir)
        {
            state.residual[ir] -= state.eigenvalue * state.orbital[ir];
        }
        SternheimerRPA::project_out_subspace(residual_projector, dot, state.residual);
    }
    return subspace;
}

std::vector<SternheimerFDHamiltonian::Complex> delta_sternheimer_perturbation_matrix_elements(
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<double>& perturbation_potential,
    const SternheimerFDHamiltonian::Vector& occupied_wavefunction,
    const double volume_element)
{
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta perturbation elements require a positive grid volume element.");
    }
    if (perturbation_potential.size() != occupied_wavefunction.size())
    {
        throw std::invalid_argument("Sternheimer delta perturbation potential size does not match wavefunction.");
    }
    validate_virtual_states(virtual_states, occupied_wavefunction.size());

    std::vector<Complex> elements(virtual_states.size(), Complex(0.0, 0.0));
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        Complex value(0.0, 0.0);
        for (std::size_t ir = 0; ir != occupied_wavefunction.size(); ++ir)
        {
            value += std::conj(virtual_states[ia].orbital[ir]) * perturbation_potential[ir]
                * occupied_wavefunction[ir];
        }
        elements[ia] = volume_element * value;
    }
    return elements;
}

SternheimerDeltaLinearResponse solve_delta_sternheimer_linear_response(
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
    const int grid_size = hamiltonian.grid().size();
    check_vector_size(rhs, static_cast<std::size_t>(grid_size), "Sternheimer delta rhs");
    if (virtual_states.size() != perturbation_matrix_elements.size())
    {
        throw std::invalid_argument(
            "Sternheimer delta linear response requires one perturbation matrix element per virtual state.");
    }
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta linear response requires a positive grid volume element.");
    }
    for (const Vector& occupied: occupied_wavefunctions)
    {
        check_vector_size(occupied, static_cast<std::size_t>(grid_size), "Sternheimer delta occupied state");
    }
    validate_virtual_states(virtual_states, static_cast<std::size_t>(grid_size));

    auto dot = [volume_element](const Vector& lhs, const Vector& rhs_vec) {
        return sternheimer_fd_grid_dot(lhs, rhs_vec, volume_element);
    };
    const std::vector<Vector> fixed_subspace = collect_fixed_subspace(occupied_wavefunctions, virtual_states);

    std::vector<Complex> denominators(virtual_states.size(), Complex(0.0, 0.0));
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        denominators[ia] = Complex(reference_eigenvalue - virtual_states[ia].eigenvalue, -omega);
        if (std::abs(denominators[ia]) < 1.0e-30)
        {
            throw std::runtime_error("Sternheimer delta linear response found a singular virtual-state denominator.");
        }
    }

    Vector projected_rhs = rhs;
    SternheimerRPA::project_out_subspace(fixed_subspace, dot, projected_rhs);
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        axpy(-perturbation_matrix_elements[ia] / denominators[ia], virtual_states[ia].residual, projected_rhs);
    }

    SternheimerRPA::LinearProblem problem;
    problem.dot = dot;
    problem.apply = [&hamiltonian, &fixed_subspace, &virtual_states, &denominators, reference_eigenvalue, omega, dot](
                        const Vector& input,
                        Vector& output) {
        Vector q_input = input;
        SternheimerRPA::project_out_subspace(fixed_subspace, dot, q_input);

        hamiltonian.apply(q_input, output);
        const Complex shift(-reference_eigenvalue, omega);
        for (std::size_t ir = 0; ir != output.size(); ++ir)
        {
            output[ir] += shift * q_input[ir];
        }
        SternheimerRPA::project_out_subspace(fixed_subspace, dot, output);
        for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
        {
            const Complex coupling = dot(virtual_states[ia].residual, q_input) / denominators[ia];
            axpy(coupling, virtual_states[ia].residual, output);
        }
    };

    SternheimerDeltaLinearResponse result;
    result.response.out_wavefunction.assign(static_cast<std::size_t>(grid_size), Complex(0.0, 0.0));
    result.solver = SternheimerRPA::solve_gmres(problem,
                                                projected_rhs,
                                                result.response.out_wavefunction,
                                                options);
    SternheimerRPA::project_out_subspace(fixed_subspace, dot, result.response.out_wavefunction);

    result.response.coefficients.assign(virtual_states.size(), Complex(0.0, 0.0));
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        const Complex residual_overlap = dot(virtual_states[ia].residual, result.response.out_wavefunction);
        result.response.coefficients[ia]
            = (perturbation_matrix_elements[ia] + residual_overlap) / denominators[ia];
    }

    result.response.reconstructed_wavefunction = result.response.out_wavefunction;
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        axpy(result.response.coefficients[ia], virtual_states[ia].orbital, result.response.reconstructed_wavefunction);
    }
    result.response.out_norm = sternheimer_fd_grid_norm(result.response.out_wavefunction, volume_element);

    Vector applied;
    hamiltonian.apply(result.response.reconstructed_wavefunction, applied);
    const Complex shift(-reference_eigenvalue, omega);
    for (std::size_t ir = 0; ir != applied.size(); ++ir)
    {
        applied[ir] += shift * result.response.reconstructed_wavefunction[ir];
    }
    SternheimerRPA::project_out_subspace(occupied_wavefunctions, dot, applied);

    Vector occupied_projected_rhs = rhs;
    SternheimerRPA::project_out_subspace(occupied_wavefunctions, dot, occupied_projected_rhs);
    for (std::size_t ir = 0; ir != applied.size(); ++ir)
    {
        applied[ir] -= occupied_projected_rhs[ir];
    }
    result.residual_norm = sternheimer_fd_grid_norm(applied, volume_element);
    result.response.reconstruction_error = result.residual_norm;
    return result;
}

SternheimerFDHamiltonian::Complex accumulate_delta_sternheimer_response(
    const SternheimerFDHamiltonian::Vector& probe_wavefunction,
    const SternheimerDeltaPostprocessResult& response,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const double volume_element)
{
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta response requires a positive grid volume element.");
    }
    if (response.coefficients.size() != virtual_states.size())
    {
        throw std::invalid_argument("Sternheimer delta response coefficient count does not match virtual states.");
    }
    check_vector_size(response.out_wavefunction, probe_wavefunction.size(), "Sternheimer delta response out state");

    Complex value = sternheimer_fd_grid_dot(probe_wavefunction, response.out_wavefunction, volume_element);
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        check_vector_size(virtual_states[ia].orbital, probe_wavefunction.size(), "Sternheimer delta response virtual");
        value += response.coefficients[ia]
            * sternheimer_fd_grid_dot(probe_wavefunction, virtual_states[ia].orbital, volume_element);
    }
    return value;
}

} // namespace ModuleRI
