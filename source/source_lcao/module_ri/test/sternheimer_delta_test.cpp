#include "source_lcao/module_ri/sternheimer_delta.h"

#include "source_lcao/module_ri/sternheimer_fd_solver.h"

#include <complex>
#include <gtest/gtest.h>
#include <vector>

namespace
{

using Complex = std::complex<double>;
using Vector = ModuleRI::SternheimerFDHamiltonian::Vector;

Vector unit_vector(const int size, const int index)
{
    Vector vec(size, Complex(0.0, 0.0));
    vec[index] = Complex(1.0, 0.0);
    return vec;
}

void expect_vector_near(const Vector& actual, const Vector& expected, const double tolerance)
{
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t ir = 0; ir != actual.size(); ++ir)
    {
        EXPECT_NEAR(actual[ir].real(), expected[ir].real(), tolerance);
        EXPECT_NEAR(actual[ir].imag(), expected[ir].imag(), tolerance);
    }
}

} // namespace

TEST(SternheimerDelta, PostprocessReconstructsStandardSolutionWithResidualCoupling)
{
    constexpr double volume_element = 1.0;
    constexpr double eps_i = -0.25;
    constexpr double eps_a = 0.75;
    constexpr double omega = 0.4;

    const Vector occupied = unit_vector(4, 0);
    const Vector eta = unit_vector(4, 1);
    const Vector residual = {Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.30, -0.10), Complex(-0.20, 0.40)};
    const Vector out = {Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.20, 0.30), Complex(-0.10, 0.25)};
    const Complex expected_coefficient(0.70, -0.20);

    Vector standard_delta = out;
    for (std::size_t ir = 0; ir != standard_delta.size(); ++ir)
    {
        standard_delta[ir] += expected_coefficient * eta[ir];
    }

    const Complex residual_overlap = ModuleRI::sternheimer_fd_grid_dot(residual, out, volume_element);
    const Complex denominator(eps_i - eps_a, -omega);
    const Complex perturbation_matrix_element = expected_coefficient * denominator - residual_overlap;

    ModuleRI::SternheimerDeltaPostprocessInput input;
    input.occupied_wavefunctions = {occupied};
    input.virtual_states = {{eta, residual, eps_a}};
    input.perturbation_matrix_elements = {perturbation_matrix_element};
    input.occupied_eigenvalue = eps_i;
    input.omega = omega;
    input.volume_element = volume_element;

    const auto result = ModuleRI::postprocess_delta_sternheimer_solution(standard_delta, input);

    ASSERT_EQ(result.coefficients.size(), 1);
    EXPECT_NEAR(result.coefficients[0].real(), expected_coefficient.real(), 1.0e-14);
    EXPECT_NEAR(result.coefficients[0].imag(), expected_coefficient.imag(), 1.0e-14);
    expect_vector_near(result.out_wavefunction, out, 1.0e-14);
    expect_vector_near(result.reconstructed_wavefunction, standard_delta, 1.0e-14);
    EXPECT_NEAR(result.reconstruction_error, 0.0, 1.0e-14);

    EXPECT_NEAR(std::abs(ModuleRI::sternheimer_fd_grid_dot(occupied, result.out_wavefunction, volume_element)),
                0.0,
                1.0e-14);
    EXPECT_NEAR(std::abs(ModuleRI::sternheimer_fd_grid_dot(eta, result.out_wavefunction, volume_element)),
                0.0,
                1.0e-14);
}

TEST(SternheimerDelta, AccumulateResponseFromDecomposedWavefunction)
{
    constexpr double volume_element = 1.0;
    const Vector eta = unit_vector(4, 1);
    const Vector residual(4, Complex(0.0, 0.0));
    const Vector out = {Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.4, -0.1), Complex(-0.3, 0.2)};
    const Complex coefficient(-0.6, 0.5);

    ModuleRI::SternheimerDeltaPostprocessResult result;
    result.out_wavefunction = out;
    result.coefficients = {coefficient};
    result.reconstructed_wavefunction = out;
    for (std::size_t ir = 0; ir != result.reconstructed_wavefunction.size(); ++ir)
    {
        result.reconstructed_wavefunction[ir] += coefficient * eta[ir];
    }

    const Vector probe = {Complex(0.0, 0.0), Complex(1.5, -0.2), Complex(-0.4, 0.1), Complex(0.7, 0.3)};
    const std::vector<ModuleRI::SternheimerDeltaVirtualState> virtual_states = {{eta, residual, 0.75}};

    const Complex decomposed
        = ModuleRI::accumulate_delta_sternheimer_response(probe, result, virtual_states, volume_element);
    const Complex direct = ModuleRI::sternheimer_fd_grid_dot(probe, result.reconstructed_wavefunction, volume_element);

    EXPECT_NEAR(decomposed.real(), direct.real(), 1.0e-14);
    EXPECT_NEAR(decomposed.imag(), direct.imag(), 1.0e-14);
}

TEST(SternheimerDelta, BuildsOrthogonalVirtualSubspaceWithHamiltonianResiduals)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid{4, 1, 1, 0.5, 1.0, 1.0, true};
    constexpr double volume_element = 0.5;
    ModuleRI::SternheimerFDHamiltonian hamiltonian(grid, std::vector<double>(grid.size(), 0.25));
    const auto states = ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 4, volume_element);

    Vector candidate = states.wavefunctions[2];
    for (std::size_t ir = 0; ir != candidate.size(); ++ir)
    {
        candidate[ir] += Complex(0.3, 0.0) * states.wavefunctions[3][ir];
    }

    ModuleRI::SternheimerDeltaSubspaceOptions options;
    options.max_virtual_states = 1;
    const auto subspace = ModuleRI::build_delta_sternheimer_subspace(hamiltonian,
                                                                     {states.wavefunctions[0]},
                                                                     {candidate},
                                                                     volume_element,
                                                                     options);

    ASSERT_EQ(subspace.virtual_states.size(), 1);
    EXPECT_EQ(subspace.accepted_candidates, 1);
    EXPECT_EQ(subspace.discarded_candidates, 0);
    EXPECT_NEAR(ModuleRI::sternheimer_fd_grid_norm(subspace.virtual_states[0].orbital, volume_element),
                1.0,
                1.0e-12);
    EXPECT_NEAR(std::abs(ModuleRI::sternheimer_fd_grid_dot(states.wavefunctions[0],
                                                           subspace.virtual_states[0].orbital,
                                                           volume_element)),
                0.0,
                1.0e-12);

    Vector h_eta;
    hamiltonian.apply(subspace.virtual_states[0].orbital, h_eta);
    for (std::size_t ir = 0; ir != h_eta.size(); ++ir)
    {
        h_eta[ir] -= subspace.virtual_states[0].eigenvalue * subspace.virtual_states[0].orbital[ir];
    }
    ModuleRI::SternheimerRPA::project_out_subspace({states.wavefunctions[0], subspace.virtual_states[0].orbital},
                                                   [volume_element](const Vector& lhs, const Vector& rhs) {
                                                       return ModuleRI::sternheimer_fd_grid_dot(lhs,
                                                                                                rhs,
                                                                                                volume_element);
                                                   },
                                                   h_eta);
    expect_vector_near(subspace.virtual_states[0].residual, h_eta, 1.0e-12);
}

TEST(SternheimerDelta, StrictProjectedSolverMatchesStandardSternheimerResponse)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid{5, 1, 1, 0.6, 1.0, 1.0, true};
    constexpr double volume_element = 0.6;
    const std::vector<double> potential = {0.10, -0.05, 0.20, -0.10, 0.05};
    ModuleRI::SternheimerFDHamiltonian hamiltonian(grid, potential);
    const auto states = ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 5, volume_element);

    Vector candidate = states.wavefunctions[2];
    for (std::size_t ir = 0; ir != candidate.size(); ++ir)
    {
        candidate[ir] += Complex(0.25, -0.10) * states.wavefunctions[3][ir];
    }
    ModuleRI::SternheimerDeltaSubspaceOptions subspace_options;
    subspace_options.max_virtual_states = 1;
    const auto subspace = ModuleRI::build_delta_sternheimer_subspace(hamiltonian,
                                                                     {states.wavefunctions[0]},
                                                                     {candidate},
                                                                     volume_element,
                                                                     subspace_options);
    ASSERT_EQ(subspace.virtual_states.size(), 1);

    const std::vector<double> perturbation = {0.7, -0.3, 0.2, 0.5, -0.4};
    Vector rhs;
    ModuleRI::SternheimerRPA::build_rhs_from_hartree_perturbation(perturbation, states.wavefunctions[0], rhs);

    ModuleRI::SternheimerRPA::SolverOptions solver_options;
    solver_options.max_iter = 80;
    solver_options.residual_tol = 1.0e-12;
    constexpr double omega = 0.45;

    const auto standard = ModuleRI::solve_sternheimer_fd_linear_response(hamiltonian,
                                                                         {states.wavefunctions[0]},
                                                                         states.eigenvalues[0],
                                                                         rhs,
                                                                         omega,
                                                                         volume_element,
                                                                         solver_options);
    ASSERT_TRUE(standard.solver.converged);

    const std::vector<Complex> perturbation_matrix_elements
        = ModuleRI::delta_sternheimer_perturbation_matrix_elements(subspace.virtual_states,
                                                                   perturbation,
                                                                   states.wavefunctions[0],
                                                                   volume_element);
    const auto delta = ModuleRI::solve_delta_sternheimer_linear_response(hamiltonian,
                                                                         {states.wavefunctions[0]},
                                                                         states.eigenvalues[0],
                                                                         rhs,
                                                                         subspace.virtual_states,
                                                                         perturbation_matrix_elements,
                                                                         omega,
                                                                         volume_element,
                                                                         solver_options);

    ASSERT_TRUE(delta.solver.converged);
    EXPECT_LT(delta.residual_norm, 1.0e-9);
    expect_vector_near(delta.response.reconstructed_wavefunction, standard.delta_wavefunction, 1.0e-8);
}
