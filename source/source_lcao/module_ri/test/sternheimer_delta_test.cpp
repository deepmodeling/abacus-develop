#include "source_lcao/module_ri/sternheimer_delta.h"

#include "source_lcao/module_ri/sternheimer_fd_solver.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <memory>
#include <string>
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

std::vector<Complex> column_major_matrix_vector_product(const std::vector<Complex>& matrix,
                                                        const std::vector<Complex>& vector)
{
    const std::size_t size = vector.size();
    std::vector<Complex> result(size, Complex(0.0, 0.0));
    for (std::size_t column = 0; column != size; ++column)
    {
        for (std::size_t row = 0; row != size; ++row)
        {
            result[row] += matrix[row + size * column] * vector[column];
        }
    }
    return result;
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
    ASSERT_EQ(result.sos_coefficients.size(), 1);
    ASSERT_EQ(result.pulay_coefficients.size(), 1);
    const Complex expected_sos_coefficient = perturbation_matrix_element / denominator;
    const Complex expected_pulay_coefficient = residual_overlap / denominator;
    EXPECT_NEAR(result.sos_coefficients[0].real(), expected_sos_coefficient.real(), 1.0e-14);
    EXPECT_NEAR(result.sos_coefficients[0].imag(), expected_sos_coefficient.imag(), 1.0e-14);
    EXPECT_NEAR(result.pulay_coefficients[0].real(), expected_pulay_coefficient.real(), 1.0e-14);
    EXPECT_NEAR(result.pulay_coefficients[0].imag(), expected_pulay_coefficient.imag(), 1.0e-14);
    EXPECT_NEAR(result.coefficients[0].real(), expected_coefficient.real(), 1.0e-14);
    EXPECT_NEAR(result.coefficients[0].imag(), expected_coefficient.imag(), 1.0e-14);
    expect_vector_near(result.out_wavefunction, out, 1.0e-14);
    Vector expected_in_sos(4, Complex(0.0, 0.0));
    Vector expected_in_pulay(4, Complex(0.0, 0.0));
    for (std::size_t ir = 0; ir != eta.size(); ++ir)
    {
        expected_in_sos[ir] = expected_sos_coefficient * eta[ir];
        expected_in_pulay[ir] = expected_pulay_coefficient * eta[ir];
    }
    expect_vector_near(result.in_sos_wavefunction, expected_in_sos, 1.0e-14);
    expect_vector_near(result.in_pulay_wavefunction, expected_in_pulay, 1.0e-14);
    expect_vector_near(result.reconstructed_wavefunction, standard_delta, 1.0e-14);
    EXPECT_NEAR(result.reconstruction_error, 0.0, 1.0e-14);

    EXPECT_NEAR(std::abs(ModuleRI::sternheimer_fd_grid_dot(occupied, result.out_wavefunction, volume_element)),
                0.0,
                1.0e-14);
    EXPECT_NEAR(std::abs(ModuleRI::sternheimer_fd_grid_dot(eta, result.out_wavefunction, volume_element)),
                0.0,
                1.0e-14);
}

TEST(SternheimerDelta, PulayOperatorTermsCloseToStoredResidualContribution)
{
    using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
    constexpr double volume_element = 0.25;
    Hamiltonian::Grid grid{4, 1, 1, 0.5, 1.0, 1.0, true};
    const std::vector<double> full_local{0.2, -0.1, 0.4, -0.3};
    const std::vector<double> fixed_local{0.1, -0.2, 0.15, -0.05};
    Hamiltonian hamiltonian(grid, full_local, 1.0);
    const Vector occupied(4, Complex(1.0, 0.0));
    const Vector eta{Complex(1.0, 0.0), Complex(0.0, 1.0),
                     Complex(-1.0, 0.0), Complex(0.0, -1.0)};
    const Vector out{Complex(1.0, 0.0), Complex(-1.0, 0.0),
                     Complex(1.0, 0.0), Complex(-1.0, 0.0)};
    const double eta_norm = ModuleRI::sternheimer_fd_grid_norm(eta, volume_element);
    const double occupied_norm = ModuleRI::sternheimer_fd_grid_norm(occupied, volume_element);
    Vector occupied_normalized = occupied;
    Vector eta_normalized = eta;
    for (Complex& value: occupied_normalized)
    {
        value /= occupied_norm;
    }
    for (Complex& value: eta_normalized)
    {
        value /= eta_norm;
    }

    constexpr double virtual_eigenvalue = 0.75;
    Vector residual;
    hamiltonian.apply(eta_normalized, residual);
    for (std::size_t ir = 0; ir != residual.size(); ++ir)
    {
        residual[ir] -= virtual_eigenvalue * eta_normalized[ir];
    }
    const auto dot = [volume_element](const Vector& lhs, const Vector& rhs) {
        return ModuleRI::sternheimer_fd_grid_dot(lhs, rhs, volume_element);
    };
    ModuleRI::SternheimerRPA::project_out_subspace(
        {occupied_normalized, eta_normalized}, dot, residual);

    const std::vector<ModuleRI::SternheimerDeltaVirtualState> states
        = {{eta_normalized, residual, virtual_eigenvalue}};
    const auto components = ModuleRI::decompose_delta_sternheimer_pulay_operator_terms(
        hamiltonian,
        fixed_local,
        {occupied_normalized},
        states,
        out,
        -0.25,
        0.4,
        volume_element);

    for (std::size_t ir = 0; ir != components.total.size(); ++ir)
    {
        const Complex sum = components.kinetic[ir] + components.fixed_local[ir]
                            + components.hxc_local[ir] + components.nonlocal[ir]
                            + components.eigenvalue[ir];
        EXPECT_NEAR(sum.real(), components.total[ir].real(), 1.0e-13);
        EXPECT_NEAR(sum.imag(), components.total[ir].imag(), 1.0e-13);
    }
}

TEST(SternheimerDelta, BuildsDirectSOSWavefunctionFromExplicitVirtualStates)
{
    constexpr double eps_i = -0.5;
    constexpr double omega = 0.25;
    const Vector first = unit_vector(3, 1);
    const Vector second = unit_vector(3, 2);
    const std::vector<ModuleRI::SternheimerDeltaVirtualState> virtual_states
        = {{first, {}, 0.5}, {second, {}, 1.5}};
    const std::vector<Complex> matrix_elements
        = {Complex(0.4, -0.2), Complex(-0.3, 0.5)};

    const Vector response = ModuleRI::build_delta_sternheimer_sos_wavefunction(
        virtual_states, matrix_elements, eps_i, omega);

    const Complex first_coefficient = matrix_elements[0] / Complex(eps_i - 0.5, -omega);
    const Complex second_coefficient = matrix_elements[1] / Complex(eps_i - 1.5, -omega);
    const Vector expected = {Complex(0.0, 0.0), first_coefficient, second_coefficient};
    expect_vector_near(response, expected, 1.0e-14);

    EXPECT_THROW(ModuleRI::build_delta_sternheimer_sos_wavefunction(
                     virtual_states, {matrix_elements[0]}, eps_i, omega),
                 std::invalid_argument);
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

TEST(SternheimerDelta, SolvesGeneralMetricSubspaceEquationByComponents)
{
    const std::vector<Complex> overlap = {Complex(1.0, 0.0), Complex(0.2, -0.1), Complex(0.2, 0.1), Complex(1.3, 0.0)};
    const std::vector<Complex> hamiltonian
        = {Complex(0.8, 0.0), Complex(0.15, -0.05), Complex(0.15, 0.05), Complex(1.4, 0.0)};
    const Complex shift(-0.25, -0.4);

    std::vector<Complex> shifted_operator(4, Complex(0.0, 0.0));
    for (std::size_t index = 0; index != shifted_operator.size(); ++index)
    {
        shifted_operator[index] = hamiltonian[index] - shift * overlap[index];
    }

    const std::vector<Complex> expected_sos = {Complex(0.3, -0.2), Complex(-0.1, 0.4)};
    const std::vector<Complex> expected_pulay = {Complex(-0.05, 0.15), Complex(0.2, -0.1)};
    const std::vector<Complex> rhs = column_major_matrix_vector_product(shifted_operator, expected_sos);
    std::vector<Complex> hamiltonian_out = column_major_matrix_vector_product(shifted_operator, expected_pulay);
    for (Complex& value: hamiltonian_out)
    {
        value = -value;
    }

    const auto result
        = ModuleRI::solve_delta_sternheimer_subspace_coefficients(hamiltonian, overlap, rhs, hamiltonian_out, shift);

    ASSERT_EQ(result.sos.size(), 2);
    ASSERT_EQ(result.pulay.size(), 2);
    ASSERT_EQ(result.total.size(), 2);
    for (std::size_t index = 0; index != 2; ++index)
    {
        EXPECT_NEAR(result.sos[index].real(), expected_sos[index].real(), 1.0e-13);
        EXPECT_NEAR(result.sos[index].imag(), expected_sos[index].imag(), 1.0e-13);
        EXPECT_NEAR(result.pulay[index].real(), expected_pulay[index].real(), 1.0e-13);
        EXPECT_NEAR(result.pulay[index].imag(), expected_pulay[index].imag(), 1.0e-13);
        EXPECT_NEAR(result.total[index].real(), (expected_sos[index] + expected_pulay[index]).real(), 1.0e-13);
        EXPECT_NEAR(result.total[index].imag(), (expected_sos[index] + expected_pulay[index]).imag(), 1.0e-13);
    }
}

TEST(SternheimerDelta, AssemblesReferenceGridHamiltonianWithAnalyticGradientsAndNonlocalProjector)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid{2, 1, 1, 1.0, 1.0, 1.0, false};
    constexpr double volume_element = 0.5;
    using Projector = ModuleRI::SternheimerFDNonlocalProjector;
    Projector::ProjectorBlock block;
    block.projectors = {{Complex(1.0, 0.0), Complex(1.0, 0.0)}};
    block.d_matrix = {{Complex(3.0, 0.0)}};
    const auto nonlocal = std::make_shared<Projector>(grid.size(), volume_element,
                                                      std::vector<Projector::ProjectorBlock>{block});
    ModuleRI::SternheimerFDHamiltonian hamiltonian(grid, {2.0, -1.0}, 0.5, nonlocal);

    ModuleRI::SternheimerDeltaGridFunction basis;
    basis.values = {Complex(1.0, 0.0), Complex(0.0, 1.0)};
    basis.gradients[0] = {Complex(2.0, 0.0), Complex(0.0, 0.0)};
    basis.gradients[1] = {Complex(0.0, 0.0), Complex(1.0, 0.0)};
    basis.gradients[2] = {Complex(0.0, 0.0), Complex(0.0, 0.0)};

    const auto matrices
        = ModuleRI::assemble_delta_sternheimer_grid_matrices(hamiltonian, {basis}, volume_element);

    ASSERT_EQ(matrices.overlap.size(), 1);
    ASSERT_EQ(matrices.hamiltonian.size(), 1);
    EXPECT_NEAR(matrices.overlap[0].real(), 1.0, 1.0e-14);
    EXPECT_NEAR(matrices.overlap[0].imag(), 0.0, 1.0e-14);
    // T=1.25, Vloc=0.5, and Vnl=1.5 in the same grid metric.
    EXPECT_NEAR(matrices.kinetic[0].real(), 1.25, 1.0e-14);
    EXPECT_NEAR(matrices.local_potential[0].real(), 0.5, 1.0e-14);
    EXPECT_NEAR(matrices.nonlocal[0].real(), 1.5, 1.0e-14);
    EXPECT_NEAR(matrices.hamiltonian[0].real(), 3.25, 1.0e-14);
    EXPECT_NEAR(matrices.hamiltonian[0].imag(), 0.0, 1.0e-14);
    EXPECT_NEAR((matrices.kinetic[0] + matrices.local_potential[0] + matrices.nonlocal[0]
                 - matrices.hamiltonian[0])
                    .real(),
                0.0,
                1.0e-14);
}

TEST(SternheimerDelta, CombinesLCAOCoefficientsWithAOValuesAndAnalyticGradients)
{
    ModuleRI::SternheimerDeltaGridFunction first;
    first.values = {Complex(1.0, 0.0), Complex(0.0, 1.0)};
    first.gradients[0] = {Complex(2.0, 0.0), Complex(0.0, 0.0)};
    first.gradients[1] = {Complex(0.0, 1.0), Complex(1.0, 0.0)};
    first.gradients[2] = {Complex(0.0, 0.0), Complex(0.0, -1.0)};

    ModuleRI::SternheimerDeltaGridFunction second;
    second.values = {Complex(0.0, -1.0), Complex(3.0, 0.0)};
    second.gradients[0] = {Complex(1.0, 1.0), Complex(0.0, 2.0)};
    second.gradients[1] = {Complex(-1.0, 0.0), Complex(0.0, 0.0)};
    second.gradients[2] = {Complex(2.0, 0.0), Complex(1.0, 1.0)};

    const Complex first_coefficient(1.0, 0.5);
    const Complex second_coefficient(-0.25, 2.0);
    const auto state = ModuleRI::linear_combination_delta_sternheimer_grid_functions(
        {first, second}, {first_coefficient, second_coefficient});

    auto combine = [first_coefficient, second_coefficient](const Vector& first_values,
                                                            const Vector& second_values) {
        Vector expected(first_values.size(), Complex(0.0, 0.0));
        for (std::size_t ir = 0; ir != expected.size(); ++ir)
        {
            expected[ir] = first_coefficient * first_values[ir]
                           + second_coefficient * second_values[ir];
        }
        return expected;
    };
    expect_vector_near(state.values, combine(first.values, second.values), 1.0e-14);
    for (int direction = 0; direction != 3; ++direction)
    {
        expect_vector_near(state.gradients[static_cast<std::size_t>(direction)],
                           combine(first.gradients[static_cast<std::size_t>(direction)],
                                   second.gradients[static_cast<std::size_t>(direction)]),
                           1.0e-14);
    }

    EXPECT_THROW(ModuleRI::linear_combination_delta_sternheimer_grid_functions({first, second},
                                                                                {first_coefficient}),
                 std::invalid_argument);
}

TEST(SternheimerDelta, AccumulatesBlochImagePhaseIntoValuesAndEveryGradient)
{
    std::vector<ModuleRI::SternheimerDeltaGridFunction> functions(1);
    functions[0].values.assign(2, Complex(0.0, 0.0));
    for (Vector& gradient : functions[0].gradients)
    {
        gradient.assign(2, Complex(0.0, 0.0));
    }

    const std::vector<double> first_values{2.0, -1.0};
    const std::array<std::vector<double>, 3> first_gradients{{
        {3.0, 4.0},
        {-2.0, 1.0},
        {0.5, -0.25},
    }};
    ModuleRI::accumulate_delta_sternheimer_bloch_samples(first_values,
                                                          first_gradients,
                                                          2,
                                                          1,
                                                          0,
                                                          0,
                                                          {0.25, 0.0, 0.0},
                                                          {0, 0, 0},
                                                          functions);

    const std::vector<double> second_values{5.0, 7.0};
    const std::array<std::vector<double>, 3> second_gradients{{
        {1.0, -3.0},
        {2.5, 6.0},
        {-4.0, 8.0},
    }};
    ModuleRI::accumulate_delta_sternheimer_bloch_samples(second_values,
                                                          second_gradients,
                                                          2,
                                                          1,
                                                          0,
                                                          0,
                                                          {0.25, 0.0, 0.0},
                                                          {1, 0, 0},
                                                          functions);

    const Complex imaginary_phase(0.0, 1.0);
    expect_vector_near(functions[0].values,
                       {Complex(2.0, 0.0) + imaginary_phase * 5.0,
                        Complex(-1.0, 0.0) + imaginary_phase * 7.0},
                       1.0e-14);
    for (std::size_t direction = 0; direction != functions[0].gradients.size(); ++direction)
    {
        const Vector expected{
            Complex(first_gradients[direction][0], 0.0)
                + imaginary_phase * second_gradients[direction][0],
            Complex(first_gradients[direction][1], 0.0)
                + imaginary_phase * second_gradients[direction][1],
        };
        expect_vector_near(functions[0].gradients[direction], expected, 1.0e-14);
    }

    std::vector<ModuleRI::SternheimerDeltaGridFunction> gamma_functions(1);
    gamma_functions[0].values.assign(2, Complex(0.0, 0.0));
    for (Vector& gradient : gamma_functions[0].gradients)
    {
        gradient.assign(2, Complex(0.0, 0.0));
    }
    ModuleRI::accumulate_delta_sternheimer_bloch_samples(first_values,
                                                          first_gradients,
                                                          2,
                                                          1,
                                                          0,
                                                          0,
                                                          {0.0, 0.0, 0.0},
                                                          {2, -1, 4},
                                                          gamma_functions);
    ModuleRI::accumulate_delta_sternheimer_bloch_samples(second_values,
                                                          second_gradients,
                                                          2,
                                                          1,
                                                          0,
                                                          0,
                                                          {0.0, 0.0, 0.0},
                                                          {-3, 5, 1},
                                                          gamma_functions);
    expect_vector_near(gamma_functions[0].values, {Complex(7.0, 0.0), Complex(6.0, 0.0)}, 1.0e-14);
    for (std::size_t direction = 0; direction != gamma_functions[0].gradients.size(); ++direction)
    {
        const Vector expected{
            Complex(first_gradients[direction][0] + second_gradients[direction][0], 0.0),
            Complex(first_gradients[direction][1] + second_gradients[direction][1], 0.0),
        };
        expect_vector_near(gamma_functions[0].gradients[direction], expected, 1.0e-14);
    }
}

TEST(SternheimerDelta, AccumulatesLCAOStatesWithoutMaterializingAOGridFunctions)
{
    std::vector<ModuleRI::SternheimerDeltaGridFunction> states(2);
    for (auto& state: states)
    {
        state.values.assign(3, Complex(0.0, 0.0));
        for (auto& gradient: state.gradients)
        {
            gradient.assign(3, Complex(0.0, 0.0));
        }
    }
    const std::vector<double> values{1.0, 3.0, 2.0, 4.0};
    const std::array<std::vector<double>, 3> gradients{{
        {10.0, 30.0, 20.0, 40.0},
        {100.0, 300.0, 200.0, 400.0},
        {1000.0, 3000.0, 2000.0, 4000.0},
    }};
    const std::vector<std::vector<Complex>> coefficients{
        {Complex(1.0, 0.0), Complex(0.0, 0.0)},
        {Complex(0.0, 0.0), Complex(0.0, 2.0)},
    };

    ModuleRI::accumulate_delta_sternheimer_lcao_state_samples(values,
                                                               gradients,
                                                               2,
                                                               2,
                                                               1,
                                                               0,
                                                               coefficients,
                                                               {0.5, 0.0, 0.0},
                                                               {1, 0, 0},
                                                               states);

    EXPECT_NEAR(std::abs(states[0].values[0]), 0.0, 1.0e-14);
    EXPECT_NEAR(states[0].values[1].real(), -1.0, 1.0e-14);
    EXPECT_NEAR(states[0].values[2].real(), -2.0, 1.0e-14);
    EXPECT_NEAR(states[1].values[1].imag(), -6.0, 1.0e-14);
    EXPECT_NEAR(states[1].values[2].imag(), -8.0, 1.0e-14);
    EXPECT_NEAR(states[0].gradients[0][1].real(), -10.0, 1.0e-14);
    EXPECT_NEAR(states[1].gradients[2][2].imag(), -8000.0, 1.0e-12);
}

TEST(SternheimerDelta, OrthonormalizesOccupiedProjectorValuesAndGradientsTogether)
{
    constexpr double volume_element = 1.0;

    ModuleRI::SternheimerDeltaGridFunction first;
    first.values = {Complex(2.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0)};
    first.gradients[0] = {Complex(4.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0)};
    first.gradients[1] = {Complex(0.0, 2.0), Complex(0.0, 0.0), Complex(0.0, 0.0)};
    first.gradients[2] = {Complex(-2.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0)};

    ModuleRI::SternheimerDeltaGridFunction second;
    second.values = {Complex(3.0, 0.0), Complex(4.0, 0.0), Complex(0.0, 0.0)};
    second.gradients[0] = {Complex(11.0, 0.0), Complex(20.0, 0.0), Complex(0.0, 0.0)};
    second.gradients[1] = {Complex(0.0, 7.0), Complex(0.0, 12.0), Complex(0.0, 0.0)};
    second.gradients[2] = {Complex(-5.0, 0.0), Complex(8.0, 0.0), Complex(0.0, 0.0)};

    const auto projector = ModuleRI::orthonormalize_delta_sternheimer_grid_functions(
        {first, second}, volume_element, 1.0e-12);

    ASSERT_EQ(projector.size(), 2);
    expect_vector_near(projector[0].values,
                       {Complex(1.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0)},
                       1.0e-14);
    expect_vector_near(projector[1].values,
                       {Complex(0.0, 0.0), Complex(1.0, 0.0), Complex(0.0, 0.0)},
                       1.0e-14);
    expect_vector_near(projector[1].gradients[0],
                       {Complex(1.25, 0.0), Complex(5.0, 0.0), Complex(0.0, 0.0)},
                       1.0e-14);
    expect_vector_near(projector[1].gradients[1],
                       {Complex(0.0, 1.0), Complex(0.0, 3.0), Complex(0.0, 0.0)},
                       1.0e-14);
    expect_vector_near(projector[1].gradients[2],
                       {Complex(-0.5, 0.0), Complex(2.0, 0.0), Complex(0.0, 0.0)},
                       1.0e-14);

    const Complex overlap00
        = ModuleRI::sternheimer_fd_grid_dot(projector[0].values, projector[0].values, volume_element);
    const Complex overlap01
        = ModuleRI::sternheimer_fd_grid_dot(projector[0].values, projector[1].values, volume_element);
    const Complex overlap11
        = ModuleRI::sternheimer_fd_grid_dot(projector[1].values, projector[1].values, volume_element);
    EXPECT_NEAR(overlap00.real(), 1.0, 1.0e-14);
    EXPECT_NEAR(std::abs(overlap01), 0.0, 1.0e-14);
    EXPECT_NEAR(overlap11.real(), 1.0, 1.0e-14);
}

TEST(SternheimerDelta, ReferenceSubspaceTransformsValuesAndGradientsTogether)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid{2, 1, 1, 1.0, 1.0, 1.0, false};
    constexpr double volume_element = 1.0;
    ModuleRI::SternheimerFDHamiltonian hamiltonian(grid, {0.0, 0.0}, 0.5);

    ModuleRI::SternheimerDeltaGridFunction occupied;
    occupied.values = {Complex(1.0, 0.0), Complex(0.0, 0.0)};
    occupied.gradients[0] = {Complex(1.0, 0.0), Complex(0.0, 0.0)};
    occupied.gradients[1] = Vector(2, Complex(0.0, 0.0));
    occupied.gradients[2] = Vector(2, Complex(0.0, 0.0));

    ModuleRI::SternheimerDeltaGridFunction candidate;
    candidate.values = {Complex(1.0, 0.0), Complex(1.0, 0.0)};
    candidate.gradients[0] = {Complex(2.0, 0.0), Complex(3.0, 0.0)};
    candidate.gradients[1] = Vector(2, Complex(0.0, 0.0));
    candidate.gradients[2] = Vector(2, Complex(0.0, 0.0));

    ModuleRI::SternheimerDeltaSubspaceOptions options;
    options.max_virtual_states = 1;
    const auto subspace = ModuleRI::build_reference_delta_sternheimer_subspace(
        hamiltonian, {occupied}, {candidate}, volume_element, options);

    ASSERT_EQ(subspace.virtual_states.size(), 1);
    ASSERT_EQ(subspace.grid_functions.size(), 1);
    expect_vector_near(subspace.grid_functions[0].values,
                       {Complex(0.0, 0.0), Complex(1.0, 0.0)},
                       1.0e-14);
    expect_vector_near(subspace.grid_functions[0].gradients[0],
                       {Complex(1.0, 0.0), Complex(3.0, 0.0)},
                       1.0e-14);
    EXPECT_NEAR(subspace.virtual_states[0].eigenvalue, 5.0, 1.0e-13);
    Vector h_virtual;
    hamiltonian.apply(subspace.virtual_states[0].orbital, h_virtual);
    const double grid_diagonal = std::real(ModuleRI::sternheimer_fd_grid_dot(
        subspace.virtual_states[0].orbital, h_virtual, volume_element));
    const double expected_difference = std::abs(grid_diagonal - subspace.virtual_states[0].eigenvalue);
    EXPECT_NEAR(subspace.full_grid_hamiltonian_max_abs_difference, expected_difference, 1.0e-13);
    EXPECT_NEAR(subspace.full_grid_hamiltonian_relative_difference,
                expected_difference / std::abs(subspace.virtual_states[0].eigenvalue),
                1.0e-13);
}

TEST(SternheimerDelta, ProductionOptionsDropGridDiagnosticsWithoutChangingVirtualState)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid{2, 1, 1, 1.0, 1.0, 1.0, false};
    constexpr double volume_element = 1.0;
    ModuleRI::SternheimerFDHamiltonian hamiltonian(grid, {0.0, 0.0}, 0.5);

    ModuleRI::SternheimerDeltaGridFunction occupied;
    occupied.values = {Complex(1.0, 0.0), Complex(0.0, 0.0)};
    occupied.gradients[0] = occupied.values;
    occupied.gradients[1] = Vector(2, Complex(0.0, 0.0));
    occupied.gradients[2] = Vector(2, Complex(0.0, 0.0));

    ModuleRI::SternheimerDeltaGridFunction candidate;
    candidate.values = {Complex(0.0, 0.0), Complex(1.0, 0.0)};
    candidate.gradients[0] = {Complex(0.0, 0.0), Complex(3.0, 0.0)};
    candidate.gradients[1] = Vector(2, Complex(0.0, 0.0));
    candidate.gradients[2] = Vector(2, Complex(0.0, 0.0));

    ModuleRI::SternheimerDeltaSubspaceOptions reference_options;
    reference_options.max_virtual_states = 1;
    const auto reference = ModuleRI::build_reference_delta_sternheimer_subspace(
        hamiltonian, {occupied}, {candidate}, volume_element, reference_options);

    ModuleRI::SternheimerDeltaSubspaceOptions production_options = reference_options;
    production_options.retain_grid_functions = false;
    production_options.evaluate_full_grid_difference = false;
    const auto production = ModuleRI::build_reference_delta_sternheimer_subspace(
        hamiltonian, {occupied}, {candidate}, volume_element, production_options);

    ASSERT_EQ(production.virtual_states.size(), 1);
    EXPECT_TRUE(production.grid_functions.empty());
    expect_vector_near(production.virtual_states[0].orbital,
                       reference.virtual_states[0].orbital,
                       1.0e-14);
    expect_vector_near(production.virtual_states[0].residual,
                       reference.virtual_states[0].residual,
                       1.0e-14);
    EXPECT_NEAR(production.virtual_states[0].eigenvalue,
                reference.virtual_states[0].eigenvalue,
                1.0e-14);
    EXPECT_EQ(production.full_grid_hamiltonian_relative_difference, -1.0);
    EXPECT_EQ(production.full_grid_hamiltonian_max_abs_difference, -1.0);
}

TEST(SternheimerDelta, ReferenceSubspacePivotsPastNearDependentCandidates)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid{4, 1, 1, 1.0, 1.0, 1.0, false};
    constexpr double volume_element = 1.0;
    ModuleRI::SternheimerFDHamiltonian hamiltonian(grid, {0.0, 1.0, 2.0, 3.0}, 0.0);

    auto grid_function = [](const Vector& values) {
        ModuleRI::SternheimerDeltaGridFunction function;
        function.values = values;
        for (Vector& gradient: function.gradients)
        {
            gradient.assign(values.size(), Complex(0.0, 0.0));
        }
        return function;
    };

    const auto occupied = grid_function(unit_vector(4, 0));
    Vector nearly_occupied = unit_vector(4, 0);
    nearly_occupied[1] = Complex(1.0e-8, 0.0);
    const std::vector<ModuleRI::SternheimerDeltaGridFunction> candidates{
        grid_function(nearly_occupied),
        grid_function(unit_vector(4, 2)),
        grid_function(unit_vector(4, 3)),
    };

    ModuleRI::SternheimerDeltaSubspaceOptions options;
    options.max_virtual_states = 2;
    options.norm_tolerance = 1.0e-10;
    const auto subspace = ModuleRI::build_reference_delta_sternheimer_subspace(
        hamiltonian, {occupied}, candidates, volume_element, options);

    ASSERT_EQ(subspace.virtual_states.size(), 2);
    EXPECT_EQ(subspace.accepted_candidates, 2);
    EXPECT_EQ(subspace.discarded_candidates, 1);
    EXPECT_NEAR(subspace.virtual_states[0].eigenvalue, 2.0, 1.0e-13);
    EXPECT_NEAR(subspace.virtual_states[1].eigenvalue, 3.0, 1.0e-13);
}

TEST(SternheimerDelta, BlockCompleteReferenceSubspaceMatchesLegacyPath)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid{5, 1, 1, 1.0, 1.0, 1.0, false};
    constexpr double volume_element = 1.0;
    ModuleRI::SternheimerFDHamiltonian hamiltonian(grid, {0.0, 0.7, 1.9, 3.2, 5.1}, 0.4);

    auto grid_function = [](const Vector& values, const double gradient_scale) {
        ModuleRI::SternheimerDeltaGridFunction function;
        function.values = values;
        for (std::size_t direction = 0; direction != function.gradients.size(); ++direction)
        {
            function.gradients[direction].resize(values.size());
            for (std::size_t ir = 0; ir != values.size(); ++ir)
            {
                const double factor = gradient_scale * static_cast<double>((direction + 1) * (ir + 1));
                function.gradients[direction][ir] = factor * values[ir];
            }
        }
        return function;
    };

    const auto occupied = grid_function(unit_vector(5, 0), 0.03);
    Vector candidate0 = unit_vector(5, 1);
    candidate0[0] = Complex(0.2, -0.1);
    Vector candidate1 = unit_vector(5, 2);
    candidate1[1] = Complex(0.3, 0.2);
    Vector candidate2 = unit_vector(5, 4);
    candidate2[2] = Complex(-0.25, 0.15);
    candidate2[3] = Complex(0.4, -0.1);
    const std::vector<ModuleRI::SternheimerDeltaGridFunction> candidates{
        grid_function(candidate0, 0.05),
        grid_function(candidate1, 0.07),
        grid_function(candidate2, 0.09),
    };

    ModuleRI::SternheimerDeltaSubspaceOptions block_options;
    block_options.max_virtual_states = static_cast<int>(candidates.size());
    block_options.evaluate_full_grid_difference = false;
    ModuleRI::SternheimerDeltaSubspaceOptions legacy_options = block_options;
    legacy_options.use_block_generalized_eigensolver = false;

    const auto block = ModuleRI::build_reference_delta_sternheimer_subspace(
        hamiltonian, {occupied}, candidates, volume_element, block_options);
    const auto legacy = ModuleRI::build_reference_delta_sternheimer_subspace(
        hamiltonian, {occupied}, candidates, volume_element, legacy_options);

    ASSERT_EQ(block.virtual_states.size(), legacy.virtual_states.size());
    EXPECT_EQ(block.accepted_candidates, legacy.accepted_candidates);
    EXPECT_EQ(block.discarded_candidates, legacy.discarded_candidates);
    for (std::size_t istate = 0; istate != block.virtual_states.size(); ++istate)
    {
        EXPECT_NEAR(block.virtual_states[istate].eigenvalue,
                    legacy.virtual_states[istate].eigenvalue,
                    1.0e-11);
        const Complex overlap = ModuleRI::sternheimer_fd_grid_dot(
            legacy.virtual_states[istate].orbital,
            block.virtual_states[istate].orbital,
            volume_element);
        EXPECT_NEAR(std::abs(overlap), 1.0, 1.0e-11);
        const Complex residual_overlap = ModuleRI::sternheimer_fd_grid_dot(
            legacy.virtual_states[istate].residual,
            block.virtual_states[istate].residual,
            volume_element);
        const double legacy_residual_norm = ModuleRI::sternheimer_fd_grid_norm(
            legacy.virtual_states[istate].residual, volume_element);
        const double block_residual_norm = ModuleRI::sternheimer_fd_grid_norm(
            block.virtual_states[istate].residual, volume_element);
        EXPECT_NEAR(block_residual_norm, legacy_residual_norm, 1.0e-11);
        EXPECT_NEAR(std::abs(residual_overlap),
                    block_residual_norm * legacy_residual_norm,
                    1.0e-11);
    }
}

TEST(SternheimerDelta, CapsRequestedVirtualStatesAtTheLCAOComplementDimension)
{
    EXPECT_EQ(ModuleRI::sternheimer_delta_virtual_state_limit(0, 176, 16), 160);
    EXPECT_EQ(ModuleRI::sternheimer_delta_virtual_state_limit(40, 176, 16), 40);
    EXPECT_EQ(ModuleRI::sternheimer_delta_virtual_state_limit(200, 176, 16), 160);
    EXPECT_EQ(ModuleRI::sternheimer_delta_virtual_state_limit(0, 16, 16), 0);
    EXPECT_THROW(ModuleRI::sternheimer_delta_virtual_state_limit(-1, 176, 16), std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_delta_virtual_state_limit(0, 15, 16), std::invalid_argument);
}

TEST(SternheimerDelta, BuildsPeriodicCenteredGradientForOccupiedGridState)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid{4, 1, 1, 0.5, 1.0, 1.0, true};
    Vector values(static_cast<std::size_t>(grid.size()), Complex(0.0, 0.0));
    for (int ix = 0; ix != grid.nx; ++ix)
    {
        const double phase = 0.5 * std::acos(-1.0) * static_cast<double>(ix);
        values[static_cast<std::size_t>(ix)] = std::exp(Complex(0.0, phase));
    }

    const auto function = ModuleRI::make_delta_sternheimer_grid_function_with_fd_gradients(values, grid);

    expect_vector_near(function.values, values, 1.0e-14);
    for (int ix = 0; ix != grid.nx; ++ix)
    {
        const Complex expected = Complex(0.0, 2.0) * values[static_cast<std::size_t>(ix)];
        EXPECT_NEAR(function.gradients[0][static_cast<std::size_t>(ix)].real(), expected.real(), 1.0e-14);
        EXPECT_NEAR(function.gradients[0][static_cast<std::size_t>(ix)].imag(), expected.imag(), 1.0e-14);
        EXPECT_NEAR(std::abs(function.gradients[1][static_cast<std::size_t>(ix)]), 0.0, 1.0e-14);
        EXPECT_NEAR(std::abs(function.gradients[2][static_cast<std::size_t>(ix)]), 0.0, 1.0e-14);
    }
}

TEST(SternheimerDelta, EnumeratesEveryPeriodicOrbitalImageInsideTheCutoff)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid{12, 12, 12, 1.0, 1.0, 1.0, true};
    const auto centered_wide
        = ModuleRI::enumerate_delta_sternheimer_periodic_images(grid, {6.0, 6.0, 6.0}, 8.0);
    EXPECT_EQ(centered_wide.size(), 27);

    const auto centered_local
        = ModuleRI::enumerate_delta_sternheimer_periodic_images(grid, {6.0, 6.0, 6.0}, 5.0);
    ASSERT_EQ(centered_local.size(), 1);
    EXPECT_EQ(centered_local[0], (std::array<int, 3>{0, 0, 0}));

    const auto boundary_crossing
        = ModuleRI::enumerate_delta_sternheimer_periodic_images(grid, {1.0, 6.0, 6.0}, 2.0);
    ASSERT_EQ(boundary_crossing.size(), 2);
    EXPECT_NE(std::find(boundary_crossing.begin(), boundary_crossing.end(), std::array<int, 3>{0, 0, 0}),
              boundary_crossing.end());
    EXPECT_NE(std::find(boundary_crossing.begin(), boundary_crossing.end(), std::array<int, 3>{1, 0, 0}),
              boundary_crossing.end());

    grid.periodic = false;
    const auto nonperiodic
        = ModuleRI::enumerate_delta_sternheimer_periodic_images(grid, {6.0, 6.0, 6.0}, 8.0);
    ASSERT_EQ(nonperiodic.size(), 1);
    EXPECT_EQ(nonperiodic[0], (std::array<int, 3>{0, 0, 0}));
}

TEST(SternheimerDelta, EnumeratesPeriodicImagesUsingNonorthogonalDualVectors)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid;
    grid.nx = 4;
    grid.ny = 4;
    grid.nz = 4;
    grid.hx = std::sqrt(2.0) / 4.0;
    grid.hy = std::sqrt(2.0) / 4.0;
    grid.hz = std::sqrt(2.0) / 4.0;
    grid.periodic = true;
    grid.lattice_vectors = {{{0.0, 1.0, 1.0},
                             {1.0, 0.0, 1.0},
                             {1.0, 1.0, 0.0}}};

    const auto images
        = ModuleRI::enumerate_delta_sternheimer_periodic_images(grid, {0.5, 0.5, 0.5}, 0.9);

    EXPECT_EQ(images.size(), 27);
    EXPECT_NE(std::find(images.begin(), images.end(), std::array<int, 3>{-1, 0, 0}), images.end());
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

TEST(SternheimerDelta, ParsesDiagnosticABlockModes)
{
    using Mode = ModuleRI::SternheimerDeltaABlockMode;

    EXPECT_EQ(ModuleRI::parse_sternheimer_delta_a_block_mode("reference_value_gradient"),
              Mode::ReferenceValueGradient);
    EXPECT_EQ(ModuleRI::parse_sternheimer_delta_a_block_mode("grid"), Mode::FullGrid);
    EXPECT_STREQ(ModuleRI::sternheimer_delta_a_block_mode_name(Mode::ReferenceValueGradient),
                 "reference_value_gradient");
    EXPECT_STREQ(ModuleRI::sternheimer_delta_a_block_mode_name(Mode::FullGrid), "grid");
    EXPECT_THROW(ModuleRI::parse_sternheimer_delta_a_block_mode("mixed"), std::invalid_argument);
}

TEST(SternheimerDelta, FullGridABlockModeMatchesStandardForGeneralProjectedSubspace)
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
    const auto occupied_function
        = ModuleRI::make_delta_sternheimer_grid_function_with_fd_gradients(states.wavefunctions[0], grid);
    const auto candidate_function
        = ModuleRI::make_delta_sternheimer_grid_function_with_fd_gradients(candidate, grid);
    ModuleRI::SternheimerDeltaSubspaceOptions subspace_options;
    subspace_options.max_virtual_states = 1;
    const auto subspace = ModuleRI::build_delta_sternheimer_subspace_by_mode(
        hamiltonian,
        {occupied_function},
        {candidate_function},
        volume_element,
        subspace_options,
        ModuleRI::SternheimerDeltaABlockMode::FullGrid);
    ASSERT_EQ(subspace.virtual_states.size(), 1);

    const std::vector<double> perturbation = {0.7, -0.3, 0.2, 0.5, -0.4};
    Vector rhs;
    ModuleRI::SternheimerRPA::build_rhs_from_hartree_perturbation(perturbation,
                                                                  states.wavefunctions[0],
                                                                  rhs);
    const auto matrix_elements = ModuleRI::delta_sternheimer_perturbation_matrix_elements(
        subspace.virtual_states, perturbation, states.wavefunctions[0], volume_element);
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
    const auto delta = ModuleRI::solve_delta_sternheimer_linear_response(hamiltonian,
                                                                         {states.wavefunctions[0]},
                                                                         states.eigenvalues[0],
                                                                         rhs,
                                                                         subspace.virtual_states,
                                                                         matrix_elements,
                                                                         omega,
                                                                         volume_element,
                                                                         solver_options);

    ASSERT_TRUE(standard.solver.converged);
    ASSERT_TRUE(delta.solver.converged);
    EXPECT_LT(delta.residual_norm, 1.0e-9);
    expect_vector_near(delta.response.reconstructed_wavefunction, standard.delta_wavefunction, 1.0e-8);
}

TEST(SternheimerDelta, FullGridABlockRetainsComplexHermitianCouplingsAtTwistedKPoint)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid{5, 1, 1, 0.6, 1.0, 1.0, true};
    grid.kpoint = {0.19, 0.0, 0.0};
    constexpr double volume_element = 0.6;
    ModuleRI::SternheimerFDHamiltonian hamiltonian(grid, {0.10, -0.05, 0.20, -0.10, 0.05});
    const auto states = ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 5, volume_element);

    const double inverse_sqrt_two = 1.0 / std::sqrt(2.0);
    Vector first_candidate(states.wavefunctions[1].size(), Complex(0.0, 0.0));
    Vector second_candidate(states.wavefunctions[1].size(), Complex(0.0, 0.0));
    for (std::size_t ir = 0; ir != first_candidate.size(); ++ir)
    {
        first_candidate[ir]
            = inverse_sqrt_two * (states.wavefunctions[1][ir] + Complex(0.0, 1.0) * states.wavefunctions[2][ir]);
        second_candidate[ir]
            = inverse_sqrt_two * (Complex(0.0, 1.0) * states.wavefunctions[1][ir] + states.wavefunctions[2][ir]);
    }
    const auto occupied_function
        = ModuleRI::make_delta_sternheimer_grid_function_with_fd_gradients(states.wavefunctions[0], grid);
    const std::vector<ModuleRI::SternheimerDeltaGridFunction> candidate_functions{
        ModuleRI::make_delta_sternheimer_grid_function_with_fd_gradients(first_candidate, grid),
        ModuleRI::make_delta_sternheimer_grid_function_with_fd_gradients(second_candidate, grid),
    };
    const auto subspace = ModuleRI::build_delta_sternheimer_subspace_by_mode(
        hamiltonian,
        {occupied_function},
        candidate_functions,
        volume_element,
        ModuleRI::SternheimerDeltaSubspaceOptions(),
        ModuleRI::SternheimerDeltaABlockMode::FullGrid);
    ASSERT_EQ(subspace.virtual_states.size(), 2);
    EXPECT_LT(subspace.full_grid_hamiltonian_relative_difference, 1.0e-10);

    const Vector perturbation = {Complex(0.70, 0.11), Complex(-0.30, 0.05), Complex(0.21, -0.08),
                                 Complex(0.48, 0.03), Complex(-0.37, 0.09)};
    Vector rhs;
    ModuleRI::SternheimerRPA::build_rhs_from_hartree_perturbation(perturbation,
                                                                  states.wavefunctions[0],
                                                                  rhs);
    const auto matrix_elements = ModuleRI::delta_sternheimer_perturbation_matrix_elements(
        subspace.virtual_states, perturbation, states.wavefunctions[0], volume_element);
    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 80;
    options.residual_tol = 1.0e-12;
    constexpr double omega = 0.45;
    const auto standard = ModuleRI::solve_sternheimer_fd_linear_response(hamiltonian,
                                                                         {states.wavefunctions[0]},
                                                                         states.eigenvalues[0],
                                                                         rhs,
                                                                         omega,
                                                                         volume_element,
                                                                         options);
    const auto delta = ModuleRI::solve_delta_sternheimer_linear_response(hamiltonian,
                                                                         {states.wavefunctions[0]},
                                                                         states.eigenvalues[0],
                                                                         rhs,
                                                                         subspace.virtual_states,
                                                                         matrix_elements,
                                                                         omega,
                                                                         volume_element,
                                                                         options);

    ASSERT_TRUE(standard.solver.converged);
    ASSERT_TRUE(delta.solver.converged);
    expect_vector_near(delta.response.reconstructed_wavefunction, standard.delta_wavefunction, 1.0e-9);
}

TEST(SternheimerDelta, EmptyABlockReducesToStandardSternheimer)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid{5, 1, 1, 0.6, 1.0, 1.0, true};
    constexpr double volume_element = 0.6;
    ModuleRI::SternheimerFDHamiltonian hamiltonian(grid, {0.10, -0.05, 0.20, -0.10, 0.05});
    const auto states = ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 5, volume_element);
    const std::vector<double> perturbation = {0.7, -0.3, 0.2, 0.5, -0.4};
    Vector rhs;
    ModuleRI::SternheimerRPA::build_rhs_from_hartree_perturbation(perturbation,
                                                                  states.wavefunctions[0],
                                                                  rhs);
    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 80;
    options.residual_tol = 1.0e-12;
    constexpr double omega = 0.45;
    const auto standard = ModuleRI::solve_sternheimer_fd_linear_response(hamiltonian,
                                                                         {states.wavefunctions[0]},
                                                                         states.eigenvalues[0],
                                                                         rhs,
                                                                         omega,
                                                                         volume_element,
                                                                         options);
    const auto delta = ModuleRI::solve_delta_sternheimer_linear_response(hamiltonian,
                                                                         {states.wavefunctions[0]},
                                                                         states.eigenvalues[0],
                                                                         rhs,
                                                                         {},
                                                                         {},
                                                                         omega,
                                                                         volume_element,
                                                                         options);

    ASSERT_TRUE(standard.solver.converged);
    ASSERT_TRUE(delta.solver.converged);
    EXPECT_TRUE(delta.response.coefficients.empty());
    expect_vector_near(delta.response.reconstructed_wavefunction, standard.delta_wavefunction, 1.0e-10);
}

TEST(SternheimerDelta, CompleteFullGridABlockReducesToSOSWithZeroOutGridResponse)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid{5, 1, 1, 0.6, 1.0, 1.0, true};
    grid.kpoint = {0.19, 0.0, 0.0};
    constexpr double volume_element = 0.6;
    ModuleRI::SternheimerFDHamiltonian hamiltonian(grid, {0.10, -0.05, 0.20, -0.10, 0.05});
    const auto states = ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 5, volume_element);
    const auto occupied_function
        = ModuleRI::make_delta_sternheimer_grid_function_with_fd_gradients(states.wavefunctions[0], grid);
    std::vector<ModuleRI::SternheimerDeltaGridFunction> candidate_functions;
    for (int state_index = 1; state_index != 5; ++state_index)
    {
        candidate_functions.push_back(ModuleRI::make_delta_sternheimer_grid_function_with_fd_gradients(
            states.wavefunctions[static_cast<std::size_t>(state_index)], grid));
    }
    ModuleRI::SternheimerDeltaSubspaceOptions subspace_options;
    const auto subspace = ModuleRI::build_delta_sternheimer_subspace_by_mode(
        hamiltonian,
        {occupied_function},
        candidate_functions,
        volume_element,
        subspace_options,
        ModuleRI::SternheimerDeltaABlockMode::FullGrid);
    ASSERT_EQ(subspace.virtual_states.size(), 4);

    const Vector perturbation = {Complex(0.70, 0.11), Complex(-0.30, 0.05), Complex(0.21, -0.08),
                                 Complex(0.48, 0.03), Complex(-0.37, 0.09)};
    Vector rhs;
    ModuleRI::SternheimerRPA::build_rhs_from_hartree_perturbation(perturbation,
                                                                  states.wavefunctions[0],
                                                                  rhs);
    const auto matrix_elements = ModuleRI::delta_sternheimer_perturbation_matrix_elements(
        subspace.virtual_states, perturbation, states.wavefunctions[0], volume_element);
    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 80;
    options.residual_tol = 1.0e-12;
    constexpr double omega = 0.45;
    const auto standard = ModuleRI::solve_sternheimer_fd_linear_response(hamiltonian,
                                                                         {states.wavefunctions[0]},
                                                                         states.eigenvalues[0],
                                                                         rhs,
                                                                         omega,
                                                                         volume_element,
                                                                         options);
    const auto delta = ModuleRI::solve_delta_sternheimer_linear_response(hamiltonian,
                                                                         {states.wavefunctions[0]},
                                                                         states.eigenvalues[0],
                                                                         rhs,
                                                                         subspace.virtual_states,
                                                                         matrix_elements,
                                                                         omega,
                                                                         volume_element,
                                                                         options);
    const Vector sos = ModuleRI::build_sternheimer_fd_complete_sos_response(states,
                                                                            1,
                                                                            0,
                                                                            rhs,
                                                                            omega,
                                                                            volume_element);

    ASSERT_TRUE(standard.solver.converged);
    ASSERT_TRUE(delta.solver.converged);
    EXPECT_LT(delta.response.out_norm, 1.0e-10);
    expect_vector_near(delta.response.reconstructed_wavefunction, standard.delta_wavefunction, 1.0e-9);
    expect_vector_near(delta.response.reconstructed_wavefunction, sos, 1.0e-9);
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
    ASSERT_EQ(delta.response.sos_coefficients.size(), subspace.virtual_states.size());
    ASSERT_EQ(delta.response.pulay_coefficients.size(), subspace.virtual_states.size());
    ASSERT_EQ(delta.response.coefficients.size(), subspace.virtual_states.size());
    for (std::size_t ia = 0; ia != subspace.virtual_states.size(); ++ia)
    {
        const Complex component_sum = delta.response.sos_coefficients[ia] + delta.response.pulay_coefficients[ia];
        EXPECT_NEAR(delta.response.coefficients[ia].real(), component_sum.real(), 1.0e-12);
        EXPECT_NEAR(delta.response.coefficients[ia].imag(), component_sum.imag(), 1.0e-12);
    }
    Vector reconstructed_from_components = delta.response.out_wavefunction;
    for (std::size_t ir = 0; ir != reconstructed_from_components.size(); ++ir)
    {
        reconstructed_from_components[ir]
            += delta.response.in_sos_wavefunction[ir] + delta.response.in_pulay_wavefunction[ir];
    }
    expect_vector_near(delta.response.reconstructed_wavefunction, reconstructed_from_components, 1.0e-12);
    expect_vector_near(delta.response.reconstructed_wavefunction, standard.delta_wavefunction, 1.0e-8);
}

TEST(SternheimerDelta, SharedFixedSubspaceMatchesCompatibilitySolver)
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
    const std::vector<Complex> perturbation_matrix_elements
        = ModuleRI::delta_sternheimer_perturbation_matrix_elements(subspace.virtual_states,
                                                                   perturbation,
                                                                   states.wavefunctions[0],
                                                                   volume_element);
    ModuleRI::SternheimerRPA::SolverOptions solver_options;
    solver_options.max_iter = 80;
    solver_options.residual_tol = 1.0e-12;
    constexpr double omega = 0.45;

    const auto compatibility = ModuleRI::solve_delta_sternheimer_linear_response(hamiltonian,
                                                                                  {states.wavefunctions[0]},
                                                                                  states.eigenvalues[0],
                                                                                  rhs,
                                                                                  subspace.virtual_states,
                                                                                  perturbation_matrix_elements,
                                                                                  omega,
                                                                                  volume_element,
                                                                                  solver_options);
    const ModuleRI::SternheimerDeltaFixedSubspace fixed_subspace
        = ModuleRI::build_delta_sternheimer_fixed_subspace({states.wavefunctions[0]},
                                                            subspace.virtual_states);
    const auto shared = ModuleRI::solve_delta_sternheimer_linear_response(hamiltonian,
                                                                           fixed_subspace,
                                                                           states.eigenvalues[0],
                                                                           rhs,
                                                                           subspace.virtual_states,
                                                                           perturbation_matrix_elements,
                                                                           omega,
                                                                           volume_element,
                                                                           solver_options);

    EXPECT_EQ(fixed_subspace.functions.size(), 2U);
    EXPECT_EQ(shared.solver.converged, compatibility.solver.converged);
    EXPECT_EQ(shared.solver.iterations, compatibility.solver.iterations);
    EXPECT_NEAR(shared.solver.relative_residual, compatibility.solver.relative_residual, 1.0e-15);
    EXPECT_NEAR(shared.residual_norm, compatibility.residual_norm, 1.0e-15);
    expect_vector_near(shared.response.out_wavefunction, compatibility.response.out_wavefunction, 1.0e-14);
    expect_vector_near(shared.response.in_sos_wavefunction, compatibility.response.in_sos_wavefunction, 1.0e-14);
    expect_vector_near(shared.response.in_pulay_wavefunction, compatibility.response.in_pulay_wavefunction, 1.0e-14);
    expect_vector_near(shared.response.reconstructed_wavefunction,
                       compatibility.response.reconstructed_wavefunction,
                       1.0e-14);
}

TEST(SternheimerDelta, BatchLinearResponseMatchesIndependentFD8Channels)
{
    using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
    using Matrix = Hamiltonian::Matrix;
    Hamiltonian::Grid grid{9, 1, 1, 0.55, 1.0, 1.0, true};
    constexpr double volume_element = 0.55;
    std::vector<double> local_potential(static_cast<std::size_t>(grid.size()));
    for (int ir = 0; ir != grid.size(); ++ir)
    {
        local_potential[static_cast<std::size_t>(ir)] = 0.08 * (ir % 4) - 0.03 * (ir % 7);
    }
    const Hamiltonian hamiltonian(grid, local_potential, 1.0, nullptr, 8);
    const auto states = ModuleRI::solve_sternheimer_fd_zero_order_dense(
        hamiltonian, grid.size(), volume_element);

    Vector candidate = states.wavefunctions[2];
    for (std::size_t ir = 0; ir != candidate.size(); ++ir)
    {
        candidate[ir] += Complex(0.2, -0.15) * states.wavefunctions[3][ir];
    }
    ModuleRI::SternheimerDeltaSubspaceOptions subspace_options;
    subspace_options.max_virtual_states = 1;
    const auto subspace = ModuleRI::build_delta_sternheimer_subspace(hamiltonian,
                                                                     {states.wavefunctions[0]},
                                                                     {candidate},
                                                                     volume_element,
                                                                     subspace_options);
    ASSERT_EQ(subspace.virtual_states.size(), 1U);
    const auto fixed_subspace = ModuleRI::build_delta_sternheimer_fixed_subspace(
        {states.wavefunctions[0]}, subspace.virtual_states);

    const std::vector<std::vector<double>> perturbations
        = {{0.7, -0.3, 0.2, 0.5, -0.4, 0.1, 0.6, -0.2, 0.35},
           {-0.2, 0.4, 0.8, -0.1, 0.3, -0.6, 0.25, 0.15, -0.45}};
    Matrix rhs(perturbations.size());
    std::vector<std::vector<Complex>> perturbation_matrix_elements(perturbations.size());
    for (std::size_t column = 0; column != perturbations.size(); ++column)
    {
        ModuleRI::SternheimerRPA::build_rhs_from_hartree_perturbation(
            perturbations[column], states.wavefunctions[0], rhs[column]);
        perturbation_matrix_elements[column]
            = ModuleRI::delta_sternheimer_perturbation_matrix_elements(subspace.virtual_states,
                                                                       perturbations[column],
                                                                       states.wavefunctions[0],
                                                                       volume_element);
    }

    ModuleRI::SternheimerRPA::SolverOptions solver_options;
    solver_options.max_iter = 120;
    solver_options.residual_tol = 1.0e-11;
    constexpr double omega = 0.45;
    std::vector<ModuleRI::SternheimerDeltaLinearResponse> scalar;
    for (std::size_t column = 0; column != rhs.size(); ++column)
    {
        scalar.push_back(ModuleRI::solve_delta_sternheimer_linear_response(hamiltonian,
                                                                           fixed_subspace,
                                                                           states.eigenvalues[0],
                                                                           rhs[column],
                                                                           subspace.virtual_states,
                                                                           perturbation_matrix_elements[column],
                                                                           omega,
                                                                           volume_element,
                                                                           solver_options));
    }

    const auto batch = ModuleRI::solve_delta_sternheimer_linear_response_batch(hamiltonian,
                                                                                fixed_subspace,
                                                                                states.eigenvalues[0],
                                                                                rhs,
                                                                                subspace.virtual_states,
                                                                                perturbation_matrix_elements,
                                                                                omega,
                                                                                volume_element,
                                                                                solver_options);

    ASSERT_EQ(batch.size(), scalar.size());
    for (std::size_t column = 0; column != scalar.size(); ++column)
    {
        EXPECT_EQ(batch[column].solver.converged, scalar[column].solver.converged);
        EXPECT_EQ(batch[column].solver.iterations, scalar[column].solver.iterations);
        EXPECT_NEAR(batch[column].solver.absolute_residual,
                    scalar[column].solver.absolute_residual,
                    1.0e-12);
        EXPECT_NEAR(batch[column].solver.relative_residual,
                    scalar[column].solver.relative_residual,
                    1.0e-12);
        EXPECT_NEAR(batch[column].residual_norm, scalar[column].residual_norm, 1.0e-11);
        expect_vector_near(batch[column].response.out_wavefunction,
                           scalar[column].response.out_wavefunction,
                           1.0e-11);
        expect_vector_near(batch[column].response.in_sos_wavefunction,
                           scalar[column].response.in_sos_wavefunction,
                           1.0e-11);
        expect_vector_near(batch[column].response.in_pulay_wavefunction,
                           scalar[column].response.in_pulay_wavefunction,
                           1.0e-11);
        expect_vector_near(batch[column].response.reconstructed_wavefunction,
                           scalar[column].response.reconstructed_wavefunction,
                           1.0e-11);
        ASSERT_EQ(batch[column].response.coefficients.size(), scalar[column].response.coefficients.size());
        for (std::size_t ia = 0; ia != scalar[column].response.coefficients.size(); ++ia)
        {
            EXPECT_NEAR(batch[column].response.coefficients[ia].real(),
                        scalar[column].response.coefficients[ia].real(),
                        1.0e-11);
            EXPECT_NEAR(batch[column].response.coefficients[ia].imag(),
                        scalar[column].response.coefficients[ia].imag(),
                        1.0e-11);
        }
    }
}

TEST(SternheimerDelta, StrictIndependentDeltaHamiltonianMatchesExplicitHybridBlockEquation)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid{5, 1, 1, 0.6, 1.0, 1.0, true};
    constexpr double volume_element = 0.6;
    const std::vector<double> potential = {0.10, -0.05, 0.20, -0.10, 0.05};
    ModuleRI::SternheimerFDHamiltonian hamiltonian(grid, potential);
    const auto states = ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 5, volume_element);

    const Vector& occupied = states.wavefunctions[0];
    constexpr double mixing_angle = 0.41;
    constexpr double mixing_phase = 0.37;
    const double cosine = std::cos(mixing_angle);
    const double sine = std::sin(mixing_angle);
    const Complex phase = std::exp(Complex(0.0, mixing_phase));

    Vector eta(grid.size(), Complex(0.0, 0.0));
    Vector out_direction(grid.size(), Complex(0.0, 0.0));
    for (std::size_t ir = 0; ir != eta.size(); ++ir)
    {
        eta[ir] = cosine * states.wavefunctions[2][ir] + sine * phase * states.wavefunctions[3][ir];
        out_direction[ir]
            = -sine * std::conj(phase) * states.wavefunctions[2][ir] + cosine * states.wavefunctions[3][ir];
    }

    auto dot = [volume_element](const Vector& lhs, const Vector& rhs) {
        return ModuleRI::sternheimer_fd_grid_dot(lhs, rhs, volume_element);
    };
    EXPECT_NEAR(std::abs(dot(occupied, eta)), 0.0, 1.0e-12);
    EXPECT_NEAR(std::abs(dot(occupied, out_direction)), 0.0, 1.0e-12);
    EXPECT_NEAR(std::abs(dot(eta, out_direction)), 0.0, 1.0e-12);

    Vector h_eta;
    Vector h_out;
    hamiltonian.apply(eta, h_eta);
    hamiltonian.apply(out_direction, h_out);
    const double fd_eta_diagonal = std::real(dot(eta, h_eta));
    const double delta_eigenvalue = fd_eta_diagonal + 0.37;
    EXPECT_GT(std::abs(delta_eigenvalue - fd_eta_diagonal), 0.3);

    Vector residual = h_eta;
    ModuleRI::SternheimerRPA::project_out_subspace({occupied, eta}, dot, residual);
    const Complex out_eta_coupling = dot(out_direction, residual);
    const Complex eta_out_coupling = std::conj(out_eta_coupling);

    const Complex rhs_eta(0.27, -0.19);
    const Complex rhs_out(-0.31, 0.23);
    Vector rhs(grid.size(), Complex(0.0, 0.0));
    for (std::size_t ir = 0; ir != rhs.size(); ++ir)
    {
        rhs[ir] = rhs_eta * eta[ir] + rhs_out * out_direction[ir];
    }

    constexpr double omega = 0.42;
    const double occupied_eigenvalue = states.eigenvalues[0];
    const Complex delta_block(delta_eigenvalue - occupied_eigenvalue, omega);
    const Complex out_block = dot(out_direction, h_out) - occupied_eigenvalue + Complex(0.0, omega);
    const Complex determinant
        = delta_block * out_block - eta_out_coupling * out_eta_coupling;
    const Complex expected_eta_coefficient
        = (rhs_eta * out_block - eta_out_coupling * rhs_out) / determinant;
    const Complex expected_out_coefficient
        = (delta_block * rhs_out - out_eta_coupling * rhs_eta) / determinant;

    ModuleRI::SternheimerRPA::SolverOptions solver_options;
    solver_options.max_iter = 80;
    solver_options.residual_tol = 1.0e-12;
    const ModuleRI::SternheimerDeltaVirtualState virtual_state{eta, residual, delta_eigenvalue};
    const Complex perturbation_matrix_element = -rhs_eta;
    const auto result = ModuleRI::solve_delta_sternheimer_linear_response(hamiltonian,
                                                                           {occupied},
                                                                           occupied_eigenvalue,
                                                                           rhs,
                                                                           {virtual_state},
                                                                           {perturbation_matrix_element},
                                                                           omega,
                                                                           volume_element,
                                                                           solver_options);

    ASSERT_TRUE(result.solver.converged);
    EXPECT_LT(result.residual_norm, 1.0e-9);
    Vector expected_out(grid.size(), Complex(0.0, 0.0));
    Vector expected_total(grid.size(), Complex(0.0, 0.0));
    for (std::size_t ir = 0; ir != expected_out.size(); ++ir)
    {
        expected_out[ir] = expected_out_coefficient * out_direction[ir];
        expected_total[ir] = expected_out[ir] + expected_eta_coefficient * eta[ir];
    }
    expect_vector_near(result.response.out_wavefunction, expected_out, 1.0e-9);
    expect_vector_near(result.response.reconstructed_wavefunction, expected_total, 1.0e-9);

    const Complex expected_sos = rhs_eta / delta_block;
    const Complex expected_pulay = -eta_out_coupling * expected_out_coefficient / delta_block;
    ASSERT_EQ(result.response.sos_coefficients.size(), 1);
    ASSERT_EQ(result.response.pulay_coefficients.size(), 1);
    EXPECT_NEAR(result.response.sos_coefficients[0].real(), expected_sos.real(), 1.0e-10);
    EXPECT_NEAR(result.response.sos_coefficients[0].imag(), expected_sos.imag(), 1.0e-10);
    EXPECT_NEAR(result.response.pulay_coefficients[0].real(), expected_pulay.real(), 1.0e-10);
    EXPECT_NEAR(result.response.pulay_coefficients[0].imag(), expected_pulay.imag(), 1.0e-10);
    EXPECT_NEAR(result.response.coefficients[0].real(), expected_eta_coefficient.real(), 1.0e-9);
    EXPECT_NEAR(result.response.coefficients[0].imag(), expected_eta_coefficient.imag(), 1.0e-9);
}

TEST(SternheimerDelta, ComplexQPerturbationMatrixElementsMatchGridIntegral)
{
    ModuleRI::SternheimerDeltaVirtualState eta0;
    eta0.orbital = {Complex(1.0, 1.0), Complex(0.5, -0.25)};
    eta0.residual = {Complex(0.0, 0.0), Complex(0.0, 0.0)};
    eta0.eigenvalue = 1.0;
    const std::vector<ModuleRI::SternheimerDeltaVirtualState> virtual_states = {eta0};
    const Vector perturbation = {Complex(2.0, -0.5), Complex(-1.0, 0.75)};
    const Vector occupied = {Complex(0.25, 0.5), Complex(-0.75, 1.0)};
    constexpr double volume_element = 0.2;

    const auto elements = ModuleRI::delta_sternheimer_perturbation_matrix_elements(
        virtual_states, perturbation, occupied, volume_element);
    const Complex expected = volume_element
        * (std::conj(eta0.orbital[0]) * perturbation[0] * occupied[0]
           + std::conj(eta0.orbital[1]) * perturbation[1] * occupied[1]);

    ASSERT_EQ(elements.size(), 1);
    EXPECT_NEAR(elements[0].real(), expected.real(), 1.0e-14);
    EXPECT_NEAR(elements[0].imag(), expected.imag(), 1.0e-14);
}
