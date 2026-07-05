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
