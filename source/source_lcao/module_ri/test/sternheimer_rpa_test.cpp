#include "source_lcao/module_ri/sternheimer_rpa.h"

#include <complex>
#include <gtest/gtest.h>
#include <vector>

namespace
{

using Complex = std::complex<double>;
using Vector = ModuleRI::SternheimerRPA::Vector;

Complex dot(const Vector& lhs, const Vector& rhs)
{
    return ModuleRI::SternheimerRPA::local_grid_dot(lhs, rhs, 1.0);
}

} // namespace

TEST(SternheimerRPA, BuildRhsFromHartreePerturbation)
{
    const std::vector<double> v_r = {2.0, -1.0};
    const Vector psi_r = {Complex(1.0, 1.0), Complex(2.0, 0.0)};

    Vector rhs;
    ModuleRI::SternheimerRPA::build_rhs_from_hartree_perturbation(v_r, psi_r, rhs);

    ASSERT_EQ(rhs.size(), 2);
    EXPECT_DOUBLE_EQ(rhs[0].real(), -2.0);
    EXPECT_DOUBLE_EQ(rhs[0].imag(), -2.0);
    EXPECT_DOUBLE_EQ(rhs[1].real(), 2.0);
    EXPECT_DOUBLE_EQ(rhs[1].imag(), 0.0);
}

TEST(SternheimerRPA, AccumulatePolarizabilityGridElement)
{
    const std::vector<double> v_r = {2.0, -1.0};
    const Vector psi_r = {Complex(1.0, 1.0), Complex(2.0, 0.0)};
    const Vector delta_psi_r = {Complex(0.5, 0.0), Complex(0.0, 1.0)};
    constexpr double grid_weight = 0.25;

    const Complex value
        = ModuleRI::SternheimerRPA::accumulate_polarizability_grid_element(v_r, psi_r, delta_psi_r, grid_weight);
    const Complex expected
        = grid_weight * (std::conj(psi_r[0]) * v_r[0] * delta_psi_r[0] + std::conj(psi_r[1]) * v_r[1] * delta_psi_r[1]);

    EXPECT_NEAR(value.real(), expected.real(), 1.0e-14);
    EXPECT_NEAR(value.imag(), expected.imag(), 1.0e-14);
}

TEST(SternheimerRPA, ProjectOutSubspace)
{
    const std::vector<Vector> occupied
        = {{Complex(1.0, 0.0), Complex(0.0, 0.0)}, {Complex(0.0, 0.0), Complex(0.0, 1.0)}};
    Vector vec = {Complex(2.0, 1.0), Complex(3.0, -4.0)};

    ModuleRI::SternheimerRPA::project_out_subspace(occupied, dot, vec);

    EXPECT_NEAR(std::abs(dot(occupied[0], vec)), 0.0, 1.0e-14);
    EXPECT_NEAR(std::abs(dot(occupied[1], vec)), 0.0, 1.0e-14);
}

TEST(SternheimerRPA, ApplyKineticPreconditioner)
{
    const std::vector<double> kinetic = {1.0, 3.0};
    const Vector input = {Complex(2.0, 0.0), Complex(0.0, 4.0)};
    Vector output;

    ModuleRI::SternheimerRPA::apply_kinetic_preconditioner(kinetic, 0.5, 2.0, 0.1, input, output);

    ASSERT_EQ(output.size(), input.size());
    EXPECT_NEAR(output[0].real(), (input[0] / Complex(0.6, 2.0)).real(), 1.0e-14);
    EXPECT_NEAR(output[0].imag(), (input[0] / Complex(0.6, 2.0)).imag(), 1.0e-14);
    EXPECT_NEAR(output[1].real(), (input[1] / Complex(2.6, 2.0)).real(), 1.0e-14);
    EXPECT_NEAR(output[1].imag(), (input[1] / Complex(2.6, 2.0)).imag(), 1.0e-14);
}

TEST(SternheimerRPA, SolveBiCGStabDiagonalComplexSystem)
{
    const Vector diagonal = {Complex(2.0, 1.0), Complex(3.0, -0.5), Complex(4.0, 2.0)};
    const Vector exact = {Complex(1.0, -0.5), Complex(-2.0, 1.0), Complex(0.25, 0.75)};
    Vector rhs(exact.size());
    for (std::size_t i = 0; i != exact.size(); ++i)
    {
        rhs[i] = diagonal[i] * exact[i];
    }

    ModuleRI::SternheimerRPA::LinearProblem problem;
    problem.apply = [&diagonal](const Vector& input, Vector& output) {
        output.resize(input.size());
        for (std::size_t i = 0; i != input.size(); ++i)
        {
            output[i] = diagonal[i] * input[i];
        }
    };
    problem.dot = dot;

    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 50;
    options.residual_tol = 1.0e-12;

    Vector solution(exact.size(), Complex(0.0, 0.0));
    const auto result = ModuleRI::SternheimerRPA::solve_bicgstab(problem, rhs, solution, options);

    EXPECT_TRUE(result.converged);
    EXPECT_LE(result.relative_residual, options.residual_tol);
    for (std::size_t i = 0; i != exact.size(); ++i)
    {
        EXPECT_NEAR(solution[i].real(), exact[i].real(), 1.0e-10);
        EXPECT_NEAR(solution[i].imag(), exact[i].imag(), 1.0e-10);
    }
}

TEST(SternheimerRPA, SolveGMRESDiagonalComplexSystem)
{
    const Vector diagonal = {Complex(2.0, 1.0), Complex(3.0, -0.5), Complex(4.0, 2.0)};
    const Vector exact = {Complex(1.0, -0.5), Complex(-2.0, 1.0), Complex(0.25, 0.75)};
    Vector rhs(exact.size());
    for (std::size_t i = 0; i != exact.size(); ++i)
    {
        rhs[i] = diagonal[i] * exact[i];
    }

    ModuleRI::SternheimerRPA::LinearProblem problem;
    problem.apply = [&diagonal](const Vector& input, Vector& output) {
        output.resize(input.size());
        for (std::size_t i = 0; i != input.size(); ++i)
        {
            output[i] = diagonal[i] * input[i];
        }
    };
    problem.dot = dot;

    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 10;
    options.residual_tol = 1.0e-12;

    Vector solution(exact.size(), Complex(0.0, 0.0));
    const auto result = ModuleRI::SternheimerRPA::solve_gmres(problem, rhs, solution, options, 3);

    EXPECT_TRUE(result.converged);
    EXPECT_LE(result.relative_residual, options.residual_tol);
    for (std::size_t i = 0; i != exact.size(); ++i)
    {
        EXPECT_NEAR(solution[i].real(), exact[i].real(), 1.0e-10);
        EXPECT_NEAR(solution[i].imag(), exact[i].imag(), 1.0e-10);
    }
}
