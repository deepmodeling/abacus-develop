#include "source_lcao/module_ri/sternheimer_fd_solver.h"

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

namespace
{

using Complex = std::complex<double>;
using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
using Vector = Hamiltonian::Vector;

} // namespace

TEST(SternheimerFDSolver, DenseZeroOrderStatesUsePhysicalGridNormalization)
{
    Hamiltonian::Grid grid{4, 1, 1, 0.5, 1.0, 1.0, true};
    const double volume_element = 0.5;
    const std::vector<double> potential(grid.size(), 0.25);
    Hamiltonian hamiltonian(grid, potential);

    const auto states = ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 3, volume_element);

    ASSERT_EQ(states.eigenvalues.size(), 3);
    ASSERT_EQ(states.wavefunctions.size(), 3);
    ASSERT_EQ(states.residual_norms.size(), 3);
    const std::vector<double> expected_eigenvalues = {0.25, 4.25, 4.25};
    for (std::size_t ib = 0; ib != states.eigenvalues.size(); ++ib)
    {
        EXPECT_NEAR(states.eigenvalues[ib], expected_eigenvalues[ib], 1.0e-12);
        EXPECT_NEAR(ModuleRI::sternheimer_fd_grid_norm(states.wavefunctions[ib], volume_element), 1.0, 1.0e-12);
        EXPECT_NEAR(states.residual_norms[ib], 0.0, 1.0e-12);
    }

    for (std::size_t ib = 0; ib != states.wavefunctions.size(); ++ib)
    {
        for (std::size_t jb = 0; jb != states.wavefunctions.size(); ++jb)
        {
            const Complex overlap
                = ModuleRI::sternheimer_fd_grid_dot(states.wavefunctions[ib], states.wavefunctions[jb], volume_element);
            const double expected = (ib == jb) ? 1.0 : 0.0;
            EXPECT_NEAR(overlap.real(), expected, 1.0e-12);
            EXPECT_NEAR(overlap.imag(), 0.0, 1.0e-12);
        }
    }
}

TEST(SternheimerFDSolver, DenseZeroOrderRejectsInvalidArguments)
{
    Hamiltonian::Grid grid{2, 1, 1, 1.0, 1.0, 1.0, true};
    Hamiltonian hamiltonian(grid, std::vector<double>(grid.size(), 0.0));

    EXPECT_THROW(ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 0, 1.0), std::invalid_argument);
    EXPECT_THROW(ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 3, 1.0), std::invalid_argument);
    EXPECT_THROW(ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 1, 0.0), std::invalid_argument);
}

TEST(SternheimerFDSolver, LinearResponseSolvesProjectedEigenmode)
{
    Hamiltonian::Grid grid{4, 1, 1, 0.5, 1.0, 1.0, true};
    const double volume_element = 0.5;
    Hamiltonian hamiltonian(grid, std::vector<double>(grid.size(), 0.25));
    const auto states = ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 4, volume_element);

    const std::vector<Vector> occupied = {states.wavefunctions[0]};
    const Vector rhs = states.wavefunctions[3];
    constexpr double omega = 0.7;
    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 50;
    options.residual_tol = 1.0e-12;

    const auto response = ModuleRI::solve_sternheimer_fd_linear_response(hamiltonian,
                                                                         occupied,
                                                                         states.eigenvalues[0],
                                                                         rhs,
                                                                         omega,
                                                                         volume_element,
                                                                         options);

    EXPECT_TRUE(response.solver.converged);
    EXPECT_NEAR(response.residual_norm, 0.0, 1.0e-10);
    EXPECT_NEAR(std::abs(ModuleRI::sternheimer_fd_grid_dot(occupied[0], response.delta_wavefunction, volume_element)),
                0.0,
                1.0e-12);

    const Complex factor = 1.0 / Complex(states.eigenvalues[3] - states.eigenvalues[0], omega);
    for (std::size_t ir = 0; ir != rhs.size(); ++ir)
    {
        const Complex expected = factor * rhs[ir];
        EXPECT_NEAR(response.delta_wavefunction[ir].real(), expected.real(), 1.0e-10);
        EXPECT_NEAR(response.delta_wavefunction[ir].imag(), expected.imag(), 1.0e-10);
    }
}
