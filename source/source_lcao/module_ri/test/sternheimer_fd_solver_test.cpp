#include "source_lcao/module_ri/sternheimer_fd_solver.h"

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

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

TEST(SternheimerFDSolver, GridDotUsesRequestedOpenMPThreadsAndMatchesSerial)
{
#ifndef _OPENMP
    GTEST_SKIP() << "OpenMP is not enabled in this build.";
#else
    constexpr std::size_t size = 16384;
    Vector lhs(size);
    Vector rhs(size);
    for (std::size_t ir = 0; ir != size; ++ir)
    {
        lhs[ir] = Complex(0.01 * static_cast<double>(ir % 23), -0.02 * static_cast<double>(ir % 19));
        rhs[ir] = Complex(-0.03 * static_cast<double>(ir % 13), 0.04 * static_cast<double>(ir % 17));
    }

    const int previous_threads = omp_get_max_threads();
    omp_set_dynamic(0);
    omp_set_num_threads(1);
    int serial_threads = 0;
    const Complex serial = ModuleRI::sternheimer_fd_grid_dot(lhs, rhs, 0.125, &serial_threads);

    omp_set_num_threads(4);
    int parallel_threads = 0;
    const Complex parallel = ModuleRI::sternheimer_fd_grid_dot(lhs, rhs, 0.125, &parallel_threads);
    omp_set_num_threads(previous_threads);

    EXPECT_EQ(serial_threads, 1);
    EXPECT_EQ(parallel_threads, 4);
    EXPECT_NEAR(parallel.real(), serial.real(), 1.0e-11);
    EXPECT_NEAR(parallel.imag(), serial.imag(), 1.0e-11);
#endif
}

TEST(SternheimerFDSolver, LanczosZeroOrderMatchesDenseLowStates)
{
    Hamiltonian::Grid grid{6, 1, 1, 0.75, 1.0, 1.0, true};
    const double volume_element = 0.75;
    const std::vector<double> potential = {0.30, -0.15, 0.05, 0.20, -0.05, 0.10};
    Hamiltonian hamiltonian(grid, potential);

    const auto dense_states = ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 3, volume_element);

    ModuleRI::SternheimerFDLanczosOptions options;
    options.max_subspace_size = 12;
    options.residual_tolerance = 1.0e-11;
    options.initial_seed = 17;
    const auto lanczos_states
        = ModuleRI::solve_sternheimer_fd_zero_order_lanczos(hamiltonian, 3, volume_element, options);

    ASSERT_EQ(lanczos_states.eigenvalues.size(), dense_states.eigenvalues.size());
    ASSERT_EQ(lanczos_states.wavefunctions.size(), dense_states.wavefunctions.size());
    for (std::size_t ib = 0; ib != dense_states.eigenvalues.size(); ++ib)
    {
        EXPECT_NEAR(lanczos_states.eigenvalues[ib], dense_states.eigenvalues[ib], 1.0e-9);
        EXPECT_NEAR(ModuleRI::sternheimer_fd_grid_norm(lanczos_states.wavefunctions[ib], volume_element),
                    1.0,
                    1.0e-10);
        EXPECT_LT(lanczos_states.residual_norms[ib], 1.0e-9);
    }
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
