#include "source_lcao/module_ri/sternheimer_periodic_solver.h"

#include <complex>
#include <gtest/gtest.h>
#include <vector>

TEST(SternheimerPeriodicSolver, StandardModeSolvesInFullGridComplement)
{
    using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
    using Complex = Hamiltonian::Complex;

    Hamiltonian::Grid grid{4, 1, 1, 0.5, 1.0, 1.0, true};
    constexpr double volume_element = 0.5;
    Hamiltonian hamiltonian(grid, std::vector<double>(grid.size(), 0.25));
    const auto states = ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 4, volume_element);
    const std::vector<Hamiltonian::Vector> occupied = {states.wavefunctions[0]};
    const Hamiltonian::Vector rhs = states.wavefunctions[3];
    constexpr double omega = 0.7;

    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 50;
    options.residual_tol = 1.0e-12;
    const auto response = ModuleRI::solve_sternheimer_periodic_linear_response(false,
                                                                                hamiltonian,
                                                                                occupied,
                                                                                states.eigenvalues[0],
                                                                                rhs,
                                                                                {},
                                                                                {},
                                                                                omega,
                                                                                volume_element,
                                                                                options);

    ASSERT_TRUE(response.solver.converged);
    EXPECT_FALSE(response.has_delta_components);
    EXPECT_LT(response.residual_norm, 1.0e-10);
    EXPECT_NEAR(response.full_grid_equation_residual_norm, response.residual_norm, 1.0e-14);
    const Complex factor = 1.0 / Complex(states.eigenvalues[3] - states.eigenvalues[0], omega);
    ASSERT_EQ(response.wavefunction.size(), rhs.size());
    for (std::size_t ir = 0; ir != rhs.size(); ++ir)
    {
        const Complex expected = factor * rhs[ir];
        EXPECT_NEAR(response.wavefunction[ir].real(), expected.real(), 1.0e-10);
        EXPECT_NEAR(response.wavefunction[ir].imag(), expected.imag(), 1.0e-10);
    }
}

TEST(SternheimerPeriodicSolver, DeltaModeReportsHybridAndFullGridResidualsSeparately)
{
    using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
    using Complex = Hamiltonian::Complex;

    Hamiltonian::Grid grid{4, 1, 1, 0.5, 1.0, 1.0, true};
    constexpr double volume_element = 0.5;
    Hamiltonian hamiltonian(grid, std::vector<double>(grid.size(), 0.25));
    const auto states = ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 4, volume_element);
    const std::vector<Hamiltonian::Vector> occupied = {states.wavefunctions[0]};
    const Hamiltonian::Vector& virtual_orbital = states.wavefunctions[2];
    Hamiltonian::Vector rhs = virtual_orbital;
    for (std::size_t ir = 0; ir != rhs.size(); ++ir)
    {
        rhs[ir] += Complex(0.25, 0.0) * states.wavefunctions[3][ir];
    }
    constexpr double omega = 0.7;
    constexpr double diagonal_shift = 0.3;

    ModuleRI::SternheimerDeltaVirtualState virtual_state;
    virtual_state.orbital = virtual_orbital;
    virtual_state.residual.assign(virtual_orbital.size(), Complex(0.0, 0.0));
    virtual_state.eigenvalue = states.eigenvalues[2] + diagonal_shift;

    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 50;
    options.residual_tol = 1.0e-12;
    const auto response = ModuleRI::solve_sternheimer_periodic_linear_response(true,
                                                                                hamiltonian,
                                                                                occupied,
                                                                                states.eigenvalues[0],
                                                                                rhs,
                                                                                {virtual_state},
                                                                                {Complex(-1.0, 0.0)},
                                                                                omega,
                                                                                volume_element,
                                                                                options);

    ASSERT_TRUE(response.solver.converged);
    EXPECT_TRUE(response.has_delta_components);
    EXPECT_LT(response.residual_norm, 1.0e-10);
    EXPECT_GT(response.full_grid_equation_residual_norm, 1.0e-2);
}
