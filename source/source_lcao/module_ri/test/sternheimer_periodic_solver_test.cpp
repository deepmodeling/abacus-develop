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
    Hamiltonian::Vector rhs = states.wavefunctions[3];
    for (std::size_t ir = 0; ir != rhs.size(); ++ir)
    {
        rhs[ir] += Complex(0.25, -0.125) * states.wavefunctions[0][ir];
    }
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
    ASSERT_EQ(response.projected_rhs.size(), states.wavefunctions[3].size());
    const Complex factor = 1.0 / Complex(states.eigenvalues[3] - states.eigenvalues[0], omega);
    ASSERT_EQ(response.wavefunction.size(), rhs.size());
    for (std::size_t ir = 0; ir != rhs.size(); ++ir)
    {
        EXPECT_NEAR(response.projected_rhs[ir].real(), states.wavefunctions[3][ir].real(), 1.0e-10);
        EXPECT_NEAR(response.projected_rhs[ir].imag(), states.wavefunctions[3][ir].imag(), 1.0e-10);
        const Complex expected = factor * states.wavefunctions[3][ir];
        EXPECT_NEAR(response.wavefunction[ir].real(), expected.real(), 1.0e-10);
        EXPECT_NEAR(response.wavefunction[ir].imag(), expected.imag(), 1.0e-10);
    }
}

TEST(SternheimerPeriodicSolver, StandardModeMatchesCompleteSameGridSOSWavefunctionAndResponse)
{
    using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
    using Complex = Hamiltonian::Complex;
    using Vector = Hamiltonian::Vector;

    Hamiltonian::Grid grid{6, 1, 1, 0.45, 1.0, 1.0, true};
    grid.kpoint = {0.17, 0.0, 0.0};
    constexpr double volume_element = 0.45;
    Hamiltonian hamiltonian(grid, {0.10, -0.07, 0.16, -0.11, 0.04, 0.13});
    const auto states = ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 6, volume_element);
    constexpr int occupied_count = 2;
    constexpr int occupied_index = 1;
    const std::vector<Vector> occupied(states.wavefunctions.begin(),
                                       states.wavefunctions.begin() + occupied_count);
    const std::vector<Vector> perturbations = {
        {Complex(0.70, 0.11), Complex(-0.30, 0.05), Complex(0.21, -0.08),
         Complex(0.48, 0.03), Complex(-0.37, 0.09), Complex(0.14, -0.06)},
        {Complex(-0.22, 0.07), Complex(0.51, -0.04), Complex(0.33, 0.12),
         Complex(-0.41, 0.02), Complex(0.19, -0.10), Complex(0.62, 0.08)}};
    constexpr double omega = 0.63;

    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 80;
    options.residual_tol = 1.0e-12;
    std::vector<Complex> standard_branch(perturbations.size() * perturbations.size(), Complex(0.0, 0.0));
    std::vector<Complex> sos_branch = standard_branch;
    for (int column = 0; column != static_cast<int>(perturbations.size()); ++column)
    {
        Vector rhs;
        ModuleRI::SternheimerRPA::build_rhs_from_hartree_perturbation(
            perturbations[static_cast<std::size_t>(column)], states.wavefunctions[occupied_index], rhs);
        const auto standard = ModuleRI::solve_sternheimer_periodic_linear_response(
            false,
            hamiltonian,
            occupied,
            states.eigenvalues[occupied_index],
            rhs,
            {},
            {},
            omega,
            volume_element,
            options);
        ASSERT_TRUE(standard.solver.converged);
        const Vector sos = ModuleRI::build_sternheimer_fd_complete_sos_response(states,
                                                                                occupied_count,
                                                                                occupied_index,
                                                                                rhs,
                                                                                omega,
                                                                                volume_element);
        ASSERT_EQ(standard.wavefunction.size(), sos.size());
        for (std::size_t ir = 0; ir != sos.size(); ++ir)
        {
            EXPECT_NEAR(standard.wavefunction[ir].real(), sos[ir].real(), 1.0e-10);
            EXPECT_NEAR(standard.wavefunction[ir].imag(), sos[ir].imag(), 1.0e-10);
        }
        ModuleRI::SternheimerRPA::accumulate_chi0_branch_column(perturbations,
                                                                states.wavefunctions[occupied_index],
                                                                standard.wavefunction,
                                                                volume_element,
                                                                1.0,
                                                                column,
                                                                standard_branch);
        ModuleRI::SternheimerRPA::accumulate_chi0_branch_column(perturbations,
                                                                states.wavefunctions[occupied_index],
                                                                sos,
                                                                volume_element,
                                                                1.0,
                                                                column,
                                                                sos_branch);
    }

    const auto standard_m
        = ModuleRI::SternheimerRPA::symmetrize_chi0_imaginary_frequency(standard_branch, perturbations.size());
    const auto sos_m
        = ModuleRI::SternheimerRPA::symmetrize_chi0_imaginary_frequency(sos_branch, perturbations.size());
    for (std::size_t index = 0; index != standard_branch.size(); ++index)
    {
        EXPECT_NEAR(standard_branch[index].real(), sos_branch[index].real(), 1.0e-10);
        EXPECT_NEAR(standard_branch[index].imag(), sos_branch[index].imag(), 1.0e-10);
        EXPECT_NEAR(standard_m[index].real(), sos_m[index].real(), 1.0e-10);
        EXPECT_NEAR(standard_m[index].imag(), sos_m[index].imag(), 1.0e-10);
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
