#include "source_lcao/module_ri/sternheimer_periodic_solver.h"
#include <limits>

#include <complex>
#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace
{
std::vector<ModuleRI::SternheimerDeltaVirtualState> excitation_test_states(
    const std::vector<double>& eigenvalues)
{
    std::vector<ModuleRI::SternheimerDeltaVirtualState> states(eigenvalues.size());
    for (std::size_t index = 0; index != eigenvalues.size(); ++index)
    {
        states[index].eigenvalue = eigenvalues[index];
        states[index].orbital = {{1.0, 0.25}};
        states[index].residual = {{0.0, -0.125}};
    }
    return states;
}
} // namespace

TEST(SternheimerExcitationGuard, AcceptsUnsortedPositiveSpectrumWithoutMutation)
{
    const std::vector<double> source = {-0.4, -0.8, -0.6};
    const auto target = excitation_test_states({0.3, -0.2, 0.1});
    const auto source_before = source;
    const auto target_before = target;
    EXPECT_DOUBLE_EQ(ModuleRI::require_positive_periodic_delta_excitations(target, source, "positive"),
                     -0.2 - (-0.4));
    EXPECT_EQ(source, source_before);
    ASSERT_EQ(target.size(), target_before.size());
    for (std::size_t index = 0; index != target.size(); ++index)
    {
        EXPECT_EQ(target[index].eigenvalue, target_before[index].eigenvalue);
        EXPECT_EQ(target[index].orbital, target_before[index].orbital);
        EXPECT_EQ(target[index].residual, target_before[index].residual);
    }
}

TEST(SternheimerExcitationGuard, AcceptsSmallPositiveGapWithoutFloor)
{
    EXPECT_DOUBLE_EQ(ModuleRI::require_positive_periodic_delta_excitations(
                         excitation_test_states({1.0e-14}), {0.0}, "small_positive"),
                     1.0e-14);
}

TEST(SternheimerExcitationGuard, RejectsZeroGap)
{
    EXPECT_THROW(ModuleRI::require_positive_periodic_delta_excitations(
                     excitation_test_states({0.3, -0.4}), {-0.8, -0.4}, "zero"),
                 ModuleRI::SternheimerExcitationError);
    EXPECT_THROW(ModuleRI::require_positive_periodic_delta_excitations(
                     excitation_test_states({-0.0}), {0.0}, "negative_zero"),
                 ModuleRI::SternheimerExcitationError);
}

TEST(SternheimerExcitationGuard, RejectsNegativeGap)
{
    EXPECT_THROW(ModuleRI::require_positive_periodic_delta_excitations(
                     excitation_test_states({0.3, -0.45}), {-0.8, -0.4}, "negative"),
                 ModuleRI::SternheimerExcitationError);
}

TEST(SternheimerExcitationGuard, UsesCompleteSourceMaximumNotFirstBand)
{
    // The first source band alone would give a positive gap of 0.2 Ry.
    EXPECT_THROW(ModuleRI::require_positive_periodic_delta_excitations(
                     excitation_test_states({-0.6, 0.3}), {-0.8, -0.4}, "complete_source"),
                 ModuleRI::SternheimerExcitationError);
}

TEST(SternheimerExcitationGuard, UsesSuppliedSourceReferenceNotTargetOccupiedReference)
{
    const auto target = excitation_test_states({-0.5, 0.3});
    const std::vector<double> target_occupied = {-0.7};
    EXPECT_DOUBLE_EQ(ModuleRI::require_positive_periodic_delta_excitations(
                         target, target_occupied, "target_reference_control"),
                     -0.5 - (-0.7));
    EXPECT_THROW(ModuleRI::require_positive_periodic_delta_excitations(
                     target, {-0.4}, "actual_source_reference"),
                 ModuleRI::SternheimerExcitationError);
}

TEST(SternheimerExcitationGuard, RejectsEmptySpectra)
{
    EXPECT_THROW(ModuleRI::require_positive_periodic_delta_excitations({}, {-0.4}, "empty_target"),
                 ModuleRI::SternheimerExcitationError);
    EXPECT_THROW(ModuleRI::require_positive_periodic_delta_excitations(
                     excitation_test_states({0.2}), {}, "empty_source"),
                 ModuleRI::SternheimerExcitationError);
    EXPECT_THROW(ModuleRI::require_positive_periodic_delta_excitations({}, {}, "both_empty"),
                 ModuleRI::SternheimerExcitationError);
}

TEST(SternheimerExcitationGuard, RejectsEveryNonfiniteSourcePosition)
{
    const double invalid_values[] = {std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::infinity(),
                                     -std::numeric_limits<double>::infinity()};
    for (const double invalid: invalid_values)
    {
        for (std::size_t index = 0; index != 3; ++index)
        {
            std::vector<double> source = {-0.8, -0.4, -0.6};
            source[index] = invalid;
            try
            {
                ModuleRI::require_positive_periodic_delta_excitations(
                    excitation_test_states({0.3, -0.2}), source, "nonfinite_source_fixture");
                FAIL() << "Accepted nonfinite source energy at " << index;
            }
            catch (const ModuleRI::SternheimerExcitationError& error)
            {
                const std::string message = error.what();
                EXPECT_NE(message.find("reason=nonfinite_source"), std::string::npos);
                EXPECT_NE(message.find("source_index=" + std::to_string(index + 1)), std::string::npos);
                EXPECT_NE(message.find("context={nonfinite_source_fixture}"), std::string::npos);
            }
        }
    }
}

TEST(SternheimerExcitationGuard, RejectsEveryNonfiniteTargetPosition)
{
    const double invalid_values[] = {std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::infinity(),
                                     -std::numeric_limits<double>::infinity()};
    for (const double invalid: invalid_values)
    {
        for (std::size_t index = 0; index != 3; ++index)
        {
            std::vector<double> target = {0.3, -0.2, 0.1};
            target[index] = invalid;
            try
            {
                ModuleRI::require_positive_periodic_delta_excitations(
                    excitation_test_states(target), {-0.8, -0.4}, "nonfinite_target_fixture");
                FAIL() << "Accepted nonfinite target energy at " << index;
            }
            catch (const ModuleRI::SternheimerExcitationError& error)
            {
                const std::string message = error.what();
                EXPECT_NE(message.find("reason=nonfinite_target"), std::string::npos);
                EXPECT_NE(message.find("target_index=" + std::to_string(index + 1)), std::string::npos);
                EXPECT_NE(message.find("context={nonfinite_target_fixture}"), std::string::npos);
            }
        }
    }
}

TEST(SternheimerExcitationGuard, RejectsOverflowFromFiniteExtrema)
{
    for (const double source: {std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()})
    {
        try
        {
            ModuleRI::require_positive_periodic_delta_excitations(
                excitation_test_states({-source}), {source}, "overflow_fixture");
            FAIL() << "Accepted nonfinite gap from finite extrema";
        }
        catch (const ModuleRI::SternheimerExcitationError& error)
        {
            EXPECT_NE(std::string(error.what()).find("reason=nonfinite_gap"), std::string::npos);
        }
    }
}

TEST(SternheimerExcitationGuard, ReportsFullPrecisionBoundsIndicesAndContext)
{
    const double source_max = -0.40123456789012345;
    const double target_min = -0.44634567890123456;
    const std::string context = "grid=24,24,48 pbe_grid=72,72,288 q=2 source_k=1 target_k=2 rank=3";
    try
    {
        ModuleRI::require_positive_periodic_delta_excitations(
            excitation_test_states({0.3, target_min, -0.42}), {-0.7, source_max, -0.5}, context);
        FAIL() << "Accepted inverted spectrum";
    }
    catch (const ModuleRI::SternheimerExcitationError& error)
    {
        const std::runtime_error& base_error = error;
        const std::string message = base_error.what();
        EXPECT_NE(message.find("context={" + context + "}"), std::string::npos);
        EXPECT_NE(message.find("units=Ry index_base=1"), std::string::npos);
        EXPECT_NE(message.find("source_count=3 target_count=3"), std::string::npos);
        EXPECT_NE(message.find("source_max_index=2"), std::string::npos);
        EXPECT_NE(message.find("target_min_index=2"), std::string::npos);
        EXPECT_NE(message.find("reason=nonpositive_gap"), std::string::npos);
        const std::string fields[] = {"source_max_Ry=", "target_min_Ry=", "min_excitation_Ry="};
        const double expected[] = {source_max, target_min, target_min - source_max};
        for (std::size_t index = 0; index != 3; ++index)
        {
            const auto position = message.find(fields[index]);
            ASSERT_NE(position, std::string::npos);
            EXPECT_EQ(std::stod(message.substr(position + fields[index].size())), expected[index]);
        }
    }
}

TEST(SternheimerExcitationGuard, ReportsFirstIndexForTiedExtrema)
{
    try
    {
        ModuleRI::require_positive_periodic_delta_excitations(
            excitation_test_states({-0.5, -0.5}), {-0.4, -0.4}, "ties");
        FAIL() << "Accepted inverted spectrum";
    }
    catch (const ModuleRI::SternheimerExcitationError& error)
    {
        const std::string message = error.what();
        EXPECT_NE(message.find("source_max_index=1"), std::string::npos);
        EXPECT_NE(message.find("target_min_index=1"), std::string::npos);
    }
}

TEST(SternheimerPeriodicSolver, SharedProjectorReservationIncludesMetadataAndRejectsOverflow)
{
    using Projectors = ModuleRI::SternheimerPeriodicResponseProjectors;
    using Vector = ModuleRI::SternheimerFDHamiltonian::Vector;
    EXPECT_EQ(Projectors::estimated_storage_bytes(10, 2, 3),
              7u * (11u * sizeof(ModuleRI::SternheimerFDHamiltonian::Complex) + sizeof(const Vector*)));
    EXPECT_EQ(Projectors::estimated_storage_bytes(10, 2, 3, false),
              7u * (sizeof(ModuleRI::SternheimerFDHamiltonian::Complex) + sizeof(const Vector*)));
    EXPECT_THROW(Projectors::estimated_storage_bytes(std::numeric_limits<std::size_t>::max(), 2, 3),
                 std::overflow_error);
    EXPECT_THROW(Projectors::estimated_storage_bytes(10, std::numeric_limits<std::size_t>::max(), 3),
                 std::overflow_error);
}

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

TEST(SternheimerPeriodicSolver, SpectralPreconditionerPreservesResponseAndReducesIterations)
{
    using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
    using Complex = Hamiltonian::Complex;
    using Vector = Hamiltonian::Vector;

    Hamiltonian::Grid grid{48, 1, 1, 0.35, 1.0, 1.0, true};
    const double pi = std::acos(-1.0);
    std::vector<double> potential(static_cast<std::size_t>(grid.size()));
    for (int ix = 0; ix != grid.nx; ++ix)
    {
        potential[static_cast<std::size_t>(ix)]
            = 0.18 * std::cos(2.0 * pi * ix / grid.nx)
              + 0.07 * std::sin(6.0 * pi * ix / grid.nx);
    }
    constexpr double volume_element = 0.35;
    Hamiltonian hamiltonian(grid, potential, 1.0, nullptr, 8);
    const auto states = ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 48, volume_element);
    const std::vector<Vector> occupied = {states.wavefunctions[0]};
    Vector rhs(static_cast<std::size_t>(grid.size()));
    for (int ix = 0; ix != grid.nx; ++ix)
    {
        rhs[static_cast<std::size_t>(ix)]
            = Complex(std::cos(10.0 * pi * ix / grid.nx),
                      0.3 * std::sin(14.0 * pi * ix / grid.nx));
    }

    ModuleRI::SternheimerRPA::SolverOptions identity_options;
    identity_options.max_iter = 150;
    identity_options.residual_tol = 1.0e-11;
    identity_options.use_fd_spectral_preconditioner = false;
    const auto identity = ModuleRI::solve_sternheimer_periodic_linear_response(false,
                                                                                hamiltonian,
                                                                                occupied,
                                                                                states.eigenvalues[0],
                                                                                rhs,
                                                                                {},
                                                                                {},
                                                                                0.08,
                                                                                volume_element,
                                                                                identity_options);

    ModuleRI::SternheimerRPA::SolverOptions preconditioned_options = identity_options;
    preconditioned_options.use_fd_spectral_preconditioner = true;
    preconditioned_options.fd_spectral_preconditioner_regularization = 0.2;
    const auto preconditioned = ModuleRI::solve_sternheimer_periodic_linear_response(false,
                                                                                      hamiltonian,
                                                                                      occupied,
                                                                                      states.eigenvalues[0],
                                                                                      rhs,
                                                                                      {},
                                                                                      {},
                                                                                      0.08,
                                                                                      volume_element,
                                                                                      preconditioned_options);

    ASSERT_TRUE(identity.solver.converged);
    ASSERT_TRUE(preconditioned.solver.converged);
    EXPECT_LT(preconditioned.solver.iterations, identity.solver.iterations);
    EXPECT_LT(preconditioned.residual_norm, 1.0e-9);
    ASSERT_EQ(preconditioned.wavefunction.size(), identity.wavefunction.size());
    for (std::size_t ir = 0; ir != identity.wavefunction.size(); ++ir)
    {
        EXPECT_NEAR(preconditioned.wavefunction[ir].real(), identity.wavefunction[ir].real(), 1.0e-8);
        EXPECT_NEAR(preconditioned.wavefunction[ir].imag(), identity.wavefunction[ir].imag(), 1.0e-8);
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

TEST(SternheimerPeriodicSolver, DeltaBatchMatchesIndependentChannels)
{
    using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
    using Complex = Hamiltonian::Complex;
    using Vector = Hamiltonian::Vector;

    Hamiltonian::Grid grid{6, 1, 1, 0.45, 1.0, 1.0, true};
    grid.kpoint = {0.13, 0.0, 0.0};
    constexpr double volume_element = 0.45;
    Hamiltonian hamiltonian(grid, {0.11, -0.05, 0.16, -0.09, 0.03, 0.14}, 1.0, nullptr, 8);
    const auto states = ModuleRI::solve_sternheimer_fd_zero_order_dense(hamiltonian, 6, volume_element);
    const std::vector<Vector> occupied = {states.wavefunctions[0]};
    std::vector<ModuleRI::SternheimerDeltaVirtualState> virtual_states;
    for (int index = 1; index != 3; ++index)
    {
        ModuleRI::SternheimerDeltaVirtualState state;
        state.orbital = states.wavefunctions[static_cast<std::size_t>(index)];
        state.residual = states.wavefunctions[4];
        for (Complex& value: state.residual)
        {
            value *= Complex(0.01 * index, -0.005 * index);
        }
        state.eigenvalue = states.eigenvalues[static_cast<std::size_t>(index)] + 0.02 * index;
        virtual_states.push_back(std::move(state));
    }
    const Hamiltonian::Matrix rhs = {
        states.wavefunctions[3],
        states.wavefunctions[5]};
    Hamiltonian::Matrix coupled_rhs = rhs;
    for (std::size_t column = 0; column != coupled_rhs.size(); ++column)
    {
        for (std::size_t ir = 0; ir != coupled_rhs[column].size(); ++ir)
        {
            coupled_rhs[column][ir] += Complex(0.3, 0.2 * (column + 1)) * states.wavefunctions[1][ir]
                                      + Complex(-0.2, 0.15) * states.wavefunctions[2][ir]
                                      + Complex(0.4, -0.1) * states.wavefunctions[4][ir];
        }
    }
    std::vector<std::vector<Complex>> perturbation_matrix_elements(rhs.size());
    for (std::size_t column = 0; column != rhs.size(); ++column)
    {
        for (const auto& state: virtual_states)
        {
            perturbation_matrix_elements[column].push_back(
                -ModuleRI::sternheimer_fd_grid_dot(state.orbital, coupled_rhs[column], volume_element));
        }
    }
    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 100;
    options.residual_tol = 1.0e-12;
    options.use_fd_spectral_preconditioner = true;
    constexpr double omega = 0.57;

    std::vector<ModuleRI::SternheimerPeriodicLinearResponse> scalar;
    for (std::size_t column = 0; column != rhs.size(); ++column)
    {
        scalar.push_back(ModuleRI::solve_sternheimer_periodic_linear_response(true,
                                                                               hamiltonian,
                                                                               occupied,
                                                                               states.eigenvalues[0],
                                                                               coupled_rhs[column],
                                                                               virtual_states,
                                                                               perturbation_matrix_elements[column],
                                                                               omega,
                                                                               volume_element,
                                                                               options));
    }
    const auto batch = ModuleRI::solve_sternheimer_periodic_linear_response_batch(
        true,
        hamiltonian,
        occupied,
        states.eigenvalues[0],
        coupled_rhs,
        virtual_states,
        perturbation_matrix_elements,
        omega,
        volume_element,
        options);

    const ModuleRI::SternheimerPeriodicResponseProjectors projectors(
        occupied, virtual_states, volume_element);
    const ModuleRI::SternheimerPeriodicResponseProjectors scalar_projectors(
        occupied, virtual_states, volume_element, false);
    EXPECT_LT(scalar_projectors.fixed.storage_bytes(), projectors.fixed.storage_bytes());
    const auto shared_batch = ModuleRI::solve_sternheimer_periodic_linear_response_batch(
        true, hamiltonian, occupied, states.eigenvalues[0], coupled_rhs, virtual_states,
        perturbation_matrix_elements, omega, volume_element, options, &projectors);
    const auto owning_subspace = ModuleRI::build_delta_sternheimer_fixed_subspace(occupied, virtual_states);
    const auto owning_batch = ModuleRI::solve_delta_sternheimer_linear_response_batch(
        hamiltonian, owning_subspace, states.eigenvalues[0], coupled_rhs, virtual_states,
        perturbation_matrix_elements, omega, volume_element, options);
    ASSERT_EQ(shared_batch.size(), batch.size());

    ASSERT_EQ(batch.size(), scalar.size());
    for (std::size_t column = 0; column != scalar.size(); ++column)
    {
        const auto shared_scalar = ModuleRI::solve_sternheimer_periodic_linear_response(
            true, hamiltonian, occupied, states.eigenvalues[0], coupled_rhs[column], virtual_states,
            perturbation_matrix_elements[column], omega, volume_element, options, &scalar_projectors);
        EXPECT_EQ(shared_batch[column].solver.converged, owning_batch[column].solver.converged);
        EXPECT_EQ(shared_batch[column].solver.iterations, owning_batch[column].solver.iterations);
        EXPECT_NEAR(shared_batch[column].residual_norm, owning_batch[column].residual_norm, 1.0e-12);
        EXPECT_NEAR(shared_scalar.residual_norm, scalar[column].residual_norm, 1.0e-12);
        EXPECT_GT(std::abs(shared_batch[column].delta_components.sos_coefficients[0]), 1.0e-4);
        EXPECT_GT(std::abs(shared_batch[column].delta_components.pulay_coefficients[0]), 1.0e-6);
        for (std::size_t state = 0; state != virtual_states.size(); ++state)
        {
            EXPECT_NEAR(std::abs(shared_batch[column].delta_components.coefficients[state]
                                - owning_batch[column].response.coefficients[state]), 0.0, 1.0e-12);
        }
        EXPECT_EQ(batch[column].solver.converged, scalar[column].solver.converged);
        EXPECT_EQ(batch[column].solver.iterations, scalar[column].solver.iterations);
        EXPECT_NEAR(batch[column].solver.relative_residual, scalar[column].solver.relative_residual, 1.0e-12);
        EXPECT_NEAR(batch[column].residual_norm, scalar[column].residual_norm, 1.0e-11);
        EXPECT_NEAR(batch[column].full_grid_equation_residual_norm,
                    scalar[column].full_grid_equation_residual_norm,
                    1.0e-11);
        ASSERT_EQ(batch[column].wavefunction.size(), scalar[column].wavefunction.size());
        for (std::size_t ir = 0; ir != scalar[column].wavefunction.size(); ++ir)
        {
            EXPECT_NEAR(std::abs(shared_batch[column].wavefunction[ir]
                                 - owning_batch[column].response.reconstructed_wavefunction[ir]), 0.0, 1.0e-12);
            EXPECT_NEAR(std::abs(shared_scalar.wavefunction[ir] - scalar[column].wavefunction[ir]), 0.0, 1.0e-12);
            EXPECT_NEAR(batch[column].projected_rhs[ir].real(), scalar[column].projected_rhs[ir].real(), 1.0e-12);
            EXPECT_NEAR(batch[column].projected_rhs[ir].imag(), scalar[column].projected_rhs[ir].imag(), 1.0e-12);
            EXPECT_NEAR(batch[column].wavefunction[ir].real(), scalar[column].wavefunction[ir].real(), 1.0e-11);
            EXPECT_NEAR(batch[column].wavefunction[ir].imag(), scalar[column].wavefunction[ir].imag(), 1.0e-11);
        }
    }
}
