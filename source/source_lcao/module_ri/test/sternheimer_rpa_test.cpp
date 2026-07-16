#include "source_lcao/module_ri/sternheimer_rpa.h"

#include <complex>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <string>
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

TEST(SternheimerRPA, ComplexQPotentialBuildsRhsAndUsesMinusQProbe)
{
    const Vector v_q = {Complex(2.0, 1.0), Complex(-1.0, 0.5)};
    const Vector psi_k = {Complex(1.0, -1.0), Complex(0.5, 2.0)};
    const Vector delta_psi_kq = {Complex(0.25, 0.5), Complex(-1.0, 0.75)};
    constexpr double grid_weight = 0.2;

    Vector rhs;
    ModuleRI::SternheimerRPA::build_rhs_from_hartree_perturbation(v_q, psi_k, rhs);
    ASSERT_EQ(rhs.size(), psi_k.size());
    for (std::size_t ir = 0; ir != rhs.size(); ++ir)
    {
        EXPECT_NEAR(rhs[ir].real(), (-v_q[ir] * psi_k[ir]).real(), 1.0e-14);
        EXPECT_NEAR(rhs[ir].imag(), (-v_q[ir] * psi_k[ir]).imag(), 1.0e-14);
    }

    const Complex value = ModuleRI::SternheimerRPA::accumulate_polarizability_grid_element(
        v_q, psi_k, delta_psi_kq, grid_weight);
    const Complex expected
        = grid_weight
          * (std::conj(psi_k[0]) * std::conj(v_q[0]) * delta_psi_kq[0]
             + std::conj(psi_k[1]) * std::conj(v_q[1]) * delta_psi_kq[1]);
    EXPECT_NEAR(value.real(), expected.real(), 1.0e-14);
    EXPECT_NEAR(value.imag(), expected.imag(), 1.0e-14);
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

TEST(SternheimerRPA, AccumulateChi0BranchColumnIncludesOccupation)
{
    const std::vector<std::vector<double>> potentials = {{2.0, -1.0}, {0.5, 3.0}};
    const Vector psi_r = {Complex(1.0, 1.0), Complex(2.0, 0.0)};
    const Vector delta_psi_r = {Complex(0.5, 0.0), Complex(0.0, 1.0)};
    constexpr double grid_weight = 0.25;
    constexpr double occupation = 2.0;

    std::vector<Complex> branch(4, Complex(0.0, 0.0));
    ModuleRI::SternheimerRPA::accumulate_chi0_branch_column(potentials,
                                                            psi_r,
                                                            delta_psi_r,
                                                            grid_weight,
                                                            occupation,
                                                            1,
                                                            branch);

    const Complex expected0 = occupation
        * ModuleRI::SternheimerRPA::accumulate_polarizability_grid_element(potentials[0],
                                                                          psi_r,
                                                                          delta_psi_r,
                                                                          grid_weight);
    const Complex expected1 = occupation
        * ModuleRI::SternheimerRPA::accumulate_polarizability_grid_element(potentials[1],
                                                                          psi_r,
                                                                          delta_psi_r,
                                                                          grid_weight);
    EXPECT_EQ(branch[0], Complex(0.0, 0.0));
    EXPECT_EQ(branch[2], Complex(0.0, 0.0));
    EXPECT_NEAR(branch[1].real(), expected0.real(), 1.0e-14);
    EXPECT_NEAR(branch[1].imag(), expected0.imag(), 1.0e-14);
    EXPECT_NEAR(branch[3].real(), expected1.real(), 1.0e-14);
    EXPECT_NEAR(branch[3].imag(), expected1.imag(), 1.0e-14);
}

TEST(SternheimerRPA, AccumulateComplexQChi0BranchColumnUsesMinusQProbe)
{
    const std::vector<Vector> potentials = {
        {Complex(2.0, 1.0), Complex(-1.0, 0.5)},
        {Complex(0.5, -0.25), Complex(3.0, 2.0)}};
    const Vector psi_k = {Complex(1.0, 1.0), Complex(2.0, -0.5)};
    const Vector delta_psi_kq = {Complex(0.5, -0.25), Complex(-0.75, 1.0)};
    constexpr double grid_weight = 0.25;
    constexpr double occupation = 1.5;

    std::vector<Complex> branch(4, Complex(0.0, 0.0));
    ModuleRI::SternheimerRPA::accumulate_chi0_branch_column(potentials,
                                                            psi_k,
                                                            delta_psi_kq,
                                                            grid_weight,
                                                            occupation,
                                                            0,
                                                            branch);

    for (std::size_t row = 0; row != potentials.size(); ++row)
    {
        const Complex expected = occupation
            * ModuleRI::SternheimerRPA::accumulate_polarizability_grid_element(
                potentials[row], psi_k, delta_psi_kq, grid_weight);
        EXPECT_NEAR(branch[2 * row].real(), expected.real(), 1.0e-14);
        EXPECT_NEAR(branch[2 * row].imag(), expected.imag(), 1.0e-14);
        EXPECT_EQ(branch[2 * row + 1], Complex(0.0, 0.0));
    }
}

TEST(SternheimerRPA, SymmetrizeImaginaryFrequencyChi0AddsAdjointBranch)
{
    const std::vector<Complex> branch = {Complex(1.0, 2.0),
                                         Complex(3.0, -4.0),
                                         Complex(5.0, 6.0),
                                         Complex(-7.0, 8.0)};

    const std::vector<Complex> chi0 = ModuleRI::SternheimerRPA::symmetrize_chi0_imaginary_frequency(branch, 2);

    ASSERT_EQ(chi0.size(), 4);
    EXPECT_EQ(chi0[0], branch[0] + std::conj(branch[0]));
    EXPECT_EQ(chi0[1], branch[1] + std::conj(branch[2]));
    EXPECT_EQ(chi0[2], branch[2] + std::conj(branch[1]));
    EXPECT_EQ(chi0[3], branch[3] + std::conj(branch[3]));
    EXPECT_NE(chi0[1], Complex(2.0 * branch[1].real(), 0.0));
}

TEST(SternheimerRPA, WriteChi0V1FileUsesAtomPairBlocks)
{
    ModuleRI::SternheimerRPA::Chi0V1Metadata metadata;
    metadata.iq = 1;
    metadata.ifrequency = 2;
    metadata.omega = 0.5;
    metadata.weight = 0.125;
    metadata.atom_naux = {1, 1};

    const std::vector<ModuleRI::SternheimerRPA::AuxiliaryChannel> channels
        = {{0, 0, 0}, {1, 1, 0}};
    const std::vector<Complex> chi0 = {Complex(1.0, 0.0),
                                       Complex(2.0, -0.5),
                                       Complex(2.0, 0.5),
                                       Complex(3.0, 0.0)};

    const std::string filename = ::testing::TempDir() + "/sternheimer_chi0_v1_test.dat";
    ModuleRI::SternheimerRPA::write_chi0_v1_file(filename, metadata, channels, chi0);

    std::ifstream in(filename.c_str(), std::ios::binary);
    ASSERT_TRUE(in.good());

    auto read_i32 = [&in]() {
        std::int32_t value = 0;
        in.read(reinterpret_cast<char*>(&value), sizeof(value));
        return value;
    };
    auto read_i64 = [&in]() {
        std::int64_t value = 0;
        in.read(reinterpret_cast<char*>(&value), sizeof(value));
        return value;
    };
    auto read_double = [&in]() {
        double value = 0.0;
        in.read(reinterpret_cast<char*>(&value), sizeof(value));
        return value;
    };

    EXPECT_EQ(read_i32(), ModuleRI::SternheimerRPA::chi0_v1_marker());
    EXPECT_EQ(read_i32(), 1);
    EXPECT_EQ(read_i32(), 2);
    EXPECT_EQ(read_i32(), 2);
    EXPECT_EQ(read_i32(), 1);
    EXPECT_EQ(read_i32(), 2);
    EXPECT_DOUBLE_EQ(read_double(), 0.5);
    EXPECT_DOUBLE_EQ(read_double(), 0.125);
    EXPECT_EQ(read_i32(), 3);
    EXPECT_EQ(read_i32(), 1);
    EXPECT_EQ(read_i32(), 1);
    EXPECT_EQ(read_i32(), 0);
    const std::int64_t first_offset = read_i64();
    EXPECT_EQ(read_i32(), 1);
    const std::int64_t second_offset = read_i64();
    EXPECT_EQ(read_i32(), 2);
    const std::int64_t third_offset = read_i64();
    EXPECT_LT(first_offset, second_offset);
    EXPECT_LT(second_offset, third_offset);

    std::remove(filename.c_str());
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

TEST(SternheimerRPA, TransitionEnergyWindowUsesOccupiedToEmptyPairsAndReturnsHartree)
{
    const std::vector<double> eigenvalues_ry = {-1.0, -0.5, 0.25, 1.5};
    const std::vector<double> occupations = {2.0, 2.0, 0.0, 0.0};

    const auto window
        = ModuleRI::SternheimerRPA::transition_energy_window_from_eigenvalues_ry(eigenvalues_ry, occupations);

    EXPECT_DOUBLE_EQ(window.emin_ha, 0.375);
    EXPECT_DOUBLE_EQ(window.emax_ha, 1.25);
}

TEST(SternheimerRPA, GreenXMinimaxFrequencyGridMatchesReference)
{
    const auto grid = ModuleRI::SternheimerRPA::generate_greenx_minimax_frequency_grid(10, 2.0, 30.0);

    const std::vector<double> expected_nodes = {0.37397426,
                                                1.17910655,
                                                2.16745299,
                                                3.50267933,
                                                5.42267468,
                                                8.30517692,
                                                12.81707194,
                                                20.35677335,
                                                34.86933737,
                                                75.79619716};
    const std::vector<double> expected_weights = {0.75729004,
                                                  0.87356340,
                                                  1.12983421,
                                                  1.57927297,
                                                  2.32153469,
                                                  3.54877702,
                                                  5.68714236,
                                                  9.94234857,
                                                  21.24030777,
                                                  80.07961331};

    ASSERT_EQ(grid.omega_ha.size(), expected_nodes.size());
    ASSERT_EQ(grid.weights_ha.size(), expected_weights.size());
    for (std::size_t i = 0; i != expected_nodes.size(); ++i)
    {
        EXPECT_NEAR(grid.omega_ha[i], expected_nodes[i], 1.0e-8);
        EXPECT_NEAR(grid.weights_ha[i], expected_weights[i], 1.0e-8);
    }
}

TEST(SternheimerRPA, ReadFrequencyGridFileAcceptsTwoAndThreeColumnRows)
{
    const std::string filename = ::testing::TempDir() + "/sternheimer_frequency_grid.dat";
    {
        std::ofstream out(filename.c_str());
        out << "# omega_Ha weight_Ha\n";
        out << "0.25 0.75\n";
        out << "2 0.50 1.25 # optional one-based index column\n";
    }

    const auto grid = ModuleRI::SternheimerRPA::read_frequency_grid_file(filename, 2);

    ASSERT_EQ(grid.omega_ha.size(), 2);
    ASSERT_EQ(grid.weights_ha.size(), 2);
    EXPECT_DOUBLE_EQ(grid.omega_ha[0], 0.25);
    EXPECT_DOUBLE_EQ(grid.weights_ha[0], 0.75);
    EXPECT_DOUBLE_EQ(grid.omega_ha[1], 0.50);
    EXPECT_DOUBLE_EQ(grid.weights_ha[1], 1.25);

    std::remove(filename.c_str());
}

TEST(SternheimerRPA, ReadFrequencyGridFileRejectsInvalidRows)
{
    const std::string count_mismatch = ::testing::TempDir() + "/sternheimer_frequency_grid_count.dat";
    {
        std::ofstream out(count_mismatch.c_str());
        out << "0.25 0.75\n";
    }
    EXPECT_THROW(ModuleRI::SternheimerRPA::read_frequency_grid_file(count_mismatch, 2), std::runtime_error);
    std::remove(count_mismatch.c_str());

    const std::string invalid_value = ::testing::TempDir() + "/sternheimer_frequency_grid_invalid.dat";
    {
        std::ofstream out(invalid_value.c_str());
        out << "-0.25 0.75\n";
    }
    EXPECT_THROW(ModuleRI::SternheimerRPA::read_frequency_grid_file(invalid_value, 1), std::invalid_argument);
    std::remove(invalid_value.c_str());
}

TEST(SternheimerRPA, FrequencyOwnerRankSupportsRankShift)
{
    EXPECT_EQ(ModuleRI::SternheimerRPA::frequency_owner_rank(0, 6, 0), 0);
    EXPECT_EQ(ModuleRI::SternheimerRPA::frequency_owner_rank(5, 6, 0), 5);
    EXPECT_EQ(ModuleRI::SternheimerRPA::frequency_owner_rank(0, 6, 1), 1);
    EXPECT_EQ(ModuleRI::SternheimerRPA::frequency_owner_rank(5, 6, 1), 0);
    EXPECT_EQ(ModuleRI::SternheimerRPA::frequency_owner_rank(0, 1, 3), 0);
    EXPECT_EQ(ModuleRI::SternheimerRPA::frequency_owner_rank(1, 4, -1), 0);
    EXPECT_THROW(ModuleRI::SternheimerRPA::frequency_owner_rank(-1, 6, 0), std::invalid_argument);
    EXPECT_THROW(ModuleRI::SternheimerRPA::frequency_owner_rank(0, 0, 0), std::invalid_argument);
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
