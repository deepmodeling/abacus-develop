#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <limits>
#include <memory>
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
using Matrix = Hamiltonian::Matrix;

void ExpectOrderEightDenseMatrixIsHermitian(const bool periodic)
{
    Hamiltonian::Grid grid{6, 5, 4, 0.4, 0.5, 0.6, periodic};
    std::vector<double> potential(grid.size());
    for (int ir = 0; ir != grid.size(); ++ir)
    {
        potential[static_cast<std::size_t>(ir)] = 0.07 * (ir % 11) - 0.13 * (ir % 5);
    }
    Hamiltonian hamiltonian(grid, potential, 1.0, nullptr, 8);

    const auto matrix = hamiltonian.dense_matrix();
    for (std::size_t row = 0; row != matrix.size(); ++row)
    {
        for (std::size_t col = 0; col != matrix.size(); ++col)
        {
            const Complex difference = matrix[row][col] - std::conj(matrix[col][row]);
            EXPECT_NEAR(difference.real(), 0.0, 1.0e-12);
            EXPECT_NEAR(difference.imag(), 0.0, 1.0e-12);
        }
    }
}

double PlaneWaveKineticMaxError(const Hamiltonian& hamiltonian, const std::array<int, 3>& modes)
{
    const auto& grid = hamiltonian.grid();
    constexpr double two_pi = 6.283185307179586476925286766559005768;
    const std::array<double, 3> reduced_wavevector{modes[0] + grid.kpoint[0],
                                                   modes[1] + grid.kpoint[1],
                                                   modes[2] + grid.kpoint[2]};
    const auto dual = ModuleRI::sternheimer_fd_grid_dual_vectors(grid);
    std::array<double, 3> cartesian_wavevector{};
    for (int direction = 0; direction != 3; ++direction)
    {
        for (int component = 0; component != 3; ++component)
        {
            cartesian_wavevector[component] += two_pi * reduced_wavevector[direction] * dual[direction][component];
        }
    }
    double kinetic_energy = 0.0;
    for (const double component: cartesian_wavevector)
    {
        kinetic_energy += component * component;
    }
    kinetic_energy *= hamiltonian.kinetic_prefactor();

    Vector psi(static_cast<std::size_t>(grid.size()));
    for (int ix = 0; ix != grid.nx; ++ix)
    {
        for (int iy = 0; iy != grid.ny; ++iy)
        {
            for (int iz = 0; iz != grid.nz; ++iz)
            {
                const double phase = two_pi
                                     * (reduced_wavevector[0] * ix / grid.nx + reduced_wavevector[1] * iy / grid.ny
                                        + reduced_wavevector[2] * iz / grid.nz);
                psi[static_cast<std::size_t>(hamiltonian.index(ix, iy, iz))]
                    = Complex(std::cos(phase), std::sin(phase));
            }
        }
    }

    Vector kinetic_psi;
    hamiltonian.apply_kinetic(psi, kinetic_psi);
    double max_error = 0.0;
    for (std::size_t ir = 0; ir != psi.size(); ++ir)
    {
        max_error = std::max(max_error, std::abs(kinetic_psi[ir] - kinetic_energy * psi[ir]));
    }
    return max_error;
}

} // namespace

TEST(SternheimerFDHamiltonian, ConstantFunctionHasOnlyLocalPotential)
{
    Hamiltonian::Grid grid{4, 3, 2, 0.2, 0.3, 0.4, true};
    const std::vector<double> potential(grid.size(), 1.25);
    Hamiltonian hamiltonian(grid, potential);
    const Vector psi(grid.size(), Complex(2.0, -1.0));

    Vector hpsi;
    hamiltonian.apply(psi, hpsi);

    ASSERT_EQ(hpsi.size(), psi.size());
    for (std::size_t ir = 0; ir != psi.size(); ++ir)
    {
        EXPECT_NEAR(hpsi[ir].real(), (1.25 * psi[ir]).real(), 1.0e-12);
        EXPECT_NEAR(hpsi[ir].imag(), (1.25 * psi[ir]).imag(), 1.0e-12);
    }
}

TEST(SternheimerFDHamiltonian, UsesABACUSRealSpaceIndexOrder)
{
    Hamiltonian::Grid grid{2, 3, 4, 1.0, 1.0, 1.0, true};
    std::vector<double> potential(grid.size(), 0.0);
    for (int ix = 0; ix != grid.nx; ++ix)
    {
        for (int iy = 0; iy != grid.ny; ++iy)
        {
            for (int iz = 0; iz != grid.nz; ++iz)
            {
                const int abacus_index = (ix * grid.ny + iy) * grid.nz + iz;
                potential[abacus_index] = static_cast<double>(100 * ix + 10 * iy + iz);
            }
        }
    }

    Hamiltonian hamiltonian(grid, potential, 0.0);
    const Vector psi(grid.size(), Complex(1.0, 0.0));
    Vector hpsi;
    hamiltonian.apply(psi, hpsi);

    for (int ix = 0; ix != grid.nx; ++ix)
    {
        for (int iy = 0; iy != grid.ny; ++iy)
        {
            for (int iz = 0; iz != grid.nz; ++iz)
            {
                const int abacus_index = (ix * grid.ny + iy) * grid.nz + iz;
                EXPECT_EQ(hamiltonian.index(ix, iy, iz), abacus_index);
                EXPECT_DOUBLE_EQ(hpsi[abacus_index].real(), potential[abacus_index]);
                EXPECT_DOUBLE_EQ(hpsi[abacus_index].imag(), 0.0);
            }
        }
    }
}

TEST(SternheimerFDHamiltonian, ApplyUsesRequestedOpenMPThreadsAndMatchesSerial)
{
#ifndef _OPENMP
    GTEST_SKIP() << "OpenMP is not enabled in this build.";
#else
    Hamiltonian::Grid grid{24, 12, 8, 0.2, 0.3, 0.4, true};
    std::vector<double> potential(static_cast<std::size_t>(grid.size()), 0.0);
    Vector psi(static_cast<std::size_t>(grid.size()));
    for (int ir = 0; ir != grid.size(); ++ir)
    {
        potential[static_cast<std::size_t>(ir)] = 0.01 * (ir % 17);
        psi[static_cast<std::size_t>(ir)] = Complex(0.03 * (ir % 11), -0.02 * (ir % 7));
    }
    Hamiltonian hamiltonian(grid, potential, 1.0, nullptr, 8);

    const int previous_threads = omp_get_max_threads();
    omp_set_dynamic(0);
    omp_set_num_threads(1);
    int serial_threads = 0;
    Vector serial;
    hamiltonian.apply(psi, serial, &serial_threads);

    omp_set_num_threads(4);
    int parallel_threads = 0;
    Vector parallel;
    hamiltonian.apply(psi, parallel, &parallel_threads);
    omp_set_num_threads(previous_threads);

    EXPECT_EQ(serial_threads, 1);
    EXPECT_EQ(parallel_threads, 4);
    ASSERT_EQ(parallel.size(), serial.size());
    for (std::size_t ir = 0; ir != serial.size(); ++ir)
    {
        EXPECT_NEAR(parallel[ir].real(), serial[ir].real(), 1.0e-13);
        EXPECT_NEAR(parallel[ir].imag(), serial[ir].imag(), 1.0e-13);
    }
#endif
}

TEST(SternheimerFDHamiltonian, PlaneWaveHasSecondOrderFiniteDifferenceKineticEnergy)
{
    constexpr int nx = 12;
    constexpr int mode = 2;
    constexpr double length = 3.0;
    const double hx = length / nx;
    Hamiltonian::Grid grid{nx, 1, 1, hx, 1.0, 1.0, true};
    const std::vector<double> potential(grid.size(), 0.0);
    Hamiltonian hamiltonian(grid, potential);

    Vector psi(grid.size());
    constexpr double pi = 3.141592653589793238462643383279502884;
    for (int ix = 0; ix != nx; ++ix)
    {
        const double phase = 2.0 * pi * mode * ix / nx;
        psi[ix] = Complex(std::cos(phase), std::sin(phase));
    }

    Vector hpsi;
    hamiltonian.apply(psi, hpsi);
    const double expected_kinetic = (1.0 - std::cos(2.0 * pi * mode / nx)) / (hx * hx);

    for (int ix = 0; ix != nx; ++ix)
    {
        EXPECT_NEAR(hpsi[ix].real(), (expected_kinetic * psi[ix]).real(), 1.0e-12);
        EXPECT_NEAR(hpsi[ix].imag(), (expected_kinetic * psi[ix]).imag(), 1.0e-12);
    }
}

TEST(SternheimerFDHamiltonian, TwistedBoundaryHasBlochPlaneWaveKineticEnergy)
{
    constexpr int nx = 10;
    constexpr int mode = 1;
    constexpr double reduced_k = 0.25;
    constexpr double length = 4.0;
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double hx = length / nx;
    Hamiltonian::Grid grid{nx, 1, 1, hx, 1.0, 1.0, true};
    grid.kpoint = {reduced_k, 0.0, 0.0};
    Hamiltonian hamiltonian(grid, std::vector<double>(grid.size(), 0.0));

    Vector psi(grid.size());
    const double phase_step = 2.0 * pi * (mode + reduced_k) / nx;
    for (int ix = 0; ix != nx; ++ix)
    {
        const double phase = phase_step * ix;
        psi[ix] = Complex(std::cos(phase), std::sin(phase));
    }

    Vector hpsi;
    hamiltonian.apply(psi, hpsi);
    const double expected_kinetic = (1.0 - std::cos(phase_step)) / (hx * hx);
    for (int ix = 0; ix != nx; ++ix)
    {
        EXPECT_NEAR(hpsi[ix].real(), (expected_kinetic * psi[ix]).real(), 1.0e-12);
        EXPECT_NEAR(hpsi[ix].imag(), (expected_kinetic * psi[ix]).imag(), 1.0e-12);
    }
    EXPECT_EQ(hamiltonian.kpoint(), grid.kpoint);
}

TEST(SternheimerFDHamiltonian, TwistedBoundaryDenseMatrixIsHermitian)
{
    Hamiltonian::Grid grid{3, 2, 2, 0.4, 0.7, 0.9, true};
    grid.kpoint = {0.25, -0.125, 0.375};
    std::vector<double> potential(grid.size(), 0.0);
    for (std::size_t ir = 0; ir != potential.size(); ++ir)
    {
        potential[ir] = 0.1 * static_cast<double>(ir);
    }
    Hamiltonian hamiltonian(grid, potential);

    const auto matrix = hamiltonian.dense_matrix();
    for (std::size_t row = 0; row != matrix.size(); ++row)
    {
        for (std::size_t col = 0; col != matrix.size(); ++col)
        {
            const Complex difference = matrix[row][col] - std::conj(matrix[col][row]);
            EXPECT_NEAR(difference.real(), 0.0, 1.0e-12);
            EXPECT_NEAR(difference.imag(), 0.0, 1.0e-12);
        }
    }
}

TEST(SternheimerFDHamiltonian, OrderEightImprovesNonorthogonalTwistedPlaneWaveKineticEnergy)
{
    Hamiltonian::Grid grid{18, 17, 16, 1.0, 1.0, 1.0, true};
    grid.kpoint = {0.125, -0.2, 0.3};
    grid.lattice_vectors = {{{0.0, 5.1, 5.1}, {5.1, 0.0, 5.1}, {5.1, 5.1, 0.0}}};
    const std::vector<double> potential(static_cast<std::size_t>(grid.size()), 0.0);
    const Hamiltonian order_two(grid, potential, 1.0, nullptr, 2);
    const Hamiltonian order_eight(grid, potential, 1.0, nullptr, 8);

    const double order_two_error = PlaneWaveKineticMaxError(order_two, {1, -1, 2});
    const double order_eight_error = PlaneWaveKineticMaxError(order_eight, {1, -1, 2});

    EXPECT_GT(order_two_error, 1.0e-4);
    EXPECT_LT(order_eight_error, 0.02 * order_two_error);
}

TEST(SternheimerFDHamiltonian, SkipsZeroMixedDerivativesOnOrthogonalGrid)
{
    Hamiltonian::Grid grid{8, 7, 6, 0.4, 0.5, 0.6, true};
    const Hamiltonian hamiltonian(grid, std::vector<double>(grid.size(), 0.0), 1.0, nullptr, 8);

    EXPECT_EQ(hamiltonian.active_mixed_derivative_pair_count(), 0);
    EXPECT_EQ(hamiltonian.cached_first_derivative_directions(), (std::array<bool, 3>{false, false, false}));
}

TEST(SternheimerFDHamiltonian, CachesOnlyDirectionsNeededByActiveMixedDerivatives)
{
    Hamiltonian::Grid grid{8, 7, 6, 1.0, 1.0, 1.0, true};
    grid.lattice_vectors = {{{2.0, 0.0, 0.0}, {1.0, 2.0, 0.0}, {0.0, 0.0, 2.0}}};
    const Hamiltonian hamiltonian(grid, std::vector<double>(grid.size(), 0.0), 1.0, nullptr, 8);

    EXPECT_EQ(hamiltonian.active_mixed_derivative_pair_count(), 1);
    EXPECT_EQ(hamiltonian.cached_first_derivative_directions(), (std::array<bool, 3>{false, true, false}));
}

TEST(SternheimerFDHamiltonian, OrderEightNonorthogonalTwistedDenseMatrixIsHermitian)
{
    Hamiltonian::Grid grid{5, 5, 5, 1.0, 1.0, 1.0, true};
    grid.kpoint = {0.25, -0.125, 0.375};
    grid.lattice_vectors = {{{0.0, 2.3, 2.1}, {2.2, 0.0, 2.4}, {2.0, 2.5, 0.0}}};
    std::vector<double> potential(static_cast<std::size_t>(grid.size()), 0.0);
    for (std::size_t ir = 0; ir != potential.size(); ++ir)
    {
        potential[ir] = 0.03 * static_cast<double>(ir % 13);
    }
    const Hamiltonian hamiltonian(grid, potential, 1.0, nullptr, 8);

    const auto matrix = hamiltonian.dense_matrix();
    for (std::size_t row = 0; row != matrix.size(); ++row)
    {
        for (std::size_t col = 0; col != matrix.size(); ++col)
        {
            const Complex difference = matrix[row][col] - std::conj(matrix[col][row]);
            EXPECT_NEAR(difference.real(), 0.0, 1.0e-11);
            EXPECT_NEAR(difference.imag(), 0.0, 1.0e-11);
        }
    }
}

TEST(SternheimerFDHamiltonian, RejectsInvalidTwistedBoundary)
{
    Hamiltonian::Grid grid{3, 1, 1, 1.0, 1.0, 1.0, false};
    grid.kpoint = {0.25, 0.0, 0.0};
    EXPECT_THROW(Hamiltonian(grid, std::vector<double>(grid.size(), 0.0)), std::invalid_argument);

    grid.periodic = true;
    grid.kpoint = {std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0};
    EXPECT_THROW(Hamiltonian(grid, std::vector<double>(grid.size(), 0.0)), std::invalid_argument);
}

TEST(SternheimerFDHamiltonian, DenseMatrixMatchesApply)
{
    Hamiltonian::Grid grid{2, 2, 1, 0.5, 0.7, 1.0, true};
    const std::vector<double> potential = {0.0, 0.1, 0.2, 0.3};
    Hamiltonian hamiltonian(grid, potential);
    const Vector psi = {Complex(1.0, 0.5), Complex(-2.0, 1.0), Complex(0.0, -1.0), Complex(3.0, 2.0)};

    Vector hpsi;
    hamiltonian.apply(psi, hpsi);
    const auto matrix = hamiltonian.dense_matrix();

    Vector matrix_hpsi(psi.size(), Complex(0.0, 0.0));
    for (std::size_t row = 0; row != psi.size(); ++row)
    {
        for (std::size_t col = 0; col != psi.size(); ++col)
        {
            matrix_hpsi[row] += matrix[row][col] * psi[col];
        }
    }

    for (std::size_t ir = 0; ir != psi.size(); ++ir)
    {
        EXPECT_NEAR(hpsi[ir].real(), matrix_hpsi[ir].real(), 1.0e-12);
        EXPECT_NEAR(hpsi[ir].imag(), matrix_hpsi[ir].imag(), 1.0e-12);
    }
}

TEST(SternheimerFDHamiltonian, ApplyIncludesNonlocalProjector)
{
    Hamiltonian::Grid grid{2, 1, 1, 1.0, 1.0, 1.0, true};
    ModuleRI::SternheimerFDNonlocalProjector::ProjectorBlock block;
    block.projectors = {{Complex(1.0, 0.0), Complex(0.0, 0.0)}};
    block.d_matrix = {{Complex(2.0, 0.0)}};
    auto nonlocal_projector = std::make_shared<ModuleRI::SternheimerFDNonlocalProjector>(
        grid.size(),
        1.0,
        std::vector<ModuleRI::SternheimerFDNonlocalProjector::ProjectorBlock>{block});
    Hamiltonian hamiltonian(grid, std::vector<double>(grid.size(), 0.0), 0.0, nonlocal_projector);

    const Vector psi = {Complex(3.0, 0.0), Complex(4.0, 0.0)};
    Vector hpsi;
    hamiltonian.apply(psi, hpsi);

    ASSERT_EQ(hpsi.size(), psi.size());
    EXPECT_NEAR(hpsi[0].real(), 6.0, 1.0e-12);
    EXPECT_NEAR(hpsi[0].imag(), 0.0, 1.0e-12);
    EXPECT_NEAR(hpsi[1].real(), 0.0, 1.0e-12);
    EXPECT_NEAR(hpsi[1].imag(), 0.0, 1.0e-12);
    ASSERT_NE(hamiltonian.nonlocal_projector(), nullptr);
}

TEST(SternheimerFDHamiltonian, OperatorComponentsSumToFullHamiltonian)
{
    Hamiltonian::Grid grid{4, 1, 1, 0.5, 1.0, 1.0, true};
    grid.kpoint = {0.25, 0.0, 0.0};
    ModuleRI::SternheimerFDNonlocalProjector::ProjectorBlock block;
    block.projectors = {{Complex(1.0, 0.0), Complex(0.0, 0.0), Complex(0.5, -0.25), Complex(0.0, 0.0)}};
    block.d_matrix = {{Complex(1.5, 0.0)}};
    auto nonlocal_projector = std::make_shared<ModuleRI::SternheimerFDNonlocalProjector>(
        grid.size(),
        0.5,
        std::vector<ModuleRI::SternheimerFDNonlocalProjector::ProjectorBlock>{block});
    Hamiltonian hamiltonian(grid, {0.1, -0.2, 0.3, -0.4}, 1.0, nonlocal_projector);
    const Vector psi = {Complex(0.2, -0.1), Complex(-0.4, 0.3), Complex(0.7, 0.2), Complex(-0.1, -0.5)};

    Vector full;
    Vector kinetic;
    Vector local;
    Vector nonlocal;
    hamiltonian.apply(psi, full);
    hamiltonian.apply_kinetic(psi, kinetic);
    hamiltonian.apply_local_potential(psi, local);
    hamiltonian.apply_nonlocal(psi, nonlocal);

    ASSERT_EQ(full.size(), psi.size());
    for (std::size_t ir = 0; ir != full.size(); ++ir)
    {
        const Complex decomposed = kinetic[ir] + local[ir] + nonlocal[ir];
        EXPECT_NEAR(full[ir].real(), decomposed.real(), 1.0e-13);
        EXPECT_NEAR(full[ir].imag(), decomposed.imag(), 1.0e-13);
    }
}

TEST(SternheimerFDHamiltonian, ApplyBatchMatchesScalarFD8TwistedSkewNonlocal)
{
    Hamiltonian::Grid grid{5, 4, 3, 0.4, 0.5, 0.6, true};
    grid.kpoint = {0.25, -0.125, 0.375};
    grid.lattice_vectors = {{{0.2, 2.1, 1.8}, {2.3, 0.1, 2.0}, {1.9, 2.4, -0.1}}};

    std::vector<double> potential(static_cast<std::size_t>(grid.size()));
    ModuleRI::SternheimerFDNonlocalProjector::ProjectorBlock block;
    block.projectors.assign(2, Vector(static_cast<std::size_t>(grid.size())));
    for (int ir = 0; ir != grid.size(); ++ir)
    {
        potential[static_cast<std::size_t>(ir)] = 0.03 * (ir % 11) - 0.02 * (ir % 7);
        block.projectors[0][static_cast<std::size_t>(ir)]
            = Complex(0.01 * (ir % 13), -0.02 * (ir % 5));
        block.projectors[1][static_cast<std::size_t>(ir)]
            = Complex(-0.015 * (ir % 9), 0.0125 * (ir % 17));
    }
    block.d_matrix = {{Complex(0.8, 0.0), Complex(0.1, 0.2)},
                      {Complex(0.1, -0.2), Complex(-0.3, 0.0)}};
    const auto nonlocal_projector = std::make_shared<ModuleRI::SternheimerFDNonlocalProjector>(
        grid.size(),
        0.125,
        std::vector<ModuleRI::SternheimerFDNonlocalProjector::ProjectorBlock>{block});
    const Hamiltonian hamiltonian(grid, potential, 1.0, nonlocal_projector, 8);

    Matrix inputs(4, Vector(static_cast<std::size_t>(grid.size())));
    for (std::size_t column = 0; column != inputs.size(); ++column)
    {
        for (int ir = 0; ir != grid.size(); ++ir)
        {
            inputs[column][static_cast<std::size_t>(ir)]
                = Complex(0.01 * static_cast<double>((ir + 3 * column) % 19),
                          -0.0075 * static_cast<double>((2 * ir + column) % 23));
        }
    }

    Matrix expected(inputs.size());
    for (std::size_t column = 0; column != inputs.size(); ++column)
    {
        hamiltonian.apply(inputs[column], expected[column]);
    }
    Matrix actual;
    hamiltonian.apply_batch(inputs, actual);

    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t column = 0; column != expected.size(); ++column)
    {
        ASSERT_EQ(actual[column].size(), expected[column].size());
        for (std::size_t ir = 0; ir != expected[column].size(); ++ir)
        {
            EXPECT_NEAR(actual[column][ir].real(), expected[column][ir].real(), 1.0e-13);
            EXPECT_NEAR(actual[column][ir].imag(), expected[column][ir].imag(), 1.0e-13);
        }
    }

    Matrix empty_output{{Complex(1.0, 0.0)}};
    hamiltonian.apply_batch({}, empty_output);
    EXPECT_TRUE(empty_output.empty());

    Matrix one_output;
    hamiltonian.apply_batch(Matrix{inputs.front()}, one_output);
    ASSERT_EQ(one_output.size(), 1U);
    EXPECT_EQ(one_output.front(), expected.front());

    Matrix invalid_inputs = inputs;
    invalid_inputs[2].pop_back();
    EXPECT_THROW(hamiltonian.apply_batch(invalid_inputs, actual), std::invalid_argument);
}

TEST(SternheimerFDHamiltonian, DenseMatrixIsHermitianForLocalRealPotential)
{
    Hamiltonian::Grid grid{3, 2, 1, 0.4, 0.7, 1.0, true};
    const std::vector<double> potential = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5};
    Hamiltonian hamiltonian(grid, potential);

    const auto matrix = hamiltonian.dense_matrix();
    for (std::size_t row = 0; row != matrix.size(); ++row)
    {
        for (std::size_t col = 0; col != matrix.size(); ++col)
        {
            const Complex difference = matrix[row][col] - std::conj(matrix[col][row]);
            EXPECT_NEAR(difference.real(), 0.0, 1.0e-12);
            EXPECT_NEAR(difference.imag(), 0.0, 1.0e-12);
        }
    }
}

TEST(SternheimerFDHamiltonian, DenseMatrixIsHermitianForOrderEightLocalRealPotentialPeriodic)
{
    ExpectOrderEightDenseMatrixIsHermitian(true);
}

TEST(SternheimerFDHamiltonian, DenseMatrixIsHermitianForOrderEightLocalRealPotentialNonperiodic)
{
    ExpectOrderEightDenseMatrixIsHermitian(false);
}

TEST(SternheimerFDHamiltonian, DenseDiagonalizationReturnsFreeParticleEigenpairs)
{
    constexpr int nx = 4;
    constexpr double hx = 1.0;
    Hamiltonian::Grid grid{nx, 1, 1, hx, 1.0, 1.0, true};
    const std::vector<double> potential(grid.size(), 0.25);
    Hamiltonian hamiltonian(grid, potential);

    const auto eigenpairs = hamiltonian.diagonalize_dense();

    const std::vector<double> expected_eigenvalues = {0.25, 1.25, 1.25, 2.25};
    ASSERT_EQ(eigenpairs.eigenvalues.size(), expected_eigenvalues.size());
    ASSERT_EQ(eigenpairs.eigenvectors.size(), expected_eigenvalues.size());
    for (std::size_t ib = 0; ib != expected_eigenvalues.size(); ++ib)
    {
        EXPECT_NEAR(eigenpairs.eigenvalues[ib], expected_eigenvalues[ib], 1.0e-12);

        Vector hpsi;
        hamiltonian.apply(eigenpairs.eigenvectors[ib], hpsi);

        double norm = 0.0;
        for (std::size_t ir = 0; ir != hpsi.size(); ++ir)
        {
            const Complex residual = hpsi[ir] - eigenpairs.eigenvalues[ib] * eigenpairs.eigenvectors[ib][ir];
            EXPECT_NEAR(residual.real(), 0.0, 1.0e-12);
            EXPECT_NEAR(residual.imag(), 0.0, 1.0e-12);
            norm += std::norm(eigenpairs.eigenvectors[ib][ir]);
        }
        EXPECT_NEAR(norm, 1.0, 1.0e-12);
    }
}

TEST(SternheimerFDHamiltonian, KineticPrefactorScalesFiniteDifferenceLaplacian)
{
    constexpr int nx = 4;
    constexpr double hx = 1.0;
    Hamiltonian::Grid grid{nx, 1, 1, hx, 1.0, 1.0, true};
    Hamiltonian hamiltonian(grid, std::vector<double>(grid.size(), 0.0), 1.0);

    const auto eigenpairs = hamiltonian.diagonalize_dense();

    const std::vector<double> expected_eigenvalues = {0.0, 2.0, 2.0, 4.0};
    ASSERT_EQ(eigenpairs.eigenvalues.size(), expected_eigenvalues.size());
    EXPECT_NEAR(hamiltonian.kinetic_prefactor(), 1.0, 1.0e-15);
    for (std::size_t ib = 0; ib != expected_eigenvalues.size(); ++ib)
    {
        EXPECT_NEAR(eigenpairs.eigenvalues[ib], expected_eigenvalues[ib], 1.0e-12);
    }
}

TEST(SternheimerFDHamiltonian, FiltersSiSymmetryToTheFD8DiscreteStencilGroup)
{
    Hamiltonian::Grid grid{15, 15, 15, 0.0, 0.0, 0.0, true};
    constexpr double half_lattice = 5.102262569990759;
    grid.lattice_vectors
        = {{{0.0, half_lattice, half_lattice}, {half_lattice, 0.0, half_lattice}, {half_lattice, half_lattice, 0.0}}};
    using Rotation = ModuleRI::SternheimerFDReducedRotation;
    const std::vector<Rotation> rotations = {
        {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}},    {{{0, -1, 0}, {1, -1, 0}, {0, -1, 1}}},
        {{{-1, 1, 0}, {-1, 0, 0}, {-1, 0, 1}}}, {{{-1, 0, 1}, {-1, 0, 0}, {-1, 1, 0}}},
        {{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}}},    {{{0, 0, -1}, {1, 0, -1}, {0, 1, -1}}},
        {{{0, 1, -1}, {1, 0, -1}, {0, 0, -1}}}, {{{-1, 0, 1}, {-1, 1, 0}, {-1, 0, 0}}},
        {{{1, -1, 0}, {0, -1, 1}, {0, -1, 0}}}, {{{0, -1, 0}, {0, -1, 1}, {1, -1, 0}}},
        {{{0, 1, -1}, {0, 0, -1}, {1, 0, -1}}}, {{{0, 0, 1}, {0, 1, 0}, {1, 0, 0}}},
        {{{1, 0, -1}, {0, 0, -1}, {0, 1, -1}}}, {{{0, 0, 1}, {1, 0, 0}, {0, 1, 0}}},
        {{{-1, 0, 0}, {-1, 0, 1}, {-1, 1, 0}}}, {{{-1, 1, 0}, {-1, 0, 1}, {-1, 0, 0}}},
        {{{1, 0, -1}, {0, 1, -1}, {0, 0, -1}}}, {{{0, -1, 1}, {1, -1, 0}, {0, -1, 0}}},
        {{{0, 0, -1}, {0, 1, -1}, {1, 0, -1}}}, {{{0, -1, 1}, {0, -1, 0}, {1, -1, 0}}},
        {{{0, 1, 0}, {0, 0, 1}, {1, 0, 0}}},    {{{-1, 0, 0}, {-1, 1, 0}, {-1, 0, 1}}},
        {{{0, 1, 0}, {1, 0, 0}, {0, 0, 1}}},    {{{1, -1, 0}, {0, -1, 0}, {0, -1, 1}}},
        {{{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}}, {{{0, 1, 0}, {-1, 1, 0}, {0, 1, -1}}},
        {{{1, -1, 0}, {1, 0, 0}, {1, 0, -1}}},  {{{1, 0, -1}, {1, 0, 0}, {1, -1, 0}}},
        {{{-1, 0, 0}, {0, 0, -1}, {0, -1, 0}}}, {{{0, 0, 1}, {-1, 0, 1}, {0, -1, 1}}},
        {{{0, -1, 1}, {-1, 0, 1}, {0, 0, 1}}},  {{{1, 0, -1}, {1, -1, 0}, {1, 0, 0}}},
        {{{-1, 1, 0}, {0, 1, -1}, {0, 1, 0}}},  {{{0, 1, 0}, {0, 1, -1}, {-1, 1, 0}}},
        {{{0, -1, 1}, {0, 0, 1}, {-1, 0, 1}}},  {{{0, 0, -1}, {0, -1, 0}, {-1, 0, 0}}},
        {{{-1, 0, 1}, {0, 0, 1}, {0, -1, 1}}},  {{{0, 0, -1}, {-1, 0, 0}, {0, -1, 0}}},
        {{{1, 0, 0}, {1, 0, -1}, {1, -1, 0}}},  {{{1, -1, 0}, {1, 0, -1}, {1, 0, 0}}},
        {{{-1, 0, 1}, {0, -1, 1}, {0, 0, 1}}},  {{{0, 1, -1}, {-1, 1, 0}, {0, 1, 0}}},
        {{{0, 0, 1}, {0, -1, 1}, {-1, 0, 1}}},  {{{0, 1, -1}, {0, 1, 0}, {-1, 1, 0}}},
        {{{0, -1, 0}, {0, 0, -1}, {-1, 0, 0}}}, {{{1, 0, 0}, {1, -1, 0}, {1, 0, -1}}},
        {{{0, -1, 0}, {-1, 0, 0}, {0, 0, -1}}}, {{{-1, 1, 0}, {0, 1, 0}, {0, 1, -1}}},
    };

    EXPECT_EQ(ModuleRI::sternheimer_fd_stencil_symmetry_indices(grid, 8, rotations),
              (std::vector<int>{0, 4, 11, 13, 20, 22, 24, 28, 35, 37, 44, 46}));
}

TEST(SternheimerFDHamiltonian, RequiresExactPeriodicGridPermutationForSpatialSymmetry)
{
    using Rotation = ModuleRI::SternheimerFDReducedRotation;
    const Hamiltonian::Grid square{4, 4, 3, 1.0, 1.0, 1.0, true};
    const Rotation identity = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    const Rotation swap_xy = {{{0, 1, 0}, {1, 0, 0}, {0, 0, 1}}};
    const Rotation invert_x = {{{-1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};

    EXPECT_TRUE(ModuleRI::sternheimer_fd_grid_rotation_is_permutation(square, identity));
    EXPECT_TRUE(ModuleRI::sternheimer_fd_grid_rotation_is_permutation(square, swap_xy));
    EXPECT_TRUE(ModuleRI::sternheimer_fd_grid_rotation_is_permutation(square, invert_x));

    const Hamiltonian::Grid rectangular{4, 5, 3, 1.0, 1.0, 1.0, true};
    EXPECT_FALSE(ModuleRI::sternheimer_fd_grid_rotation_is_permutation(rectangular, swap_xy));
}

TEST(SternheimerFDHamiltonian, RejectsHexagonalC3WhenTheFD8StencilIsNotCovariant)
{
    using Rotation = ModuleRI::SternheimerFDReducedRotation;
    Hamiltonian::Grid grid{24, 24, 48, 0.0, 0.0, 0.0, true};
    constexpr double lattice_a = 2.5;
    grid.lattice_vectors = {{{lattice_a, 0.0, 0.0},
                             {-0.5 * lattice_a, std::sqrt(3.0) * 0.5 * lattice_a, 0.0},
                             {0.0, 0.0, 20.0}}};
    const Rotation identity = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    const Rotation c3 = {{{-1, -1, 0}, {1, 0, 0}, {0, 0, 1}}};

    // The C3 maps the periodic tensor grid exactly, but the mixed-derivative
    // finite-difference stencil assigns different weights to its mapped offsets.
    EXPECT_TRUE(ModuleRI::sternheimer_fd_grid_rotation_is_permutation(grid, c3));
    EXPECT_EQ(ModuleRI::sternheimer_fd_stencil_symmetry_indices(grid, 8, {identity, c3}),
              (std::vector<int>{0}));
}

TEST(SternheimerFDHamiltonian, FiltersStencilSymmetryWhenGridPermutationFails)
{
    using Rotation = ModuleRI::SternheimerFDReducedRotation;
    const Hamiltonian::Grid rectangular{4, 5, 3, 1.0, 1.0, 1.0, true};
    const Rotation identity = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    const Rotation swap_xy = {{{0, 1, 0}, {1, 0, 0}, {0, 0, 1}}};

    EXPECT_EQ(ModuleRI::sternheimer_fd_second_order_stencil_symmetry_indices(
                  rectangular, {identity, swap_xy}),
              (std::vector<int>{0}));
}
