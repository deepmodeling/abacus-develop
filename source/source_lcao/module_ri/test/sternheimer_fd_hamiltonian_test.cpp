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
    Hamiltonian hamiltonian(grid, potential);

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
    block.projectors = {{Complex(1.0, 0.0), Complex(0.0, 0.0),
                         Complex(0.5, -0.25), Complex(0.0, 0.0)}};
    block.d_matrix = {{Complex(1.5, 0.0)}};
    auto nonlocal_projector = std::make_shared<ModuleRI::SternheimerFDNonlocalProjector>(
        grid.size(),
        0.5,
        std::vector<ModuleRI::SternheimerFDNonlocalProjector::ProjectorBlock>{block});
    Hamiltonian hamiltonian(grid, {0.1, -0.2, 0.3, -0.4}, 1.0, nonlocal_projector);
    const Vector psi = {Complex(0.2, -0.1), Complex(-0.4, 0.3),
                        Complex(0.7, 0.2), Complex(-0.1, -0.5)};

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

TEST(SternheimerFDHamiltonian, FiltersSiSymmetryToTheSecondOrderDiscreteStencilGroup)
{
    Hamiltonian::Grid grid{15, 15, 15, 0.0, 0.0, 0.0, true};
    constexpr double half_lattice = 5.102262569990759;
    grid.lattice_vectors = {{{0.0, half_lattice, half_lattice},
                             {half_lattice, 0.0, half_lattice},
                             {half_lattice, half_lattice, 0.0}}};
    using Rotation = ModuleRI::SternheimerFDReducedRotation;
    const std::vector<Rotation> rotations = {
        {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}},
        {{{0, -1, 0}, {1, -1, 0}, {0, -1, 1}}},
        {{{-1, 1, 0}, {-1, 0, 0}, {-1, 0, 1}}},
        {{{-1, 0, 1}, {-1, 0, 0}, {-1, 1, 0}}},
        {{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}}},
        {{{0, 0, -1}, {1, 0, -1}, {0, 1, -1}}},
        {{{0, 1, -1}, {1, 0, -1}, {0, 0, -1}}},
        {{{-1, 0, 1}, {-1, 1, 0}, {-1, 0, 0}}},
        {{{1, -1, 0}, {0, -1, 1}, {0, -1, 0}}},
        {{{0, -1, 0}, {0, -1, 1}, {1, -1, 0}}},
        {{{0, 1, -1}, {0, 0, -1}, {1, 0, -1}}},
        {{{0, 0, 1}, {0, 1, 0}, {1, 0, 0}}},
        {{{1, 0, -1}, {0, 0, -1}, {0, 1, -1}}},
        {{{0, 0, 1}, {1, 0, 0}, {0, 1, 0}}},
        {{{-1, 0, 0}, {-1, 0, 1}, {-1, 1, 0}}},
        {{{-1, 1, 0}, {-1, 0, 1}, {-1, 0, 0}}},
        {{{1, 0, -1}, {0, 1, -1}, {0, 0, -1}}},
        {{{0, -1, 1}, {1, -1, 0}, {0, -1, 0}}},
        {{{0, 0, -1}, {0, 1, -1}, {1, 0, -1}}},
        {{{0, -1, 1}, {0, -1, 0}, {1, -1, 0}}},
        {{{0, 1, 0}, {0, 0, 1}, {1, 0, 0}}},
        {{{-1, 0, 0}, {-1, 1, 0}, {-1, 0, 1}}},
        {{{0, 1, 0}, {1, 0, 0}, {0, 0, 1}}},
        {{{1, -1, 0}, {0, -1, 0}, {0, -1, 1}}},
        {{{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}},
        {{{0, 1, 0}, {-1, 1, 0}, {0, 1, -1}}},
        {{{1, -1, 0}, {1, 0, 0}, {1, 0, -1}}},
        {{{1, 0, -1}, {1, 0, 0}, {1, -1, 0}}},
        {{{-1, 0, 0}, {0, 0, -1}, {0, -1, 0}}},
        {{{0, 0, 1}, {-1, 0, 1}, {0, -1, 1}}},
        {{{0, -1, 1}, {-1, 0, 1}, {0, 0, 1}}},
        {{{1, 0, -1}, {1, -1, 0}, {1, 0, 0}}},
        {{{-1, 1, 0}, {0, 1, -1}, {0, 1, 0}}},
        {{{0, 1, 0}, {0, 1, -1}, {-1, 1, 0}}},
        {{{0, -1, 1}, {0, 0, 1}, {-1, 0, 1}}},
        {{{0, 0, -1}, {0, -1, 0}, {-1, 0, 0}}},
        {{{-1, 0, 1}, {0, 0, 1}, {0, -1, 1}}},
        {{{0, 0, -1}, {-1, 0, 0}, {0, -1, 0}}},
        {{{1, 0, 0}, {1, 0, -1}, {1, -1, 0}}},
        {{{1, -1, 0}, {1, 0, -1}, {1, 0, 0}}},
        {{{-1, 0, 1}, {0, -1, 1}, {0, 0, 1}}},
        {{{0, 1, -1}, {-1, 1, 0}, {0, 1, 0}}},
        {{{0, 0, 1}, {0, -1, 1}, {-1, 0, 1}}},
        {{{0, 1, -1}, {0, 1, 0}, {-1, 1, 0}}},
        {{{0, -1, 0}, {0, 0, -1}, {-1, 0, 0}}},
        {{{1, 0, 0}, {1, -1, 0}, {1, 0, -1}}},
        {{{0, -1, 0}, {-1, 0, 0}, {0, 0, -1}}},
        {{{-1, 1, 0}, {0, 1, 0}, {0, 1, -1}}},
    };

    EXPECT_EQ(ModuleRI::sternheimer_fd_second_order_stencil_symmetry_indices(grid, rotations),
              (std::vector<int>{0, 4, 11, 13, 20, 22, 24, 28, 35, 37, 44, 46}));
}
