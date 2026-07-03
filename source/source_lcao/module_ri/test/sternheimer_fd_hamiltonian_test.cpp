#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <memory>
#include <vector>

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
