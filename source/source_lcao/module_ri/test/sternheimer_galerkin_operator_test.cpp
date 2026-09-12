#include "source_lcao/module_ri/sternheimer_galerkin_operator.h"

#include "source_base/module_external/lapack_connector.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace
{

using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
using Galerkin = ModuleRI::SternheimerGalerkinOperator;
using Projector = ModuleRI::SternheimerFDNonlocalProjector;
using Grid = Hamiltonian::Grid;
using Complex = Hamiltonian::Complex;
using Vector = Hamiltonian::Vector;
using Matrix = Hamiltonian::Matrix;
constexpr double two_pi = 6.283185307179586476925286766559;

enum class Term { Total, Kinetic, Local, Nonlocal };

Grid MakeGrid(const int nx, const int ny, const int nz)
{
    Grid grid{nx, ny, nz, 4.0 / nx, std::sqrt(12.74) / ny, std::sqrt(9.2) / nz, true};
    grid.lattice_vectors = {{{4.0, 0.0, 0.0}, {0.7, 3.5, 0.0}, {0.2, 0.4, 3.0}}};
    grid.kpoint = {0.21, -0.17, 0.13};
    return grid;
}

std::array<double, 3> Position(const Grid& grid, const int index)
{
    return {{static_cast<double>(index / (grid.ny * grid.nz)) / grid.nx,
             static_cast<double>((index / grid.nz) % grid.ny) / grid.ny,
             static_cast<double>(index % grid.nz) / grid.nz}};
}

double VolumeElement(const Grid& grid)
{
    const auto lattice = ModuleRI::sternheimer_fd_grid_lattice_vectors(grid);
    const auto& a = lattice[0];
    const auto& b = lattice[1];
    const auto& c = lattice[2];
    return std::abs(a[0] * (b[1] * c[2] - b[2] * c[1])
                    - a[1] * (b[0] * c[2] - b[2] * c[0])
                    + a[2] * (b[0] * c[1] - b[1] * c[0])) / grid.size();
}

Complex Inner(const Vector& left, const Vector& right)
{
    Complex value{};
    for (std::size_t i = 0; i != left.size(); ++i)
    {
        value += std::conj(left[i]) * right[i];
    }
    return value;
}

Vector Multiply(const Matrix& matrix, const Vector& vector)
{
    Vector result(matrix.size());
    for (std::size_t row = 0; row != matrix.size(); ++row)
    {
        for (std::size_t col = 0; col != vector.size(); ++col)
        {
            result[row] += matrix[row][col] * vector[col];
        }
    }
    return result;
}

void ExpectNear(const Vector& actual, const Vector& expected, const double tolerance = 5.0e-11)
{
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i != actual.size(); ++i)
    {
        EXPECT_NEAR(std::abs(actual[i] - expected[i]), 0.0, tolerance) << "index " << i;
    }
}

Vector Probe(const Grid& grid, const double offset = 0.0)
{
    Vector vector(grid.size());
    for (int i = 0; i != grid.size(); ++i)
    {
        vector[i] = Complex(std::cos(0.37 * i + offset), std::sin(0.23 * i - offset));
    }
    return vector;
}

void ApplyFine(const Hamiltonian& fine, const Term term, const Vector& input, Vector& output)
{
    switch (term)
    {
    case Term::Total: fine.apply(input, output); break;
    case Term::Kinetic: fine.apply_kinetic(input, output); break;
    case Term::Local: fine.apply_local_potential(input, output); break;
    case Term::Nonlocal: fine.apply_nonlocal(input, output); break;
    }
}

void ApplyCoarse(Galerkin& op, const Term term, const Vector& input, Vector& output)
{
    switch (term)
    {
    case Term::Total: op.apply(input, output); break;
    case Term::Kinetic: op.apply_kinetic(input, output); break;
    case Term::Local: op.apply_local_potential(input, output); break;
    case Term::Nonlocal: op.apply_nonlocal(input, output); break;
    }
}

void ApplyBatch(Galerkin& op, const Term term, const Matrix& input, Matrix& output)
{
    switch (term)
    {
    case Term::Total: op.apply_batch(input, output); break;
    case Term::Kinetic: op.apply_kinetic_batch(input, output); break;
    case Term::Local: op.apply_local_potential_batch(input, output); break;
    case Term::Nonlocal: op.apply_nonlocal_batch(input, output); break;
    }
}

// Direct Fourier sums are independent of the FFT transfer implementation.
Matrix ExplicitInterpolation(const Grid& coarse, const Grid& fine)
{
    Matrix interpolation(fine.size(), Vector(coarse.size()));
    const std::array<int, 3> dims{{coarse.nx, coarse.ny, coarse.nz}};
    for (int row = 0; row != fine.size(); ++row)
    {
        const auto rf = Position(fine, row);
        for (int col = 0; col != coarse.size(); ++col)
        {
            const auto rc = Position(coarse, col);
            for (int mode = 0; mode != coarse.size(); ++mode)
            {
                const std::array<int, 3> indices{{mode / (coarse.ny * coarse.nz),
                                                  (mode / coarse.nz) % coarse.ny, mode % coarse.nz}};
                double phase = 0.0;
                for (int axis = 0; axis != 3; ++axis)
                {
                    const int signed_mode = indices[axis] <= dims[axis] / 2
                                                ? indices[axis] : indices[axis] - dims[axis];
                    phase += (signed_mode + coarse.kpoint[axis]) * (rf[axis] - rc[axis]);
                }
                interpolation[row][col] += std::polar(1.0 / coarse.size(), two_pi * phase);
            }
        }
    }
    return interpolation;
}

Vector ExplicitAdjoint(const Matrix& interpolation, const Grid& coarse, const Grid& fine, const Vector& vector)
{
    Vector result(coarse.size());
    const double weight = VolumeElement(fine) / VolumeElement(coarse);
    for (int col = 0; col != coarse.size(); ++col)
    {
        for (int row = 0; row != fine.size(); ++row)
        {
            result[col] += weight * std::conj(interpolation[row][col]) * vector[row];
        }
    }
    return result;
}

Matrix ExplicitGalerkin(const Hamiltonian& fine, const Grid& coarse, const Term term)
{
    const auto interpolation = ExplicitInterpolation(coarse, fine.grid());
    Matrix fine_matrix(fine.grid().size(), Vector(fine.grid().size()));
    Vector basis(fine.grid().size());
    Vector image;
    for (int col = 0; col != fine.grid().size(); ++col)
    {
        std::fill(basis.begin(), basis.end(), Complex{});
        basis[col] = 1.0;
        ApplyFine(fine, term, basis, image);
        for (int row = 0; row != fine.grid().size(); ++row)
        {
            fine_matrix[row][col] = image[row];
        }
    }
    Matrix result(coarse.size(), Vector(coarse.size()));
    for (int col = 0; col != coarse.size(); ++col)
    {
        Vector lifted(fine.grid().size());
        for (int row = 0; row != fine.grid().size(); ++row)
        {
            lifted[row] = interpolation[row][col];
        }
        const auto restricted = ExplicitAdjoint(interpolation, coarse, fine.grid(), Multiply(fine_matrix, lifted));
        for (int row = 0; row != coarse.size(); ++row)
        {
            result[row][col] = restricted[row];
        }
    }
    return result;
}

std::shared_ptr<const Projector> MakeProjector(const Grid& grid)
{
    Projector::ProjectorBlock block;
    block.projectors.assign(2, Vector(grid.size()));
    block.d_matrix = {{Complex(0.6, 0.0), Complex(0.1, 0.05)},
                      {Complex(0.1, -0.05), Complex(0.4, 0.0)}};
    for (int i = 0; i != grid.size(); ++i)
    {
        const auto r = Position(grid, i);
        const Complex phase = std::polar(1.0, two_pi * (grid.kpoint[0] * r[0]
                                                       + grid.kpoint[1] * r[1] + grid.kpoint[2] * r[2]));
        block.projectors[0][i] = phase * (0.3 + 0.2 * std::cos(two_pi * r[0]));
        block.projectors[1][i] = phase * Complex(0.2 * std::sin(two_pi * r[1]), 0.15 * std::cos(two_pi * r[2]));
    }
    return std::make_shared<const Projector>(grid.size(), VolumeElement(grid),
                                              std::vector<Projector::ProjectorBlock>{block});
}

std::shared_ptr<const Hamiltonian> MakeFine(const Grid& grid)
{
    std::vector<double> potential(grid.size());
    for (int i = 0; i != grid.size(); ++i)
    {
        const auto r = Position(grid, i);
        potential[i] = 1.2 + 0.2 * std::cos(two_pi * r[0]) + 0.1 * std::sin(two_pi * (r[1] + r[2]));
    }
    return std::make_shared<const Hamiltonian>(grid, potential, 1.0, MakeProjector(grid), 8);
}

Matrix MatrixProduct(const Matrix& left, const Matrix& right)
{
    Matrix result(left.size(), Vector(right.front().size()));
    for (std::size_t row = 0; row != left.size(); ++row)
    {
        for (std::size_t col = 0; col != right.front().size(); ++col)
        {
            for (std::size_t inner = 0; inner != right.size(); ++inner)
            {
                result[row][col] += left[row][inner] * right[inner][col];
            }
        }
    }
    return result;
}

Matrix DirectInverse(const Matrix& matrix)
{
    const int n = static_cast<int>(matrix.size());
    Vector a(n * n);
    Vector inverse(n * n);
    for (int col = 0; col != n; ++col)
    {
        inverse[col + n * col] = 1.0;
        for (int row = 0; row != n; ++row) a[row + n * col] = matrix[row][col];
    }
    std::vector<int> pivots(n);
    int info = 0;
    zgesv_(&n, &n, a.data(), &n, pivots.data(), inverse.data(), &n, &info);
    if (info != 0) throw std::runtime_error("Test direct complex inverse failed.");
    Matrix result(n, Vector(n));
    for (int row = 0; row != n; ++row)
    {
        for (int col = 0; col != n; ++col) result[row][col] = inverse[row + n * col];
    }
    return result;
}

std::vector<double> HermitianEigenvalues(const Matrix& matrix)
{
    const int n = static_cast<int>(matrix.size());
    Vector a(n * n);
    for (int col = 0; col != n; ++col)
    {
        for (int row = 0; row != n; ++row) a[row + n * col] = matrix[row][col];
    }
    const char jobz = 'N';
    const char uplo = 'U';
    const int lwork = std::max(1, 2 * n - 1);
    Vector work(lwork);
    std::vector<double> rwork(std::max(1, 3 * n - 2));
    std::vector<double> values(n);
    int info = 0;
    zheev_(&jobz, &uplo, &n, a.data(), &n, values.data(), work.data(), &lwork, rwork.data(), &info);
    if (info != 0) throw std::runtime_error("Test Hermitian diagonalization failed.");
    return values;
}

Matrix CoarseMatrix(Galerkin& op)
{
    const int n = op.grid().size();
    Matrix matrix(n, Vector(n));
    for (int col = 0; col != n; ++col)
    {
        Vector basis(n);
        Vector image;
        basis[col] = 1.0;
        op.apply(basis, image);
        for (int row = 0; row != n; ++row) matrix[row][col] = image[row];
    }
    return matrix;
}

Vector BlochState(const Grid& grid, const int xmode = 0)
{
    Vector state(grid.size());
    const double amplitude = 1.0 / std::sqrt(VolumeElement(grid) * grid.size());
    for (int i = 0; i != grid.size(); ++i)
    {
        const auto r = Position(grid, i);
        state[i] = std::polar(amplitude, two_pi * ((xmode + grid.kpoint[0]) * r[0]
                                                   + grid.kpoint[1] * r[1] + grid.kpoint[2] * r[2]));
    }
    return state;
}

Matrix VirtualProjector(const Vector& occupied, const double dv)
{
    Matrix pc(occupied.size(), Vector(occupied.size()));
    for (std::size_t row = 0; row != occupied.size(); ++row)
    {
        for (std::size_t col = 0; col != occupied.size(); ++col)
        {
            pc[row][col] = (row == col ? 1.0 : 0.0) - dv * occupied[row] * std::conj(occupied[col]);
        }
    }
    return pc;
}

Matrix ResponseOperatorWithOccupiedFill(Matrix h, const Matrix& pc, const double occupied_energy, const double occupied_fill)
{
    for (std::size_t i = 0; i != h.size(); ++i) h[i][i] -= occupied_energy;
    auto a = MatrixProduct(MatrixProduct(pc, h), pc);
    // Only fill the unused occupied block. The virtual reference differences
    // must already be positive; no excited-state energy adjustment is allowed.
    for (std::size_t row = 0; row != a.size(); ++row)
    {
        for (std::size_t col = 0; col != a.size(); ++col)
        {
            a[row][col] += occupied_fill * ((row == col ? 1.0 : 0.0) - pc[row][col]);
        }
    }
    return a;
}

void ExpectUnshiftedVirtualSpectrumPositive(Matrix h,
                                           const Matrix& pc,
                                           const Vector& occupied,
                                           const double dv,
                                           const double occupied_energy)
{
    EXPECT_NEAR(dv * Inner(occupied, occupied).real(), 1.0, 1.0e-12);
    ExpectNear(Multiply(pc, occupied), Vector(occupied.size()), 1.0e-12);
    const auto pc2 = MatrixProduct(pc, pc);
    for (std::size_t row = 0; row != pc.size(); ++row)
    {
        ExpectNear(pc2[row], pc[row], 1.0e-12);
        for (std::size_t col = 0; col != pc.size(); ++col)
        {
            EXPECT_NEAR(std::abs(pc[row][col] - std::conj(pc[col][row])), 0.0, 1.0e-12);
        }
    }
    for (std::size_t i = 0; i != h.size(); ++i) h[i][i] -= occupied_energy;
    const auto unshifted = MatrixProduct(MatrixProduct(pc, h), pc);
    const auto eigenvalues = HermitianEigenvalues(unshifted);
    // These fixtures have one occupied state: exactly one null mode is allowed,
    // and all remaining (virtual) reference differences must be positive.
    ASSERT_GE(eigenvalues.size(), 2U);
    EXPECT_NEAR(eigenvalues.front(), 0.0, 1.0e-10);
    EXPECT_GT(eigenvalues[1], 1.0e-8);
}

Matrix FinePerturbationProducts(const Grid& fine, const Vector& source)
{
    Matrix products(3, Vector(fine.size()));
    for (int i = 0; i != fine.size(); ++i)
    {
        const auto r = Position(fine, i);
        products[0][i] = Complex(std::cos(two_pi * r[0]), 0.2 * std::sin(two_pi * r[1])) * source[i];
        products[1][i] = Complex(0.3 * std::cos(two_pi * r[0]) + std::sin(two_pi * r[1]),
                                  0.1 * std::cos(two_pi * r[2])) * source[i];
        products[2][i] = Complex(0.2 + 0.1 * std::cos(two_pi * (r[0] + r[1])),
                                  0.4 * std::sin(two_pi * r[0])) * source[i];
    }
    return products;
}

Matrix DirectResponse(const Matrix& a, const Matrix& rhs_vectors, const double dv, const double omega)
{
    Matrix plus = a;
    Matrix minus = a;
    for (std::size_t i = 0; i != a.size(); ++i)
    {
        plus[i][i] += Complex(0.0, omega);
        minus[i][i] -= Complex(0.0, omega);
    }
    const auto plus_inverse = DirectInverse(plus);
    const auto minus_inverse = DirectInverse(minus);
    Matrix chi(rhs_vectors.size(), Vector(rhs_vectors.size()));
    // Synthetic two-sided response: identical projected vertices on left/right,
    // chi = -B^dagger [(A+i*w)^-1 + (A-i*w)^-1] B, with the grid inner product.
    // This tests the FullGrid FD8 reference construction, not old weak-A or RPA acceptance.
    for (std::size_t col = 0; col != rhs_vectors.size(); ++col)
    {
        const auto xp = Multiply(plus_inverse, rhs_vectors[col]);
        const auto xm = Multiply(minus_inverse, rhs_vectors[col]);
        ExpectNear(Multiply(plus, xp), rhs_vectors[col], 1.0e-11);
        ExpectNear(Multiply(minus, xm), rhs_vectors[col], 1.0e-11);
        for (std::size_t row = 0; row != rhs_vectors.size(); ++row)
        {
            chi[row][col] = -dv * (Inner(rhs_vectors[row], xp) + Inner(rhs_vectors[row], xm));
        }
    }
    return chi;
}

} // namespace

TEST(SternheimerGalerkinOperator, EveryComponentMatchesExplicitFullMatrixAndIsHermitian)
{
    const auto coarse = MakeGrid(3, 2, 2);
    const auto fine = MakeFine(MakeGrid(6, 5, 4));
    Galerkin op(fine, coarse);
    EXPECT_EQ(op.fine_hamiltonian().finite_difference_order(), 8);
    EXPECT_EQ(&op.fine_hamiltonian(), fine.get());
    ASSERT_GT(fine->active_mixed_derivative_pair_count(), 0);
    Matrix total(coarse.size(), Vector(coarse.size()));
    Matrix sum(coarse.size(), Vector(coarse.size()));
    for (const auto term: {Term::Total, Term::Kinetic, Term::Local, Term::Nonlocal})
    {
        const auto expected = ExplicitGalerkin(*fine, coarse, term);
        Matrix actual(coarse.size(), Vector(coarse.size()));
        for (int col = 0; col != coarse.size(); ++col)
        {
            Vector basis(coarse.size());
            basis[col] = 1.0;
            Vector image;
            ApplyCoarse(op, term, basis, image);
            ASSERT_EQ(image.size(), basis.size());
            for (int row = 0; row != coarse.size(); ++row)
            {
                actual[row][col] = image[row];
                if (term == Term::Total) total[row][col] = image[row];
                else sum[row][col] += image[row];
            }
        }
        for (int row = 0; row != coarse.size(); ++row)
        {
            ExpectNear(actual[row], expected[row]);
            for (int col = 0; col != coarse.size(); ++col)
            {
                EXPECT_NEAR(std::abs(actual[row][col] - std::conj(actual[col][row])), 0.0, 5.0e-11);
            }
        }
    }
    double imaginary_max = 0.0;
    for (int row = 0; row != coarse.size(); ++row)
    {
        ExpectNear(total[row], sum[row]);
        for (const auto value: total[row]) imaginary_max = std::max(imaginary_max, std::abs(value.imag()));
    }
    EXPECT_GT(imaginary_max, 1.0e-3);
}

TEST(SternheimerGalerkinOperator, PositiveFineOperatorBoundsCoarseRayleighQuotients)
{
    const auto coarse = MakeGrid(3, 2, 2);
    const auto fine = MakeFine(MakeGrid(6, 5, 4));
    const auto eigenpairs = fine->diagonalize_dense();
    ASSERT_GT(eigenpairs.eigenvalues.front(), 0.0);
    const auto interpolation = ExplicitInterpolation(coarse, fine->grid());
    Galerkin op(fine, coarse);
    for (int trial = 0; trial != 5; ++trial)
    {
        const auto psi = Probe(coarse, 0.3 * trial);
        const auto lifted = Multiply(interpolation, psi);
        Vector coarse_image;
        Vector fine_image;
        op.apply(psi, coarse_image);
        fine->apply(lifted, fine_image);
        const double coarse_norm = Inner(psi, psi).real();
        const double fine_norm = Inner(lifted, lifted).real();
        const Complex rayleigh = Inner(psi, coarse_image) / coarse_norm;
        EXPECT_NEAR(coarse_norm * VolumeElement(coarse), fine_norm * VolumeElement(fine->grid()), 1.0e-10);
        EXPECT_NEAR(rayleigh.imag(), 0.0, 1.0e-11);
        EXPECT_NEAR(rayleigh.real(), (Inner(lifted, fine_image) / fine_norm).real(), 1.0e-10);
        EXPECT_GE(rayleigh.real(), eigenpairs.eigenvalues.front() - 1.0e-10);
        EXPECT_LE(rayleigh.real(), eigenpairs.eigenvalues.back() + 1.0e-10);
        op.apply_kinetic(psi, coarse_image);
        EXPECT_GE(Inner(psi, coarse_image).real(), -1.0e-10);
    }
}

TEST(SternheimerGalerkinOperator, ConstantFinePotentialShiftIsExactlyCoarseIdentityShift)
{
    const auto coarse = MakeGrid(3, 2, 2);
    const auto fine = MakeFine(MakeGrid(6, 5, 4));
    auto potential = fine->local_potential();
    constexpr double shift = 1.73;
    for (auto& value: potential) value += shift;
    const Hamiltonian shifted(fine->grid(), potential, 1.0, MakeProjector(fine->grid()), 8);
    Galerkin first(fine, coarse);
    Galerkin second(shifted, coarse);
    const auto psi = Probe(coarse);
    for (const auto term: {Term::Total, Term::Kinetic, Term::Local, Term::Nonlocal})
    {
        Vector before;
        Vector after;
        ApplyCoarse(first, term, psi, before);
        ApplyCoarse(second, term, psi, after);
        for (std::size_t i = 0; i != psi.size(); ++i)
        {
            const auto expected = before[i] + ((term == Term::Total || term == Term::Local) ? shift * psi[i] : Complex{});
            EXPECT_NEAR(std::abs(after[i] - expected), 0.0, 5.0e-11);
        }
    }
}

TEST(SternheimerGalerkinOperator, NonlocalRestrictionKeepsOriginalDAndDoesNotNormalizeProjectors)
{
    const auto coarse = MakeGrid(3, 2, 2);
    const auto fine = MakeFine(MakeGrid(6, 5, 4));
    const auto original_blocks = fine->nonlocal_projector()->blocks();
    auto coarse_blocks = original_blocks;
    const auto interpolation = ExplicitInterpolation(coarse, fine->grid());
    for (auto& block: coarse_blocks)
    {
        for (auto& beta: block.projectors)
        {
            beta = ExplicitAdjoint(interpolation, coarse, fine->grid(), beta);
        }
    }
    const Projector coarse_projector(coarse.size(), VolumeElement(coarse), coarse_blocks);
    Galerkin op(fine, coarse);
    const auto psi = Probe(coarse);
    Vector actual;
    Vector expected;
    op.apply_nonlocal(psi, actual);
    coarse_projector.apply(psi, expected);
    ExpectNear(actual, expected);
    ASSERT_EQ(fine->nonlocal_projector()->blocks().size(), original_blocks.size());
    for (std::size_t i = 0; i != original_blocks.size(); ++i)
    {
        EXPECT_EQ(fine->nonlocal_projector()->blocks()[i].d_matrix, original_blocks[i].d_matrix);
        EXPECT_EQ(fine->nonlocal_projector()->blocks()[i].projectors, original_blocks[i].projectors);
    }
}

TEST(SternheimerGalerkinOperator, PointSampledLocalPotentialAliasesButGalerkinDoesNot)
{
    const auto coarse = MakeGrid(3, 1, 1);
    const auto fine_grid = MakeGrid(9, 1, 1);
    std::vector<double> potential(fine_grid.size());
    for (int i = 0; i != fine_grid.size(); ++i)
    {
        potential[i] = std::cos(two_pi * 3.0 * i / fine_grid.nx);
    }
    const Hamiltonian fine(fine_grid, potential, 1.0, nullptr, 8);
    Galerkin op(fine, coarse);
    const auto psi = Probe(coarse);
    Vector actual;
    op.apply_local_potential(psi, actual);
    ExpectNear(actual, Vector(coarse.size()), 1.0e-12);
    // Every coarse point samples cos(6*pi*x) as +1 although P V_f P is zero.
    const Hamiltonian sampled(coarse, std::vector<double>(coarse.size(), 1.0), 1.0, nullptr, 8);
    Vector aliased;
    sampled.apply_local_potential(psi, aliased);
    ExpectNear(aliased, psi);
    EXPECT_GT(std::abs(Inner(psi, aliased) - Inner(psi, actual)), 1.0);
}

TEST(SternheimerGalerkinOperator, SameGridExactlyMatchesFullGridFineFD8NotHybridResponse)
{
    const auto fine = MakeFine(MakeGrid(3, 2, 2));
    Galerkin op(fine, fine->grid());
    for (const auto term: {Term::Total, Term::Kinetic, Term::Local, Term::Nonlocal})
    {
        auto psi = Probe(fine->grid());
        Vector expected;
        ApplyFine(*fine, term, psi, expected);
        ApplyCoarse(op, term, psi, psi);
        EXPECT_EQ(psi, expected);
        const Matrix input{Probe(fine->grid()), Probe(fine->grid(), 0.4)};
        Matrix output;
        ApplyBatch(op, term, input, output);
        ASSERT_EQ(output.size(), input.size());
        for (std::size_t i = 0; i != input.size(); ++i)
        {
            ApplyFine(*fine, term, input[i], expected);
            EXPECT_EQ(output[i], expected);
        }
        Matrix empty;
        ApplyBatch(op, term, empty, empty);
        EXPECT_TRUE(empty.empty());
    }
}

TEST(SternheimerGalerkinOperator, SeparateWorkersApplyConcurrentlyWithSharedImmutableFineHamiltonian)
{
#ifndef _OPENMP
    GTEST_SKIP() << "OpenMP is not enabled in this build.";
#else
    const auto coarse = MakeGrid(3, 2, 2);
    const auto fine = MakeFine(MakeGrid(6, 5, 4));
    std::array<std::unique_ptr<Galerkin>, 2> workers;
    for (auto& worker: workers) worker.reset(new Galerkin(fine, coarse));
    const Matrix inputs{Probe(coarse), Probe(coarse, 0.4)};
    for (const auto term: {Term::Total, Term::Kinetic, Term::Local, Term::Nonlocal})
    {
        Matrix expected(2);
        Matrix actual(2);
        for (int worker = 0; worker != 2; ++worker)
        {
            ApplyCoarse(*workers[worker], term, inputs[worker], expected[worker]);
        }
#pragma omp parallel for num_threads(2) schedule(static)
        for (int worker = 0; worker < 2; ++worker)
        {
            for (int repeat = 0; repeat != 4; ++repeat)
            {
                ApplyCoarse(*workers[worker], term, inputs[worker], actual[worker]);
            }
        }
        for (int worker = 0; worker != 2; ++worker) ExpectNear(actual[worker], expected[worker]);
    }
#endif
}

TEST(SternheimerGalerkinOperator, OwnsFineLifetimeAndPreservesZeroNonlocalComponent)
{
    const auto coarse = MakeGrid(3, 2, 2);
    const auto fine_grid = MakeGrid(6, 5, 4);
    auto fine = MakeFine(fine_grid);
    const std::weak_ptr<const Hamiltonian> weak = fine;
    Galerkin shared(fine, coarse);
    fine.reset();
    EXPECT_FALSE(weak.expired());
    Vector output;
    EXPECT_NO_THROW(shared.apply(Probe(coarse), output));
    Galerkin copied(Hamiltonian(fine_grid, std::vector<double>(fine_grid.size(), 1.0), 1.0, nullptr, 8), coarse);
    EXPECT_EQ(copied.fine_hamiltonian().finite_difference_order(), 8);
    EXPECT_EQ(copied.fine_hamiltonian().nonlocal_projector(), nullptr);
    copied.apply_nonlocal(Probe(coarse), output);
    EXPECT_EQ(output, Vector(coarse.size()));
}

TEST(SternheimerGalerkinOperator, StreamedComponentBatchesMatchScalarAndUseIndependentScratch)
{
    const auto coarse = MakeGrid(3, 2, 2);
    const auto fine = MakeFine(MakeGrid(6, 5, 4));
    Galerkin first(fine, coarse);
    Galerkin second(fine, coarse);
    const auto workspace = first.workspace_bytes();
    const Matrix inputs{Probe(coarse), Probe(coarse, 0.4), Probe(coarse, 0.8)};
    for (const auto term: {Term::Total, Term::Kinetic, Term::Local, Term::Nonlocal})
    {
        Matrix outputs;
        ApplyBatch(first, term, inputs, outputs);
        ASSERT_EQ(outputs.size(), inputs.size());
        for (std::size_t i = 0; i != inputs.size(); ++i)
        {
            Vector expected;
            ApplyCoarse(second, term, inputs[i], expected);
            ExpectNear(outputs[i], expected);
        }
        Matrix in_place = inputs;
        ApplyBatch(first, term, in_place, in_place);
        for (std::size_t i = 0; i != inputs.size(); ++i) ExpectNear(in_place[i], outputs[i]);
        EXPECT_EQ(first.workspace_bytes(), workspace);
    }
}

TEST(SternheimerGalerkinOperator, RejectsMalformedInputsBeforeChangingOutputs)
{
    const auto coarse = MakeGrid(3, 2, 2);
    Galerkin op(MakeFine(MakeGrid(6, 5, 4)), coarse);
    const Vector sentinel{Complex(17.0, -2.0)};
    for (const auto term: {Term::Total, Term::Kinetic, Term::Local, Term::Nonlocal})
    {
        Vector output = sentinel;
        EXPECT_THROW(ApplyCoarse(op, term, Vector(coarse.size() - 1), output), std::invalid_argument);
        EXPECT_EQ(output, sentinel);
        Matrix batch_output{sentinel};
        EXPECT_THROW(ApplyBatch(op, term, Matrix{Probe(coarse), Vector(1)}, batch_output), std::invalid_argument);
        EXPECT_EQ(batch_output, Matrix{sentinel});
    }
}

TEST(SternheimerGalerkinOperator, RejectsInvalidReferenceOrGridAndPreflightsWorkspaceBudget)
{
    const auto coarse = MakeGrid(3, 2, 2);
    const auto fine = MakeFine(MakeGrid(6, 5, 4));
    EXPECT_THROW((Galerkin(std::shared_ptr<const Hamiltonian>{}, coarse)), std::invalid_argument);
    const Hamiltonian second_order(fine->grid(), fine->local_potential(), 1.0, nullptr, 2);
    EXPECT_THROW((Galerkin(second_order, coarse)), std::invalid_argument);
    auto invalid = coarse;
    invalid.kpoint[0] += 0.1;
    EXPECT_THROW((Galerkin(fine, invalid)), std::invalid_argument);
    invalid = coarse;
    invalid.lattice_vectors[0][0] += 0.1;
    EXPECT_THROW((Galerkin(fine, invalid)), std::invalid_argument);
    invalid = coarse;
    invalid.periodic = false;
    EXPECT_THROW((Galerkin(fine, invalid)), std::invalid_argument);
    invalid = coarse;
    invalid.nx = 7;
    EXPECT_THROW((Galerkin(fine, invalid)), std::invalid_argument);
    invalid = coarse;
    invalid.nx = 0;
    EXPECT_THROW((Galerkin(fine, invalid)), std::invalid_argument);
    invalid = coarse;
    invalid.nx = std::numeric_limits<int>::max();
    EXPECT_THROW((Galerkin(fine, invalid)), std::invalid_argument);
    const auto bytes = Galerkin::workspace_bytes_required(*fine, coarse);
    EXPECT_GT(bytes, 2 * fine->grid().size() * sizeof(Complex));
    EXPECT_THROW((Galerkin(fine, coarse, bytes - 1)), std::length_error);
    Galerkin exact_budget(fine, coarse, bytes);
    EXPECT_EQ(exact_budget.workspace_bytes(), bytes);
    Vector output;
    EXPECT_NO_THROW(exact_budget.apply(Probe(coarse), output));
}

TEST(SternheimerGalerkinOperator, OccupiedFillLeavesUnshiftedVirtualBlockUnchanged)
{
    const Matrix h{{-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 3.0}};
    const auto pc = VirtualProjector(Vector{1.0, 0.0, 0.0}, 1.0);
    const auto a = ResponseOperatorWithOccupiedFill(h, pc, -1.0, 0.8);
    EXPECT_EQ(a, (Matrix{{0.8, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 4.0}}));
    const auto other_fill = ResponseOperatorWithOccupiedFill(h, pc, -1.0, 1.6);
    EXPECT_EQ(MatrixProduct(MatrixProduct(pc, a), pc), MatrixProduct(MatrixProduct(pc, other_fill), pc));
}

TEST(SternheimerGalerkinOperator, OccupiedFillCannotHideNegativeVirtualGap)
{
    const Matrix h{{-1.0, 0.0, 0.0}, {0.0, -1.2, 0.0}, {0.0, 0.0, 1.0}};
    const auto pc = VirtualProjector(Vector{1.0, 0.0, 0.0}, 1.0);
    const auto a = ResponseOperatorWithOccupiedFill(h, pc, -1.0, 0.8);
    EXPECT_NEAR(a[1][1].real(), -0.2, 1.0e-14);
    EXPECT_LT(HermitianEigenvalues(a).front(), -0.19);
}

TEST(SternheimerGalerkinOperator, PositiveAAndResponseSignDoNotValidateAnInvalidPc)
{
    const Matrix h{{0.0, 0.0}, {0.0, 2.0}};
    const Vector occupied{1.0, 0.0};
    const auto correct_pc = VirtualProjector(occupied, 1.0);
    // Using an incorrectly normalized occupied state makes Pc indefinite and
    // non-idempotent, although the filled A can still be strictly positive.
    const auto bad_pc = VirtualProjector(Vector{std::sqrt(1.5), 0.0}, 1.0);
    const auto bad_pc2 = MatrixProduct(bad_pc, bad_pc);
    EXPECT_LT(HermitianEigenvalues(bad_pc).front(), -0.49);
    EXPECT_GT(std::abs(bad_pc2[0][0] - bad_pc[0][0]), 0.74);
    const auto bad_a = ResponseOperatorWithOccupiedFill(h, bad_pc, 0.0, 0.8);
    ASSERT_GT(HermitianEigenvalues(bad_a).front(), 0.0);
    const auto bad_rhs = Multiply(bad_pc, occupied);
    EXPECT_GT(std::abs(Inner(occupied, bad_rhs)), 0.49);
    const auto bad_chi = DirectResponse(bad_a, Matrix{bad_rhs}, 1.0, 0.37);
    EXPECT_LT(bad_chi[0][0].real(), -0.1);
    EXPECT_NEAR(bad_chi[0][0].imag(), 0.0, 1.0e-14);
    const auto correct_a = ResponseOperatorWithOccupiedFill(h, correct_pc, 0.0, 0.8);
    const auto correct_chi = DirectResponse(correct_a, Matrix{Multiply(correct_pc, occupied)}, 1.0, 0.37);
    EXPECT_EQ(correct_chi, (Matrix{{0.0}}));
}

TEST(SternheimerGalerkinOperator, FineProductProjectedDirectResponseIsHermitianNonpositiveAndDielectricPositive)
{
    const auto coarse = MakeGrid(3, 2, 2);
    const auto fine_grid = MakeGrid(6, 5, 4);
    const Hamiltonian fine(fine_grid, std::vector<double>(fine_grid.size(), 1.2), 1.0, nullptr, 8);
    Galerkin op(fine, coarse);
    ModuleRI::SternheimerGridTransfer transfer(coarse, fine_grid);
    const auto source = BlochState(fine_grid);
    Vector occupied;
    Vector hsource;
    transfer.restrict_adjoint(source, occupied);
    fine.apply(source, hsource);
    const double occupied_energy = VolumeElement(fine_grid) * Inner(source, hsource).real();
    for (std::size_t i = 0; i != source.size(); ++i)
    {
        EXPECT_NEAR(std::abs(hsource[i] - occupied_energy * source[i]), 0.0, 1.0e-11);
    }
    const double dv = VolumeElement(coarse);
    EXPECT_NEAR(dv * Inner(occupied, occupied).real(), 1.0, 1.0e-12);
    const auto pc = VirtualProjector(occupied, dv);
    const auto coarse_matrix = CoarseMatrix(op);
    ExpectUnshiftedVirtualSpectrumPositive(coarse_matrix, pc, occupied, dv, occupied_energy);
    const auto a = ResponseOperatorWithOccupiedFill(coarse_matrix, pc, occupied_energy, 0.8);
    ASSERT_GT(HermitianEigenvalues(a).front(), 0.0);
    const auto products = FinePerturbationProducts(fine_grid, source);
    const auto explicit_j = ExplicitInterpolation(coarse, fine_grid);
    Matrix rhs;
    for (const auto& product: products)
    {
        Vector restricted;
        transfer.restrict_adjoint(product, restricted);
        ExpectNear(restricted, ExplicitAdjoint(explicit_j, coarse, fine_grid, product), 1.0e-12);
        rhs.push_back(Multiply(pc, restricted));
        EXPECT_NEAR(std::abs(dv * Inner(occupied, rhs.back())), 0.0, 1.0e-12);
    }
    const auto chi = DirectResponse(a, rhs, dv, 0.37);
    Matrix epsilon(chi.size(), Vector(chi.size()));
    const std::array<double, 3> coulomb{{0.4, 1.3, 2.1}};
    double off_diagonal_imaginary = 0.0;
    for (std::size_t row = 0; row != chi.size(); ++row)
    {
        EXPECT_LT(chi[row][row].real(), -1.0e-5);
        for (std::size_t col = 0; col != chi.size(); ++col)
        {
            EXPECT_NEAR(std::abs(chi[row][col] - std::conj(chi[col][row])), 0.0, 1.0e-12);
            epsilon[row][col] = (row == col ? 1.0 : 0.0) - std::sqrt(coulomb[row] * coulomb[col]) * chi[row][col];
            off_diagonal_imaginary = std::max(off_diagonal_imaginary, std::abs(chi[row][col].imag()));
        }
    }
    EXPECT_GT(off_diagonal_imaginary, 1.0e-4);
    EXPECT_LE(HermitianEigenvalues(chi).back(), 1.0e-12);
    EXPECT_GE(HermitianEigenvalues(epsilon).front(), 1.0 - 1.0e-12);
}

TEST(SternheimerGalerkinOperator, SameGridDirectResponseEqualsFullGridFineFD8DirectSolve)
{
    const auto grid = MakeGrid(3, 2, 2);
    const Hamiltonian fine(grid, std::vector<double>(grid.size(), 1.2), 1.0, nullptr, 8);
    Galerkin op(fine, grid);
    ModuleRI::SternheimerGridTransfer identity(grid, grid);
    const auto source = BlochState(grid);
    const double dv = VolumeElement(grid);
    const auto fine_matrix = fine.dense_matrix();
    const auto coarse_matrix = CoarseMatrix(op);
    ASSERT_EQ(coarse_matrix, fine_matrix);
    const double occupied_energy = dv * Inner(source, Multiply(fine_matrix, source)).real();
    const auto pc = VirtualProjector(source, dv);
    ExpectUnshiftedVirtualSpectrumPositive(fine_matrix, pc, source, dv, occupied_energy);
    ExpectUnshiftedVirtualSpectrumPositive(coarse_matrix, pc, source, dv, occupied_energy);
    const auto fine_a = ResponseOperatorWithOccupiedFill(fine_matrix, pc, occupied_energy, 0.8);
    const auto coarse_a = ResponseOperatorWithOccupiedFill(coarse_matrix, pc, occupied_energy, 0.8);
    ASSERT_GT(HermitianEigenvalues(fine_a).front(), 0.0);
    Matrix coarse_rhs;
    Matrix fine_rhs;
    for (const auto& product: FinePerturbationProducts(grid, source))
    {
        Vector restricted;
        identity.restrict_adjoint(product, restricted);
        EXPECT_EQ(restricted, product);
        coarse_rhs.push_back(Multiply(pc, restricted));
        fine_rhs.push_back(Multiply(pc, product));
    }
    const auto fine_chi = DirectResponse(fine_a, fine_rhs, dv, 0.37);
    const auto coarse_chi = DirectResponse(coarse_a, coarse_rhs, dv, 0.37);
    for (std::size_t row = 0; row != fine_chi.size(); ++row) ExpectNear(coarse_chi[row], fine_chi[row], 1.0e-13);
}

TEST(SternheimerGalerkinOperator, RestrictingPerturbationAndSourceSeparatelyLosesRetainedProduct)
{
    auto coarse = MakeGrid(3, 1, 1);
    auto fine = MakeGrid(9, 1, 1);
    ModuleRI::SternheimerGridTransfer fields(coarse, fine);
    const auto source = BlochState(fine, 1);
    Vector potential(fine.size());
    Vector product(fine.size());
    for (int i = 0; i != fine.size(); ++i)
    {
        potential[i] = std::polar(1.0, -2.0 * two_pi * i / fine.nx);
        product[i] = potential[i] * source[i];
    }
    Vector occupied;
    Vector correct;
    fields.restrict_adjoint(source, occupied);
    fields.restrict_adjoint(product, correct);
    const auto pc = VirtualProjector(occupied, VolumeElement(coarse));
    correct = Multiply(pc, correct);
    // A scalar potential has k=0, while the source carries the Bloch twist.
    // Source mode +1 and potential mode -2 produce retained mode -1. R(V) is
    // zero, so (R V)(R source) misses a unit-norm Pc R(V source) entirely.
    coarse.kpoint = {0.0, 0.0, 0.0};
    fine.kpoint = {0.0, 0.0, 0.0};
    ModuleRI::SternheimerGridTransfer scalars(coarse, fine);
    Vector restricted_potential;
    scalars.restrict_adjoint(potential, restricted_potential);
    Vector wrong(occupied.size());
    for (std::size_t i = 0; i != wrong.size(); ++i) wrong[i] = restricted_potential[i] * occupied[i];
    wrong = Multiply(pc, wrong);
    const double dv = VolumeElement(coarse);
    EXPECT_NEAR(dv * Inner(correct, correct).real(), 1.0, 1.0e-12);
    EXPECT_LT(dv * Inner(wrong, wrong).real(), 1.0e-24);
    EXPECT_GT(dv * (Inner(correct, correct) - Inner(wrong, wrong)).real(), 0.99);
}
