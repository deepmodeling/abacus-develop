#include "source_lcao/module_ri/sternheimer_galerkin_audit.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace
{

using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
using Grid = Hamiltonian::Grid;
using Complex = Hamiltonian::Complex;
using Vector = Hamiltonian::Vector;
using Matrix = Hamiltonian::Matrix;
using Projector = ModuleRI::SternheimerFDNonlocalProjector;
using Result = ModuleRI::SternheimerGalerkinAuditMatrices;
constexpr std::size_t test_budget = 64ULL * 1024 * 1024;

std::vector<const Vector*> views(const Matrix& states)
{
    std::vector<const Vector*> result;
    for (const auto& state : states)
    {
        result.push_back(&state);
    }
    return result;
}

Result audit(const std::shared_ptr<const Hamiltonian>& fine, const Matrix& states, const Grid& coarse,
             const std::size_t budget = test_budget)
{
    return ModuleRI::audit_sternheimer_galerkin_matrices(fine, views(states), coarse, budget);
}

Grid make_grid(const int nx, const int ny, const int nz)
{
    Grid grid{nx, ny, nz, 4.0 / nx, std::sqrt(12.74) / ny, std::sqrt(9.2) / nz, true};
    grid.lattice_vectors = {{{4.0, 0.0, 0.0}, {0.7, 3.5, 0.0}, {0.2, 0.4, 3.0}}};
    grid.kpoint = {0.21, -0.17, 0.13};
    return grid;
}

double dvolume(const Grid& grid)
{
    return 42.0 / grid.size();
}

std::array<double, 3> reduced_position(const Grid& grid, const int index)
{
    return {static_cast<double>(index / (grid.ny * grid.nz)) / grid.nx,
            static_cast<double>((index / grid.nz) % grid.ny) / grid.ny,
            static_cast<double>(index % grid.nz) / grid.nz};
}

Vector mode(const Grid& grid, const std::array<int, 3>& frequencies)
{
    Vector values(grid.size());
    for (int i = 0; i != grid.size(); ++i)
    {
        const auto r = reduced_position(grid, i);
        double phase = 0.0;
        for (int axis = 0; axis != 3; ++axis)
        {
            phase += (frequencies[axis] + grid.kpoint[axis]) * r[axis];
        }
        values[i] = std::polar(1.0, 2.0 * std::acos(-1.0) * phase);
    }
    return values;
}

Vector probe(const Grid& grid, const double seed)
{
    Vector values(grid.size());
    for (int i = 0; i != grid.size(); ++i)
    {
        values[i] = Complex(std::sin(0.31 * i + seed), std::cos(0.17 * i - 2.0 * seed));
    }
    return values;
}

Complex dot(const Vector& left, const Vector& right, const Grid& grid)
{
    Complex result(0.0, 0.0);
    for (std::size_t i = 0; i != left.size(); ++i)
    {
        result += std::conj(left[i]) * right[i];
    }
    return dvolume(grid) * result;
}

std::shared_ptr<const Hamiltonian> make_fine(const Grid& grid, const int order = 8)
{
    std::vector<double> potential(grid.size());
    Projector::ProjectorBlock block;
    block.projectors = {mode(grid, {0, 0, 0}), mode(grid, {1, 0, 0})};
    block.d_matrix = {{Complex(0.08, 0.0), Complex(0.015, 0.02)},
                      {Complex(0.015, -0.02), Complex(-0.04, 0.0)}};
    for (int i = 0; i != grid.size(); ++i)
    {
        const auto r = reduced_position(grid, i);
        potential[i] = -0.7 + 0.23 * std::cos(2.0 * std::acos(-1.0) * (r[0] - r[1]));
    }
    auto nonlocal = std::make_shared<const Projector>(grid.size(), dvolume(grid),
                                                      std::vector<Projector::ProjectorBlock>{block});
    return std::make_shared<const Hamiltonian>(grid, potential, 1.0, nonlocal, order);
}

// Construct J by direct Fourier sums, independently of the FFT transfer helper.
Matrix interpolation(const Grid& coarse, const Grid& fine)
{
    Matrix j(fine.size(), Vector(coarse.size(), Complex(0.0, 0.0)));
    for (int mx = -((coarse.nx - 1) / 2); mx <= coarse.nx / 2; ++mx)
    {
        for (int my = -((coarse.ny - 1) / 2); my <= coarse.ny / 2; ++my)
        {
            for (int mz = -((coarse.nz - 1) / 2); mz <= coarse.nz / 2; ++mz)
            {
                const Vector c = mode(coarse, {mx, my, mz});
                const Vector f = mode(fine, {mx, my, mz});
                for (int row = 0; row != fine.size(); ++row)
                {
                    for (int col = 0; col != coarse.size(); ++col)
                    {
                        j[row][col] += f[row] * std::conj(c[col]) / static_cast<double>(coarse.size());
                    }
                }
            }
        }
    }
    return j;
}

Vector lift(const Matrix& j, const Vector& coarse)
{
    Vector fine(j.size(), Complex(0.0, 0.0));
    for (std::size_t row = 0; row != j.size(); ++row)
    {
        for (std::size_t col = 0; col != coarse.size(); ++col)
        {
            fine[row] += j[row][col] * coarse[col];
        }
    }
    return fine;
}

Vector restrict_direct(const Matrix& j, const Vector& fine, const Grid& coarse_grid, const Grid& fine_grid)
{
    Vector coarse(coarse_grid.size(), Complex(0.0, 0.0));
    for (std::size_t col = 0; col != coarse.size(); ++col)
    {
        for (std::size_t row = 0; row != fine.size(); ++row)
        {
            coarse[col] += std::conj(j[row][col]) * fine[row] * (dvolume(fine_grid) / dvolume(coarse_grid));
        }
    }
    return coarse;
}

void check_matrices(const Result& result, const Hamiltonian& fine, const Matrix& lifted)
{
    const std::size_t n = lifted.size();
    ASSERT_EQ(result.dimension, n);
    ASSERT_EQ(result.overlap.size(), n * n);
    EXPECT_NEAR(result.filtered_matrix_identity_max_abs_error, 0.0, 3.0e-9);
    using Apply = void (Hamiltonian::*)(const Vector&, Vector&) const;
    const std::array<Apply, 4> apply{{&Hamiltonian::apply_kinetic, &Hamiltonian::apply_local_potential,
                                     &Hamiltonian::apply_nonlocal, &Hamiltonian::apply}};
    const std::array<const Vector*, 4> actual{{&result.kinetic, &result.local_potential,
                                              &result.nonlocal, &result.hamiltonian}};
    for (std::size_t term = 0; term != apply.size(); ++term)
    {
        ASSERT_EQ(actual[term]->size(), n * n);
        for (std::size_t col = 0; col != n; ++col)
        {
            Vector image;
            (fine.*apply[term])(lifted[col], image);
            for (std::size_t row = 0; row != n; ++row)
            {
                EXPECT_NEAR(std::abs((*actual[term])[row + n * col] - dot(lifted[row], image, fine.grid())),
                            0.0, 3.0e-9) << "term " << term << " row " << row << " col " << col;
            }
        }
    }
    for (std::size_t col = 0; col != n; ++col)
    {
        for (std::size_t row = 0; row != n; ++row)
        {
            const std::size_t i = row + n * col;
            EXPECT_NEAR(std::abs(result.overlap[i] - dot(lifted[row], lifted[col], fine.grid())), 0.0, 2.0e-10);
            EXPECT_NEAR(std::abs(result.hamiltonian[i] - result.kinetic[i]
                                - result.local_potential[i] - result.nonlocal[i]), 0.0, 3.0e-9);
            EXPECT_NEAR(std::abs(result.hamiltonian[i] - std::conj(result.hamiltonian[col + n * row])),
                        0.0, 3.0e-9);
        }
    }
}

} // namespace

TEST(SternheimerGalerkinAudit, MatchesExplicitBlochTransferAndFineComponentForms)
{
    const Grid coarse = make_grid(2, 3, 2);
    const auto fine = make_fine(make_grid(5, 4, 3));
    const Matrix states{probe(fine->grid(), 0.1), probe(fine->grid(), 0.8), probe(fine->grid(), -0.4)};
    const Matrix original = states;
    const auto result = audit(fine, states, coarse);
    const Matrix j = interpolation(coarse, fine->grid());
    Matrix lifted;
    for (std::size_t i = 0; i != states.size(); ++i)
    {
        const Vector projected = restrict_direct(j, states[i], coarse, fine->grid());
        EXPECT_NEAR(result.fine_norm_squared[i], dot(states[i], states[i], fine->grid()).real(), 2.0e-11);
        EXPECT_NEAR(result.projected_norm_squared[i], dot(projected, projected, coarse).real(), 2.0e-11);
        lifted.push_back(lift(j, projected));
        Vector difference = lifted.back();
        for (std::size_t point = 0; point != difference.size(); ++point)
        {
            difference[point] -= states[i][point];
        }
        EXPECT_NEAR(result.reconstruction_relative_error[i],
                    std::sqrt(dot(difference, difference, fine->grid()).real() / result.fine_norm_squared[i]),
                    2.0e-11);
        EXPECT_NEAR(result.fine_norm_squared[i] - result.projected_norm_squared[i],
                    dot(difference, difference, fine->grid()).real(), 2.0e-10);
    }
    check_matrices(result, *fine, lifted);
    EXPECT_EQ(states, original);
}

TEST(SternheimerGalerkinAudit, SameGridRetainsOriginalUnnormalizedStates)
{
    const auto fine = make_fine(make_grid(4, 3, 3));
    const Matrix states{probe(fine->grid(), 0.1), probe(fine->grid(), 1.1)};
    const auto result = audit(fine, states, fine->grid());
    EXPECT_EQ(result.fine_norm_squared, result.projected_norm_squared);
    EXPECT_EQ(result.reconstruction_relative_error, std::vector<double>(states.size(), 0.0));
    check_matrices(result, *fine, states);
}

TEST(SternheimerGalerkinAudit, KeepsZeroDependentScaledStatesInOriginalOrder)
{
    const auto fine = make_fine(make_grid(5, 4, 3));
    const Grid coarse = make_grid(3, 3, 2);
    const Vector state = mode(fine->grid(), {1, 0, 0});
    Vector scaled = state;
    const Complex scale(2.0, -0.7);
    for (auto& value : scaled)
    {
        value *= scale;
    }
    const Matrix states{Vector(state.size(), Complex(0.0, 0.0)), scaled, state, scaled};
    const auto result = audit(fine, states, coarse);
    ASSERT_EQ(result.dimension, 4U);
    EXPECT_DOUBLE_EQ(result.projected_norm_squared[0], 0.0);
    EXPECT_DOUBLE_EQ(result.reconstruction_relative_error[0], 0.0);
    EXPECT_NEAR(result.projected_norm_squared[1], std::norm(scale) * result.projected_norm_squared[2], 2.0e-10);
    EXPECT_NEAR(result.projected_norm_squared[1], result.projected_norm_squared[3], 2.0e-11);
    for (std::size_t i = 0; i != 4; ++i)
    {
        EXPECT_NEAR(std::abs(result.overlap[i * 4]), 0.0, 2.0e-11);
        EXPECT_NEAR(std::abs(result.hamiltonian[i]), 0.0, 2.0e-11);
    }
    EXPECT_NEAR(std::abs(result.overlap[1 + 4 * 2] - std::conj(scale) * result.overlap[2 + 4 * 2]),
                0.0, 2.0e-10);
}

TEST(SternheimerGalerkinAudit, FineOnlyModeRemainsAZeroColumnRatherThanBeingTruncated)
{
    const auto fine = make_fine(make_grid(7, 5, 5));
    const Matrix states{mode(fine->grid(), {1, 0, 0}), mode(fine->grid(), {2, 0, 0})};
    const auto result = audit(fine, states, make_grid(3, 3, 3));
    ASSERT_EQ(result.dimension, 2U);
    EXPECT_GT(result.fine_norm_squared[1], 1.0);
    EXPECT_NEAR(result.projected_norm_squared[1], 0.0, 2.0e-11);
    EXPECT_NEAR(result.reconstruction_relative_error[1], 1.0, 2.0e-11);
    EXPECT_NEAR(std::abs(result.overlap[3]), 0.0, 2.0e-11);
    const auto unfiltered = audit(fine, states, fine->grid());
    EXPECT_GT(unfiltered.overlap[3].real(), 1.0);
    EXPECT_GT(std::abs(unfiltered.hamiltonian[3] - result.hamiltonian[3]), 1.0);
    EXPECT_NEAR(result.filtered_matrix_identity_max_abs_error, 0.0, 3.0e-9);
}

TEST(SternheimerGalerkinAudit, MultipleStateBlocksAndTailMatchDirectForms)
{
    const auto fine = make_fine(make_grid(5, 4, 3));
    const Grid coarse = make_grid(3, 3, 2);
    Matrix states;
    Matrix lifted;
    const Matrix j = interpolation(coarse, fine->grid());
    for (int state = 0; state != 17; ++state)
    {
        states.push_back(probe(fine->grid(), 0.13 * state));
        lifted.push_back(lift(j, restrict_direct(j, states.back(), coarse, fine->grid())));
    }
    check_matrices(audit(fine, states, coarse), *fine, lifted);
}

TEST(SternheimerGalerkinAudit, GridContractionHandlesPartialFinalTile)
{
    const Grid grid = make_grid(21, 21, 19);
    auto fine = std::make_shared<const Hamiltonian>(grid, std::vector<double>(grid.size(), -0.7), 1.0, nullptr, 8);
    const Matrix states{mode(grid, {1, 0, 0}), mode(grid, {0, 1, 1})};
    const auto result = audit(fine, states, grid);
    check_matrices(result, *fine, states);
}

TEST(SternheimerGalerkinAudit, RejectsInsufficientBudgetAndAcceptsExactEstimate)
{
    const auto fine = make_fine(make_grid(5, 4, 3));
    const Grid coarse = make_grid(3, 3, 2);
    const Matrix states{probe(fine->grid(), 0.1), probe(fine->grid(), 0.8)};
    const auto required = ModuleRI::sternheimer_galerkin_audit_workspace_bytes(*fine, coarse, states.size());
    ASSERT_GT(required, 0U);
    EXPECT_THROW(audit(fine, states, coarse, required - 1), std::length_error);
    const auto result = audit(fine, states, coarse, required);
    EXPECT_EQ(result.workspace_bytes, required);
}

TEST(SternheimerGalerkinAudit, BudgetCountsFilteredCacheButIdentityBorrowsFullFineStates)
{
    const auto fine = make_fine(make_grid(31, 31, 31));
    const Grid coarse = make_grid(3, 3, 3);
    const std::size_t count = 176;
    const std::size_t fine_cache_bytes = sizeof(Complex) * fine->grid().size() * count;
    const auto projected = ModuleRI::sternheimer_galerkin_audit_workspace_bytes(*fine, coarse, count);
    const auto identity = ModuleRI::sternheimer_galerkin_audit_workspace_bytes(*fine, fine->grid(), count);
    EXPECT_GT(projected, fine_cache_bytes);
    EXPECT_LT(identity, fine_cache_bytes);
    // Repeated borrowed views represent a valid dependent basis, not permission
    // to deduplicate columns or bypass the full-cache budget check.
    const Vector state = mode(fine->grid(), {0, 0, 0});
    const std::vector<const Vector*> borrowed(count, &state);
    EXPECT_THROW(ModuleRI::audit_sternheimer_galerkin_matrices(fine, borrowed, coarse, projected - 1),
                 std::length_error);
}

TEST(SternheimerGalerkinAudit, RejectsMissingOrMalformedStates)
{
    const auto fine = make_fine(make_grid(5, 4, 3));
    const Grid coarse = make_grid(3, 3, 2);
    const Matrix valid{probe(fine->grid(), 0.1)};
    EXPECT_THROW(audit(nullptr, valid, coarse), std::invalid_argument);
    EXPECT_THROW(audit(fine, Matrix{}, coarse), std::invalid_argument);
    EXPECT_THROW(ModuleRI::audit_sternheimer_galerkin_matrices(fine, {nullptr}, coarse, test_budget),
                 std::invalid_argument);
    Matrix malformed = valid;
    malformed.push_back(Vector(fine->grid().size() - 1));
    EXPECT_THROW(audit(fine, malformed, coarse), std::invalid_argument);
    malformed = valid;
    malformed[0][0] = Complex(std::numeric_limits<double>::quiet_NaN(), 0.0);
    EXPECT_THROW(audit(fine, malformed, coarse), std::invalid_argument);
    malformed[0][0] = Complex(0.0, std::numeric_limits<double>::infinity());
    EXPECT_THROW(audit(fine, malformed, coarse), std::invalid_argument);
}

TEST(SternheimerGalerkinAudit, RejectsInvalidGridAndReference)
{
    const auto fine = make_fine(make_grid(5, 4, 3));
    const Grid coarse = make_grid(3, 3, 2);
    const Matrix states{probe(fine->grid(), 0.1)};
    Grid mismatched = coarse;
    mismatched.kpoint[0] += 0.01;
    EXPECT_THROW(audit(fine, states, mismatched), std::invalid_argument);
    mismatched = coarse;
    mismatched.lattice_vectors[1][0] += 0.1;
    EXPECT_THROW(audit(fine, states, mismatched), std::invalid_argument);
    mismatched = coarse;
    mismatched.periodic = false;
    EXPECT_THROW(audit(fine, states, mismatched), std::invalid_argument);
    mismatched = coarse;
    mismatched.nx = 0;
    EXPECT_THROW(audit(fine, states, mismatched), std::invalid_argument);
    EXPECT_THROW(audit(fine, states, make_grid(6, 3, 2)), std::invalid_argument);
    EXPECT_THROW(audit(make_fine(fine->grid(), 2), states, coarse), std::invalid_argument);
}

TEST(SternheimerGalerkinAudit, EstimatorRejectsInvalidCountsAndOverflowWithoutAllocation)
{
    const auto fine = make_fine(make_grid(5, 4, 3));
    const Grid coarse = make_grid(3, 3, 2);
    EXPECT_THROW(ModuleRI::sternheimer_galerkin_audit_workspace_bytes(*fine, coarse, 0), std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_galerkin_audit_workspace_bytes(
                     *fine, coarse, static_cast<std::size_t>(std::numeric_limits<int>::max()) + 1),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_galerkin_audit_workspace_bytes(
                     *fine, coarse, static_cast<std::size_t>(std::numeric_limits<int>::max())),
                 std::length_error);
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
