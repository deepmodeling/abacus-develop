#include "source_lcao/module_ri/sternheimer_grid_transfer.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{

using Transfer = ModuleRI::SternheimerGridTransfer;
using Grid = Transfer::Grid;
using Vector = Transfer::Vector;
using Complex = ModuleRI::SternheimerFDHamiltonian::Complex;
using Lattice = ModuleRI::SternheimerFDLatticeVectors;

const Lattice oblique_cell{{{5.0, 0.0, 0.0}, {-2.5, 4.3, 0.0}, {0.4, -0.2, 8.0}}};

std::size_t size(const Grid& grid)
{
    return static_cast<std::size_t>(grid.nx) * grid.ny * grid.nz;
}

Grid make_grid(const std::array<int, 3>& dimensions,
               const ModuleRI::SternheimerReducedKPoint& kpoint = {0.0, 0.0, 0.0})
{
    Grid grid;
    grid.nx = dimensions[0];
    grid.ny = dimensions[1];
    grid.nz = dimensions[2];
    grid.lattice_vectors = oblique_cell;
    grid.kpoint = kpoint;
    const auto length = [](const std::array<double, 3>& a) {
        return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    };
    grid.hx = length(oblique_cell[0]) / grid.nx;
    grid.hy = length(oblique_cell[1]) / grid.ny;
    grid.hz = length(oblique_cell[2]) / grid.nz;
    return grid;
}

double volume(const Grid& grid)
{
    const auto a = ModuleRI::sternheimer_fd_grid_lattice_vectors(grid);
    return std::abs(a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                    - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                    + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]));
}

Complex dot(const Vector& left, const Vector& right, const Grid& grid)
{
    Complex result(0.0, 0.0);
    for (std::size_t i = 0; i != left.size(); ++i)
    {
        result += std::conj(left[i]) * right[i];
    }
    return result * (volume(grid) / static_cast<double>(size(grid)));
}

Vector samples(const Grid& grid, const double seed)
{
    Vector result(size(grid));
    for (std::size_t i = 0; i != result.size(); ++i)
    {
        const double x = static_cast<double>(i);
        result[i] = Complex(std::sin(0.37 * x + seed), std::cos(0.23 * x - 2.0 * seed));
    }
    return result;
}

Vector mode(const Grid& grid, const std::array<int, 3>& frequencies, const double phase = 0.0)
{
    Vector result(size(grid));
    const double two_pi = 2.0 * std::acos(-1.0);
    for (int ix = 0; ix != grid.nx; ++ix)
    {
        for (int iy = 0; iy != grid.ny; ++iy)
        {
            for (int iz = 0; iz != grid.nz; ++iz)
            {
                const double angle = phase + two_pi
                    * ((frequencies[0] + grid.kpoint[0]) * ix / grid.nx
                       + (frequencies[1] + grid.kpoint[1]) * iy / grid.ny
                       + (frequencies[2] + grid.kpoint[2]) * iz / grid.nz);
                const std::size_t index = (static_cast<std::size_t>(ix) * grid.ny + iy) * grid.nz + iz;
                result[index] = std::exp(Complex(0.0, angle));
            }
        }
    }
    return result;
}

void expect_near(const Vector& actual, const Vector& expected, const double tolerance = 2.0e-12)
{
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i != actual.size(); ++i)
    {
        EXPECT_NEAR(std::abs(actual[i] - expected[i]), 0.0, tolerance) << "index " << i;
    }
}

// Independently solve a_i dot q = 2 pi (m_i + k_i) for the triangular test cell.
std::array<double, 3> wavevector(const Grid& grid, const std::array<int, 3>& frequencies)
{
    const double two_pi = 2.0 * std::acos(-1.0);
    const double qx = two_pi * (frequencies[0] + grid.kpoint[0]) / 5.0;
    const double qy = (two_pi * (frequencies[1] + grid.kpoint[1]) + 2.5 * qx) / 4.3;
    const double qz = (two_pi * (frequencies[2] + grid.kpoint[2]) - 0.4 * qx + 0.2 * qy) / 8.0;
    return {qx, qy, qz};
}

std::array<Vector, 3> direct_gradients(const Grid& coarse, const Grid& fine, const Vector& input)
{
    std::array<Vector, 3> expected;
    for (auto& field : expected)
    {
        field.resize(size(fine));
    }
    for (int mx = -(coarse.nx - 1) / 2; mx <= coarse.nx / 2; ++mx)
    {
        for (int my = -(coarse.ny - 1) / 2; my <= coarse.ny / 2; ++my)
        {
            for (int mz = -(coarse.nz - 1) / 2; mz <= coarse.nz / 2; ++mz)
            {
                const std::array<int, 3> frequencies{mx, my, mz};
                const Vector fine_mode = mode(fine, frequencies);
                const Complex coefficient = dot(mode(coarse, frequencies), input, coarse) / volume(coarse);
                const auto q = wavevector(coarse, frequencies);
                for (int axis = 0; axis != 3; ++axis)
                {
                    for (std::size_t i = 0; i != size(fine); ++i)
                    {
                        expected[axis][i] += Complex(0.0, q[axis]) * coefficient * fine_mode[i];
                    }
                }
            }
        }
    }
    return expected;
}

} // namespace

TEST(SternheimerGridTransfer, SameGridIsExactIdentityIncludingInPlace)
{
    const Grid grid = make_grid({4, 5, 3}, {0.21, -0.17, 0.08});
    Transfer transfer(grid, grid);
    const Vector input = samples(grid, 0.31);
    Vector output;
    transfer.interpolate(input, output);
    EXPECT_EQ(output, input);
    transfer.restrict_adjoint(input, output);
    EXPECT_EQ(output, input);
    transfer.interpolate(output, output);
    transfer.restrict_adjoint(output, output);
    EXPECT_EQ(output, input);
    EXPECT_EQ(transfer.coarse_grid().nx, grid.nx);
    EXPECT_EQ(transfer.fine_grid().kpoint, grid.kpoint);
}

TEST(SternheimerGridTransfer, PreservesWeightedNormOnObliqueBlochGrid)
{
    const Grid coarse = make_grid({4, 5, 3}, {0.21, -0.17, 0.08});
    const Grid fine = make_grid({9, 8, 7}, coarse.kpoint);
    Transfer transfer(coarse, fine);
    const Vector input = samples(coarse, 0.31);
    Vector output;
    transfer.interpolate(input, output);
    EXPECT_NEAR(dot(output, output, fine).real(), dot(input, input, coarse).real(), 2.0e-10);
}

TEST(SternheimerGridTransfer, RestrictionIsVolumeWeightedAdjoint)
{
    const Grid coarse = make_grid({5, 4, 3}, {-0.31, 0.13, 0.27});
    const Grid fine = make_grid({8, 9, 7}, coarse.kpoint);
    Transfer transfer(coarse, fine);
    const Vector x = samples(coarse, 0.19);
    const Vector y = samples(fine, 0.73);
    Vector jx;
    Vector adjoint_y;
    transfer.interpolate(x, jx);
    transfer.restrict_adjoint(y, adjoint_y);
    EXPECT_NEAR(std::abs(dot(jx, y, fine) - dot(x, adjoint_y, coarse)), 0.0, 2.0e-11);
}

TEST(SternheimerGridTransfer, RoundTripRecoversArbitraryComplexCoarseVector)
{
    const Grid coarse = make_grid({4, 3, 5}, {0.17, -0.11, 0.08});
    const Grid fine = make_grid({7, 8, 9}, coarse.kpoint);
    Transfer transfer(coarse, fine);
    const Vector input = samples(coarse, 0.41);
    Vector lifted;
    Vector recovered;
    transfer.interpolate(input, lifted);
    transfer.restrict_adjoint(lifted, recovered);
    expect_near(recovered, input);
}

TEST(SternheimerGridTransfer, MatchesIndependentDirectFourierSum)
{
    const Grid coarse = make_grid({3, 4, 2}, {0.17, -0.11, 0.08});
    const Grid fine = make_grid({5, 7, 4}, coarse.kpoint);
    Transfer transfer(coarse, fine);
    const Vector input = samples(coarse, 0.37);
    Vector expected(size(fine), Complex(0.0, 0.0));
    for (int mx = -(coarse.nx - 1) / 2; mx <= coarse.nx / 2; ++mx)
    {
        for (int my = -(coarse.ny - 1) / 2; my <= coarse.ny / 2; ++my)
        {
            for (int mz = -(coarse.nz - 1) / 2; mz <= coarse.nz / 2; ++mz)
            {
                const Vector coarse_mode = mode(coarse, {mx, my, mz});
                const Vector fine_mode = mode(fine, {mx, my, mz});
                const Complex coefficient = dot(coarse_mode, input, coarse) / volume(coarse);
                for (std::size_t i = 0; i != expected.size(); ++i)
                {
                    expected[i] += coefficient * fine_mode[i];
                }
            }
        }
    }
    Vector actual;
    transfer.interpolate(input, actual);
    expect_near(actual, expected);
}

TEST(SternheimerGridTransfer, InterpolatesAnisotropicObliqueBlochMode)
{
    const Grid coarse = make_grid({5, 4, 3}, {0.17, -0.11, 0.08});
    const Grid fine = make_grid({8, 7, 3}, coarse.kpoint);
    Transfer transfer(coarse, fine);
    Vector actual;
    transfer.interpolate(mode(coarse, {-2, 1, -1}, 0.47), actual);
    expect_near(actual, mode(fine, {-2, 1, -1}, 0.47));
}

TEST(SternheimerGridTransfer, RemovesFineOnlyModesWithoutAliasing)
{
    const Grid coarse = make_grid({4, 3, 3}, {0.17, -0.11, 0.08});
    const Grid fine = make_grid({9, 7, 7}, coarse.kpoint);
    Transfer transfer(coarse, fine);
    for (const auto& frequencies : std::vector<std::array<int, 3>>{{3, 0, 0}, {0, -2, 0}, {0, 0, 2}})
    {
        Vector actual;
        transfer.restrict_adjoint(mode(fine, frequencies, 0.47), actual);
        expect_near(actual, Vector(size(coarse), Complex(0.0, 0.0)));
    }
}

TEST(SternheimerGridTransfer, PositiveNyquistRetainsComplexPhaseWithoutPairMerging)
{
    const Grid coarse = make_grid({4, 6, 4}, {0.21, 0.08, -0.17});
    const Grid fine = make_grid({10, 12, 9}, coarse.kpoint);
    Transfer transfer(coarse, fine);
    const Vector input = mode(coarse, {2, 3, 2}, 0.63);
    Vector actual;
    transfer.interpolate(input, actual);
    expect_near(actual, mode(fine, {2, 3, 2}, 0.63));
    transfer.restrict_adjoint(actual, actual);
    expect_near(actual, input);
    transfer.restrict_adjoint(mode(fine, {-2, 0, 0}, 0.63), actual);
    expect_near(actual, Vector(size(coarse), Complex(0.0, 0.0)));
}

TEST(SternheimerGridTransfer, IntegerShiftedBlochLabelsPreserveCommonPhysicalMode)
{
    const Grid coarse = make_grid({6, 6, 4}, {0.17, -0.11, 0.08});
    const Grid fine = make_grid({11, 10, 9}, coarse.kpoint);
    Grid shifted_coarse = coarse;
    shifted_coarse.kpoint[0] += 1.0;
    shifted_coarse.kpoint[1] -= 1.0;
    Grid shifted_fine = fine;
    shifted_fine.kpoint = shifted_coarse.kpoint;
    Transfer original(coarse, fine);
    Transfer shifted(shifted_coarse, shifted_fine);
    Vector original_result;
    Vector shifted_result;
    original.interpolate(mode(coarse, {1, -1, 1}, 0.29), original_result);
    shifted.interpolate(mode(shifted_coarse, {0, 0, 1}, 0.29), shifted_result);
    expect_near(shifted_result, original_result);
    expect_near(shifted_result, mode(fine, {1, -1, 1}, 0.29));
}

TEST(SternheimerGridTransfer, SupportsInPlaceTransferAndIndependentInstances)
{
    const Grid coarse = make_grid({3, 4, 5}, {0.1, 0.2, -0.3});
    const Grid fine = make_grid({7, 8, 9}, coarse.kpoint);
    Transfer first(coarse, fine);
    Transfer second(coarse, fine);
    const Vector input = samples(coarse, 0.31);
    Vector in_place = input;
    Vector other;
    first.interpolate(in_place, in_place);
    second.interpolate(samples(coarse, 0.79), other);
    first.restrict_adjoint(in_place, in_place);
    expect_near(in_place, input);
    second.restrict_adjoint(other, other);
    expect_near(other, samples(coarse, 0.79));
}

TEST(SternheimerGridTransfer, SeparateWorkersOwnIndependentPlans)
{
    std::vector<Vector> results(4);
#ifdef _OPENMP
#pragma omp parallel for num_threads(4)
#endif
    for (int worker = 0; worker < 4; ++worker)
    {
        const Grid coarse = make_grid({3 + worker, 4, 3}, {0.13 * worker, -0.11, 0.08});
        const Grid fine = make_grid({8 + worker, 7, 6}, coarse.kpoint);
        Transfer transfer(coarse, fine);
        Vector lifted;
        transfer.interpolate(samples(coarse, 0.31), lifted);
        transfer.restrict_adjoint(lifted, results[static_cast<std::size_t>(worker)]);
    }
    for (int worker = 0; worker != 4; ++worker)
    {
        const Grid coarse = make_grid({3 + worker, 4, 3}, {0.13 * worker, -0.11, 0.08});
        expect_near(results[static_cast<std::size_t>(worker)], samples(coarse, 0.31));
    }
}

TEST(SternheimerGridTransfer, ReportsNumericWorkspaceWithoutAllocatingIt)
{
    const Grid coarse = make_grid({3, 4, 5});
    const Grid fine = make_grid({7, 8, 9});
    const std::size_t expected = (size(coarse) + size(fine)) * (2 * sizeof(double) + sizeof(Complex))
                                 + size(coarse) * sizeof(std::size_t);
    EXPECT_EQ(Transfer::workspace_bytes_required(coarse, fine), expected);
    EXPECT_EQ(Transfer::workspace_bytes_required(coarse, coarse), 0U);
    Grid invalid = fine;
    invalid.nx = 2;
    EXPECT_THROW(Transfer::workspace_bytes_required(coarse, invalid), std::invalid_argument);
}

TEST(SternheimerGridTransfer, AcceptsImplicitAndExplicitIdenticalOrthogonalCells)
{
    Grid coarse{4, 3, 2, 1.0, 2.0, 3.0, true};
    Grid fine{8, 9, 4, 0.5, 2.0 / 3.0, 1.5, true};
    fine.lattice_vectors = {{{4.0, 0.0, 0.0}, {0.0, 6.0, 0.0}, {0.0, 0.0, 6.0}}};
    Transfer transfer(coarse, fine);
    Vector actual;
    transfer.interpolate(mode(coarse, {1, 1, 0}), actual);
    expect_near(actual, mode(fine, {1, 1, 0}));
}

TEST(SternheimerGridTransfer, RejectsMismatchedCellKPointAndNonperiodicGrids)
{
    const Grid coarse = make_grid({4, 4, 4}, {0.17, -0.11, 0.08});
    const Grid fine = make_grid({8, 8, 8}, coarse.kpoint);
    Grid invalid = fine;
    invalid.lattice_vectors[1][0] += 0.01;
    EXPECT_THROW((Transfer(coarse, invalid)), std::invalid_argument);
    invalid = fine;
    invalid.kpoint[0] += 1.0;
    EXPECT_THROW((Transfer(coarse, invalid)), std::invalid_argument);
    invalid = fine;
    invalid.periodic = false;
    EXPECT_THROW((Transfer(coarse, invalid)), std::invalid_argument);
    invalid = coarse;
    invalid.periodic = false;
    EXPECT_THROW((Transfer(invalid, fine)), std::invalid_argument);
}

TEST(SternheimerGridTransfer, RejectsInvalidDimensionsAndNonfiniteGeometry)
{
    const Grid coarse = make_grid({4, 4, 4});
    const Grid fine = make_grid({8, 8, 8});
    Grid invalid = fine;
    invalid.nx = 3;
    EXPECT_THROW((Transfer(coarse, invalid)), std::invalid_argument);
    invalid = coarse;
    invalid.ny = 0;
    EXPECT_THROW((Transfer(invalid, fine)), std::invalid_argument);
    invalid.ny = -1;
    EXPECT_THROW((Transfer(invalid, fine)), std::invalid_argument);
    invalid = fine;
    invalid.nx = std::numeric_limits<int>::max();
    EXPECT_THROW((Transfer(coarse, invalid)), std::invalid_argument);
    invalid = fine;
    invalid.hz = 0.0;
    EXPECT_THROW((Transfer(coarse, invalid)), std::invalid_argument);
    invalid = fine;
    invalid.kpoint[1] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW((Transfer(coarse, invalid)), std::invalid_argument);
    invalid = fine;
    invalid.lattice_vectors[0][0] = std::numeric_limits<double>::infinity();
    EXPECT_THROW((Transfer(coarse, invalid)), std::invalid_argument);
    invalid = fine;
    invalid.lattice_vectors[2] = invalid.lattice_vectors[0];
    EXPECT_THROW((Transfer(coarse, invalid)), std::invalid_argument);
}

TEST(SternheimerGridTransfer, RejectsWrongVectorSizesBeforeIdentityOrFFT)
{
    const Grid coarse = make_grid({3, 4, 5});
    const Grid fine = make_grid({7, 8, 9});
    Transfer transfer(coarse, fine);
    Transfer identity(coarse, coarse);
    Vector output;
    EXPECT_THROW(transfer.interpolate(Vector(size(coarse) - 1), output), std::invalid_argument);
    EXPECT_THROW(transfer.restrict_adjoint(Vector(size(fine) + 1), output), std::invalid_argument);
    EXPECT_THROW(identity.interpolate(Vector(), output), std::invalid_argument);
    EXPECT_THROW(identity.restrict_adjoint(Vector(), output), std::invalid_argument);
}

TEST(SternheimerGridTransfer, MoveTransfersOwnedPlans)
{
    const Grid coarse = make_grid({3, 4, 5});
    const Grid fine = make_grid({7, 8, 9});
    Transfer original(coarse, fine);
    Transfer moved(std::move(original));
    Transfer assigned(coarse, coarse);
    assigned = std::move(moved);
    Vector actual;
    assigned.interpolate(mode(coarse, {1, -1, 2}), actual);
    expect_near(actual, mode(fine, {1, -1, 2}));
}

TEST(SternheimerGridTransfer, GradientsMatchIndependentObliqueBlochFourierSum)
{
    const Grid coarse = make_grid({3, 4, 2}, {0.17, -0.11, 0.08});
    const Grid fine = make_grid({5, 7, 4}, coarse.kpoint);
    const Vector input = samples(coarse, 0.37);
    Transfer transfer(coarse, fine);
    Vector values;
    Vector original_values;
    std::array<Vector, 3> gradients;
    transfer.interpolate(input, original_values);
    transfer.interpolate_with_gradients(input, values, gradients);
    expect_near(values, original_values);
    const auto expected = direct_gradients(coarse, fine, input);
    for (int axis = 0; axis != 3; ++axis)
    {
        expect_near(gradients[axis], expected[axis]);
    }
}

TEST(SternheimerGridTransfer, GradientRestrictionIsWeightedAdjointAndMatchesDirectSum)
{
    const Grid coarse = make_grid({3, 4, 2}, {-0.31, 0.13, 0.27});
    const Grid fine = make_grid({5, 7, 4}, coarse.kpoint);
    const Vector input = samples(coarse, 0.19);
    const std::array<Vector, 3> fields{{samples(fine, 0.73), samples(fine, -0.41), samples(fine, 0.17)}};
    Transfer transfer(coarse, fine);
    Vector values;
    Vector adjoint;
    std::array<Vector, 3> gradients;
    transfer.interpolate_with_gradients(input, values, gradients);
    transfer.restrict_gradient_adjoint(fields, adjoint);
    Complex left(0.0, 0.0);
    for (int axis = 0; axis != 3; ++axis)
    {
        left += dot(gradients[axis], fields[axis], fine);
    }
    EXPECT_NEAR(std::abs(left - dot(input, adjoint, coarse)), 0.0, 3.0e-11);
    Vector expected(size(coarse));
    for (int mx = -(coarse.nx - 1) / 2; mx <= coarse.nx / 2; ++mx)
    {
        for (int my = -(coarse.ny - 1) / 2; my <= coarse.ny / 2; ++my)
        {
            for (int mz = -(coarse.nz - 1) / 2; mz <= coarse.nz / 2; ++mz)
            {
                const std::array<int, 3> frequencies{mx, my, mz};
                const auto q = wavevector(coarse, frequencies);
                const Vector fine_mode = mode(fine, frequencies);
                const Vector coarse_mode = mode(coarse, frequencies);
                Complex coefficient(0.0, 0.0);
                for (int axis = 0; axis != 3; ++axis)
                {
                    coefficient += Complex(0.0, -q[axis]) * dot(fine_mode, fields[axis], fine) / volume(fine);
                }
                for (std::size_t i = 0; i != expected.size(); ++i)
                {
                    expected[i] += coefficient * coarse_mode[i];
                }
            }
        }
    }
    expect_near(adjoint, expected);
}

TEST(SternheimerGridTransfer, SameGridGradientsAreDerivativesWhileValuesRemainExact)
{
    const Grid grid = make_grid({3, 4, 2}, {0.21, -0.17, 0.08});
    const Vector input = samples(grid, 0.31);
    const std::array<Vector, 3> fields{{samples(grid, 0.1), samples(grid, -0.3), samples(grid, 0.7)}};
    Transfer transfer(grid, grid);
    Vector values;
    Vector adjoint;
    std::array<Vector, 3> gradients;
    transfer.interpolate_with_gradients(input, values, gradients);
    EXPECT_EQ(values, input);
    const auto expected = direct_gradients(grid, grid, input);
    Complex left(0.0, 0.0);
    for (int axis = 0; axis != 3; ++axis)
    {
        expect_near(gradients[axis], expected[axis]);
        left += dot(gradients[axis], fields[axis], grid);
    }
    transfer.restrict_gradient_adjoint(fields, adjoint);
    EXPECT_NEAR(std::abs(left - dot(input, adjoint, grid)), 0.0, 3.0e-11);
    transfer.interpolate(input, values);
    EXPECT_EQ(values, input);
    transfer.restrict_adjoint(input, values);
    EXPECT_EQ(values, input);
}

TEST(SternheimerGridTransfer, GradientNyquistUsesPositiveComplexModeWithoutTRClaim)
{
    const Grid coarse = make_grid({4, 6, 4}, {0.21, 0.08, -0.17});
    const Grid fine = make_grid({9, 10, 7}, coarse.kpoint);
    const std::array<int, 3> frequencies{2, 3, 2};
    const auto q = wavevector(coarse, frequencies);
    Transfer transfer(coarse, fine);
    Vector values;
    std::array<Vector, 3> gradients;
    transfer.interpolate_with_gradients(mode(coarse, frequencies, 0.63), values, gradients);
    const Vector fine_mode = mode(fine, frequencies, 0.63);
    for (int axis = 0; axis != 3; ++axis)
    {
        Vector expected = fine_mode;
        for (auto& value : expected)
        {
            value *= Complex(0.0, q[axis]);
        }
        expect_near(gradients[axis], expected, 5.0e-12);
    }
    std::array<Vector, 3> fields{{mode(fine, {-2, 0, 0}), mode(fine, {3, 0, 0}), mode(fine, {0, -4, 0})}};
    Vector adjoint;
    transfer.restrict_gradient_adjoint(fields, adjoint);
    expect_near(adjoint, Vector(size(coarse)));
}

TEST(SternheimerGridTransfer, GradientIntegerShiftedBlochLabelsPreserveCommonPhysicalMode)
{
    const Grid coarse = make_grid({6, 6, 4}, {0.17, -0.11, 0.08});
    const Grid fine = make_grid({11, 10, 9}, coarse.kpoint);
    Grid shifted_coarse = coarse;
    shifted_coarse.kpoint[0] += 1.0;
    shifted_coarse.kpoint[1] -= 1.0;
    Grid shifted_fine = fine;
    shifted_fine.kpoint = shifted_coarse.kpoint;
    Transfer original(coarse, fine);
    Transfer shifted(shifted_coarse, shifted_fine);
    Vector values;
    Vector shifted_values;
    std::array<Vector, 3> gradients;
    std::array<Vector, 3> shifted_gradients;
    original.interpolate_with_gradients(mode(coarse, {1, -1, 1}, 0.29), values, gradients);
    shifted.interpolate_with_gradients(mode(shifted_coarse, {0, 0, 1}, 0.29), shifted_values, shifted_gradients);
    expect_near(shifted_values, values);
    for (int axis = 0; axis != 3; ++axis)
    {
        expect_near(shifted_gradients[axis], gradients[axis]);
    }
}

TEST(SternheimerGridTransfer, GradientAPIsSupportInputOutputAliasingButRejectOverlappingOutputs)
{
    for (const bool same_grid : {false, true})
    {
        const Grid coarse = make_grid({3, 4, 2}, {0.17, -0.11, 0.08});
        const Grid fine = same_grid ? coarse : make_grid({5, 7, 4}, coarse.kpoint);
        Transfer transfer(coarse, fine);
        const Vector input = samples(coarse, 0.37);
        Vector values;
        std::array<Vector, 3> expected;
        transfer.interpolate_with_gradients(input, values, expected);
        Vector in_place = input;
        std::array<Vector, 3> gradients;
        transfer.interpolate_with_gradients(in_place, in_place, gradients);
        expect_near(in_place, values);
        gradients[1] = input;
        transfer.interpolate_with_gradients(gradients[1], in_place, gradients);
        expect_near(in_place, values);
        for (int axis = 0; axis != 3; ++axis)
        {
            expect_near(gradients[axis], expected[axis]);
        }
        Vector expected_adjoint;
        transfer.restrict_gradient_adjoint(expected, expected_adjoint);
        transfer.restrict_gradient_adjoint(gradients, gradients[2]);
        expect_near(gradients[2], expected_adjoint);
        const auto snapshot = gradients;
        EXPECT_THROW(transfer.interpolate_with_gradients(input, gradients[0], gradients), std::invalid_argument);
        EXPECT_EQ(gradients, snapshot);
    }
}

TEST(SternheimerGridTransfer, GradientValidationPrecedesAnyOutputChanges)
{
    for (const bool same_grid : {false, true})
    {
        const Grid coarse = make_grid({3, 4, 2});
        const Grid fine = same_grid ? coarse : make_grid({5, 7, 4});
        Transfer transfer(coarse, fine);
        const Vector sentinel(2, Complex(2.0, -3.0));
        Vector values = sentinel;
        std::array<Vector, 3> gradients{{sentinel, sentinel, sentinel}};
        EXPECT_THROW(transfer.interpolate_with_gradients(Vector(size(coarse) - 1), values, gradients),
                     std::invalid_argument);
        Vector bad = samples(coarse, 0.1);
        bad[0] = Complex(std::numeric_limits<double>::quiet_NaN(), 0.0);
        EXPECT_THROW(transfer.interpolate_with_gradients(bad, values, gradients), std::invalid_argument);
        std::array<Vector, 3> fields{{samples(fine, 0.1), samples(fine, 0.2), samples(fine, 0.3)}};
        fields[2].pop_back();
        EXPECT_THROW(transfer.restrict_gradient_adjoint(fields, values), std::invalid_argument);
        fields[2] = samples(fine, 0.3);
        fields[2][0] = Complex(0.0, std::numeric_limits<double>::infinity());
        EXPECT_THROW(transfer.restrict_gradient_adjoint(fields, values), std::invalid_argument);
        EXPECT_EQ(values, sentinel);
        EXPECT_EQ(gradients, (std::array<Vector, 3>{{sentinel, sentinel, sentinel}}));
    }
}

TEST(SternheimerGridTransfer, GradientWorkspaceIncludesLazyIdentityBuffersAndEnforcesBudget)
{
    const Grid coarse = make_grid({3, 4, 2}, {0.17, -0.11, 0.08});
    for (const bool same_grid : {false, true})
    {
        const Grid fine = same_grid ? coarse : make_grid({5, 7, 4}, coarse.kpoint);
        const std::size_t expected = (size(coarse) + size(fine)) * (2 * sizeof(double) + sizeof(Complex))
                                     + size(coarse) * sizeof(std::size_t);
        EXPECT_EQ(Transfer::gradient_workspace_bytes_required(coarse, fine), expected);
        EXPECT_EQ(Transfer::workspace_bytes_required(coarse, fine), same_grid ? 0U : expected);
        Transfer transfer(coarse, fine);
        const Vector sentinel(2, Complex(1.0, 2.0));
        Vector values = sentinel;
        std::array<Vector, 3> gradients{{sentinel, sentinel, sentinel}};
        const Vector input = samples(coarse, 0.1);
        EXPECT_THROW(transfer.interpolate_with_gradients(input, values, gradients, expected - 1), std::length_error);
        EXPECT_EQ(values, sentinel);
        EXPECT_EQ(gradients, (std::array<Vector, 3>{{sentinel, sentinel, sentinel}}));
        transfer.interpolate_with_gradients(input, values, gradients, expected);
        EXPECT_THROW(transfer.restrict_gradient_adjoint(gradients, values, expected - 1), std::length_error);
        EXPECT_NO_THROW(transfer.restrict_gradient_adjoint(gradients, values, expected));
    }
    Grid invalid = coarse;
    invalid.periodic = false;
    EXPECT_THROW(Transfer::gradient_workspace_bytes_required(coarse, invalid), std::invalid_argument);
}

TEST(SternheimerGridTransfer, GradientContextsAreIndependentAcrossWorkers)
{
    std::array<std::array<Vector, 3>, 4> results;
#ifdef _OPENMP
#pragma omp parallel for num_threads(4)
#endif
    for (int worker = 0; worker < 4; ++worker)
    {
        const Grid coarse = make_grid({3, 4, 2}, {0.13 * worker, -0.11, 0.08});
        const Grid fine = worker % 2 == 0 ? coarse : make_grid({5, 7, 4}, coarse.kpoint);
        Transfer transfer(coarse, fine);
        Vector values;
        transfer.interpolate_with_gradients(samples(coarse, 0.31), values, results[worker]);
    }
    for (int worker = 0; worker != 4; ++worker)
    {
        const Grid coarse = make_grid({3, 4, 2}, {0.13 * worker, -0.11, 0.08});
        const Grid fine = worker % 2 == 0 ? coarse : make_grid({5, 7, 4}, coarse.kpoint);
        const auto expected = direct_gradients(coarse, fine, samples(coarse, 0.31));
        for (int axis = 0; axis != 3; ++axis)
        {
            expect_near(results[worker][axis], expected[axis]);
        }
    }
}

TEST(SternheimerGridTransfer, SpectralNegativeLaplacianEqualsGradientAdjointGradient)
{
    const Grid coarse = make_grid({3, 4, 2}, {0.21, -0.17, 0.08});
    for (const bool same_grid : {false, true})
    {
        const Grid fine = same_grid ? coarse : make_grid({5, 7, 4}, coarse.kpoint);
        Transfer transfer(coarse, fine);
        const Vector input = samples(coarse, 0.31);
        Vector values;
        Vector actual;
        Vector composed;
        std::array<Vector, 3> gradients;
        transfer.interpolate_with_gradients(input, values, gradients);
        transfer.restrict_gradient_adjoint(gradients, composed);
        transfer.apply_negative_laplacian(input, actual);
        expect_near(actual, composed, 6.0e-12);
        Complex gradient_form(0.0, 0.0);
        for (const auto& gradient : gradients)
        {
            gradient_form += dot(gradient, gradient, fine);
        }
        EXPECT_NEAR(std::abs(dot(input, actual, coarse) - gradient_form), 0.0, 3.0e-10);
        Vector in_place = input;
        transfer.apply_negative_laplacian(in_place, in_place);
        expect_near(in_place, actual);
    }
}

TEST(SternheimerGridTransfer, SpectralNegativeLaplacianHasAnalyticBlochNyquistEigenvalue)
{
    const Grid coarse = make_grid({4, 6, 4}, {0.21, 0.08, -0.17});
    const Grid fine = make_grid({9, 10, 7}, coarse.kpoint);
    const std::array<int, 3> frequencies{2, 3, 2};
    const auto q = wavevector(coarse, frequencies);
    Vector expected = mode(coarse, frequencies, 0.63);
    const Vector input = expected;
    for (auto& value : expected)
    {
        value *= q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
    }
    Transfer transfer(coarse, fine);
    Vector actual;
    transfer.apply_negative_laplacian(input, actual);
    expect_near(actual, expected, 8.0e-12);
    const auto required = Transfer::gradient_workspace_bytes_required(coarse, fine);
    EXPECT_THROW(transfer.apply_negative_laplacian(input, actual, required - 1), std::length_error);
    EXPECT_NO_THROW(transfer.apply_negative_laplacian(input, actual, required));
    const Vector snapshot = actual;
    EXPECT_THROW(transfer.apply_negative_laplacian(Vector(), actual), std::invalid_argument);
    Vector bad = input;
    bad[0] = Complex(std::numeric_limits<double>::infinity(), 0.0);
    EXPECT_THROW(transfer.apply_negative_laplacian(bad, actual), std::invalid_argument);
    EXPECT_EQ(actual, snapshot);
}

TEST(SternheimerGridTransfer, SpectralOperationsRejectArithmeticOverflow)
{
    const Grid coarse = make_grid({4, 4, 4}, {0.0, 0.0, 0.0});
    const Grid fine = make_grid({6, 6, 6}, coarse.kpoint);
    Transfer transfer(coarse, fine);
    const Vector input(size(coarse), Complex(1.0e308, 0.0));
    Vector output;
    std::array<Vector, 3> gradients;
    EXPECT_THROW(transfer.interpolate_with_gradients(input, output, gradients), std::overflow_error);
    EXPECT_THROW(transfer.apply_negative_laplacian(input, output), std::overflow_error);
    const Vector field(size(fine), Complex(1.0e308, 0.0));
    const std::array<Vector, 3> fields{{field, field, field}};
    EXPECT_THROW(transfer.restrict_gradient_adjoint(fields, output), std::overflow_error);
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
