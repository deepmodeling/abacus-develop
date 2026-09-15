#include "source_lcao/module_ri/sternheimer_weak_preconditioner.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{

using Preconditioner = ModuleRI::SternheimerWeakSpectralPreconditioner;
using Grid = Preconditioner::Grid;
using Vector = Preconditioner::Vector;
using Complex = Preconditioner::Complex;
using Lattice = std::array<std::array<double, 3>, 3>;

constexpr double alpha = 1.0;
constexpr double epsilon = 0.37;
constexpr double omega = 0.29;
constexpr double regularization = 0.02;

Grid make_grid(const std::array<int, 3>& dimensions,
               const Lattice& lattice,
               const std::array<double, 3>& kpoint)
{
    Grid grid;
    grid.nx = dimensions[0];
    grid.ny = dimensions[1];
    grid.nz = dimensions[2];
    grid.periodic = true;
    grid.lattice_vectors = lattice;
    grid.kpoint = kpoint;
    const auto length = [](const std::array<double, 3>& value) {
        return std::sqrt(value[0] * value[0] + value[1] * value[1]
                         + value[2] * value[2]);
    };
    grid.hx = length(lattice[0]) / grid.nx;
    grid.hy = length(lattice[1]) / grid.ny;
    grid.hz = length(lattice[2]) / grid.nz;
    return grid;
}

Vector mode(const Grid& grid, const std::array<int, 3>& frequencies,
            const double phase = 0.0)
{
    Vector result(static_cast<std::size_t>(grid.size()));
    const double two_pi = 2.0 * std::acos(-1.0);
    for (int ix = 0; ix != grid.nx; ++ix)
    {
        for (int iy = 0; iy != grid.ny; ++iy)
        {
            for (int iz = 0; iz != grid.nz; ++iz)
            {
                const double angle
                    = phase
                      + two_pi
                            * ((frequencies[0] + grid.kpoint[0]) * ix / grid.nx
                               + (frequencies[1] + grid.kpoint[1]) * iy / grid.ny
                               + (frequencies[2] + grid.kpoint[2]) * iz / grid.nz);
                const auto index
                    = (static_cast<std::size_t>(ix) * grid.ny + iy) * grid.nz + iz;
                result[index] = std::exp(Complex(0.0, angle));
            }
        }
    }
    return result;
}

void expect_scaled(const Vector& actual,
                   const Vector& input,
                   const Complex factor,
                   const double tolerance = 3.0e-12)
{
    ASSERT_EQ(actual.size(), input.size());
    for (std::size_t i = 0; i != input.size(); ++i)
    {
        EXPECT_NEAR(std::abs(actual[i] - factor * input[i]), 0.0, tolerance)
            << "index " << i;
    }
}

double squared_norm(const std::array<double, 3>& value)
{
    return value[0] * value[0] + value[1] * value[1] + value[2] * value[2];
}

std::array<double, 3> orthogonal_wavevector(const Grid& grid,
                                            const std::array<int, 3>& frequencies)
{
    const double two_pi = 2.0 * std::acos(-1.0);
    return {two_pi * (frequencies[0] + grid.kpoint[0]) / 6.0,
            two_pi * (frequencies[1] + grid.kpoint[1]) / 7.0,
            two_pi * (frequencies[2] + grid.kpoint[2]) / 8.0};
}

std::array<double, 3> skew_wavevector(const Grid& grid,
                                      const std::array<int, 3>& frequencies)
{
    const double two_pi = 2.0 * std::acos(-1.0);
    const double qx = two_pi * (frequencies[0] + grid.kpoint[0]) / 5.0;
    const double qy
        = (two_pi * (frequencies[1] + grid.kpoint[1]) + 2.5 * qx) / 4.3;
    const double qz
        = (two_pi * (frequencies[2] + grid.kpoint[2]) - 0.4 * qx + 0.2 * qy)
          / 8.0;
    return {qx, qy, qz};
}

Complex inverse_denominator(const std::array<double, 3>& wavevector,
                            const double signed_omega = omega)
{
    return 1.0
           / Complex(alpha * squared_norm(wavevector) - epsilon + regularization,
                     signed_omega);
}

} // namespace

TEST(SternheimerWeakPreconditioner, MatchesOrthogonalBlochFourierMode)
{
    const Lattice lattice{{{6.0, 0.0, 0.0},
                           {0.0, 7.0, 0.0},
                           {0.0, 0.0, 8.0}}};
    const Grid grid = make_grid({5, 4, 6}, lattice, {0.17, -0.11, 0.08});
    const std::array<int, 3> frequencies{2, -1, 2};
    const Vector input = mode(grid, frequencies, 0.43);
    const Preconditioner preconditioner(grid, alpha, epsilon, omega, regularization);

    Vector actual;
    preconditioner.apply(input, actual);

    expect_scaled(actual, input,
                  inverse_denominator(orthogonal_wavevector(grid, frequencies)));
}

TEST(SternheimerWeakPreconditioner, MatchesSkewBlochFourierMode)
{
    const Lattice lattice{{{5.0, 0.0, 0.0},
                           {-2.5, 4.3, 0.0},
                           {0.4, -0.2, 8.0}}};
    const Grid grid = make_grid({5, 4, 3}, lattice, {-0.21, 0.13, 0.07});
    const std::array<int, 3> frequencies{-2, 1, -1};
    const Vector input = mode(grid, frequencies, -0.31);
    const Preconditioner preconditioner(grid, alpha, epsilon, omega, regularization);

    Vector actual;
    preconditioner.apply(input, actual);

    expect_scaled(actual, input, inverse_denominator(skew_wavevector(grid, frequencies)));
}

TEST(SternheimerWeakPreconditioner, SignedFrequencyConjugacyAndInPlaceApply)
{
    const Lattice lattice{{{6.0, 0.0, 0.0},
                           {0.0, 7.0, 0.0},
                           {0.0, 0.0, 8.0}}};
    const Grid grid = make_grid({4, 4, 4}, lattice, {0.0, 0.0, 0.0});
    Vector positive = mode(grid, {0, 0, 0});
    Vector negative = positive;
    const Preconditioner plus(grid, alpha, -0.3, omega, 0.0);
    const Preconditioner minus(grid, alpha, -0.3, -omega, 0.0);

    plus.apply(positive, positive);
    minus.apply(negative, negative);

    ASSERT_EQ(positive.size(), negative.size());
    for (std::size_t i = 0; i != positive.size(); ++i)
    {
        EXPECT_NEAR(std::abs(negative[i] - std::conj(positive[i])), 0.0, 2.0e-12);
    }
}

TEST(SternheimerWeakPreconditioner, RejectsInvalidConstruction)
{
    const Lattice lattice{{{6.0, 0.0, 0.0},
                           {0.0, 7.0, 0.0},
                           {0.0, 0.0, 8.0}}};
    Grid grid = make_grid({4, 4, 4}, lattice, {0.0, 0.0, 0.0});
    EXPECT_THROW(Preconditioner(grid, alpha, 0.0, 0.0, 0.0), std::domain_error);
    EXPECT_THROW(Preconditioner(grid, 0.0, epsilon, omega, regularization),
                 std::invalid_argument);
    EXPECT_THROW(Preconditioner(grid, alpha,
                                std::numeric_limits<double>::quiet_NaN(), omega,
                                regularization),
                 std::invalid_argument);
    EXPECT_THROW(Preconditioner(grid, alpha, epsilon,
                                std::numeric_limits<double>::infinity(),
                                regularization),
                 std::invalid_argument);
    EXPECT_THROW(Preconditioner(grid, alpha, epsilon, omega, -0.1),
                 std::invalid_argument);
    grid.periodic = false;
    EXPECT_THROW(Preconditioner(grid, alpha, epsilon, omega, regularization),
                 std::invalid_argument);
    grid.periodic = true;
    grid.nx = 0;
    EXPECT_THROW(Preconditioner(grid, alpha, epsilon, omega, regularization),
                 std::invalid_argument);
}

TEST(SternheimerWeakPreconditioner, RejectsInvalidInputWithoutChangingOutput)
{
    const Lattice lattice{{{6.0, 0.0, 0.0},
                           {0.0, 7.0, 0.0},
                           {0.0, 0.0, 8.0}}};
    const Grid grid = make_grid({4, 4, 4}, lattice, {0.0, 0.0, 0.0});
    const Preconditioner preconditioner(grid, alpha, epsilon, omega, regularization);
    Vector output{Complex(7.0, -3.0)};
    EXPECT_THROW(preconditioner.apply(Vector(3, 0.0), output),
                 std::invalid_argument);
    EXPECT_EQ(output, Vector({Complex(7.0, -3.0)}));
    Vector input(static_cast<std::size_t>(grid.size()), 0.0);
    input[7] = Complex(std::numeric_limits<double>::quiet_NaN(), 0.0);
    EXPECT_THROW(preconditioner.apply(input, output), std::invalid_argument);
    EXPECT_EQ(output, Vector({Complex(7.0, -3.0)}));
}
