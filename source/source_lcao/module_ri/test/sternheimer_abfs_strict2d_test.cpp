#include <gtest/gtest.h>

#if __has_include("source_lcao/module_ri/sternheimer_abfs_strict2d.h")
#include "source_lcao/module_ri/sternheimer_abfs_strict2d.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

namespace
{
using Complex = std::complex<double>;
using Grid = ModuleRI::SternheimerFDHamiltonian::Grid;
using Channel = ModuleRI::SternheimerABFBlochGridChannel;
using QPoint = ModuleRI::SternheimerReducedKPoint;
constexpr double two_pi = 6.283185307179586476925286766559;

Grid skew_grid(const int nx = 4, const int ny = 3, const int nz = 5)
{
    Grid grid;
    grid.nx = nx;
    grid.ny = ny;
    grid.nz = nz;
    grid.lattice_vectors = {{{4.2, 0.8, 0.0}, {1.1, 3.7, 0.0}, {0.0, 0.0, 6.5}}};
    grid.hx = std::hypot(4.2, 0.8) / nx;
    grid.hy = std::hypot(1.1, 3.7) / ny;
    grid.hz = 6.5 / nz;
    return grid;
}

int index(const Grid& grid, const int x, const int y, const int z)
{
    return (x * grid.ny + y) * grid.nz + z;
}

// Independent Cartesian 2x2 inverse, not the production grid-dual helper.
double wave_number(const Grid& grid, const double gx, const double gy)
{
    const auto& a = grid.lattice_vectors[0];
    const auto& b = grid.lattice_vectors[1];
    const double determinant = a[0] * b[1] - a[1] * b[0];
    return two_pi * std::hypot((b[1] * gx - a[1] * gy) / determinant,
                              (a[0] * gy - b[0] * gx) / determinant);
}

Channel density(const Grid& grid, const int seed)
{
    Channel channel;
    channel.channel_index = seed;
    channel.atom_index = 7;
    channel.label = "complex_fine_grid_density";
    channel.max_abs = -1.0;
    channel.potential_r.resize(grid.nx * grid.ny * grid.nz);
    for (std::size_t i = 0; i < channel.potential_r.size(); ++i)
    {
        const double t = static_cast<double>(i + 1);
        channel.potential_r[i] = Complex(std::sin((0.31 + seed * 0.17) * t),
                                         std::cos((0.47 + seed * 0.11) * t));
    }
    return channel;
}

// Explicit planar Fourier projection of FULL Bloch fields and direct z sum.
// This oracle uses neither FFTW nor recurrences nor production geometry helpers.
std::vector<Complex> direct_potential(const Channel& rho, const Grid& grid, const QPoint& q)
{
    std::vector<Complex> result(rho.potential_r.size(), Complex{});
    const double dz = std::abs(grid.lattice_vectors[2][2]) / grid.nz;
    for (int mx = -grid.nx / 2; mx <= (grid.nx - 1) / 2; ++mx)
    {
        for (int my = -grid.ny / 2; my <= (grid.ny - 1) / 2; ++my)
        {
            const double Q = wave_number(grid, mx + q[0], my + q[1]);
            std::vector<Complex> coefficients(grid.nz, Complex{});
            for (int z = 0; z < grid.nz; ++z)
            {
                for (int x = 0; x < grid.nx; ++x)
                {
                    for (int y = 0; y < grid.ny; ++y)
                    {
                        const double angle = two_pi * ((mx + q[0]) * x / grid.nx
                                                       + (my + q[1]) * y / grid.ny);
                        coefficients[z] += rho.potential_r[index(grid, x, y, z)]
                                           * std::exp(Complex(0.0, -angle)) / double(grid.nx * grid.ny);
                    }
                }
            }
            for (int z = 0; z < grid.nz; ++z)
            {
                Complex potential{};
                for (int source_z = 0; source_z < grid.nz; ++source_z)
                {
                    potential += two_pi * dz / Q * std::exp(-Q * dz * std::abs(z - source_z))
                                 * coefficients[source_z];
                }
                for (int x = 0; x < grid.nx; ++x)
                {
                    for (int y = 0; y < grid.ny; ++y)
                    {
                        const double angle = two_pi * ((mx + q[0]) * x / grid.nx
                                                       + (my + q[1]) * y / grid.ny);
                        result[index(grid, x, y, z)] += potential * std::exp(Complex(0.0, angle));
                    }
                }
            }
        }
    }
    return result;
}

Complex inner_product(const Channel& a, const Channel& b, const double dv)
{
    Complex result{};
    for (std::size_t i = 0; i < a.potential_r.size(); ++i)
    {
        result += std::conj(a.potential_r[i]) * b.potential_r[i] * dv;
    }
    return result;
}
} // namespace

TEST(SternheimerABFSStrict2D, MatchesDirectDiscreteSumOnSkewOddEvenGrids)
{
    for (const Grid grid: {skew_grid(), skew_grid(3, 4, 1), skew_grid(1, 1, 3)})
    {
        // Includes negative q, a zone edge, and an unfolded reciprocal shift.
        for (const QPoint q: {QPoint{0.0, 0.25, 0.0}, QPoint{0.5, -0.23, 0.0}, QPoint{1.17, -1.21, 0.0}})
        {
            std::vector<Channel> channels{density(grid, 0), density(grid, 1)};
            const auto original = channels;
            ModuleRI::solve_sternheimer_abf_strict2d_coulomb_in_place(channels, grid, q);
            for (std::size_t c = 0; c < channels.size(); ++c)
            {
                const auto expected = direct_potential(original[c], grid, q);
                double max_abs = 0.0;
                for (std::size_t i = 0; i < expected.size(); ++i)
                {
                    EXPECT_NEAR(std::abs(channels[c].potential_r[i] - expected[i]), 0.0,
                                3.0e-12 * std::max(1.0, std::abs(expected[i])));
                    max_abs = std::max(max_abs, std::abs(expected[i]));
                }
                EXPECT_NEAR(channels[c].max_abs, max_abs, 3.0e-12 * std::max(1.0, max_abs));
                EXPECT_EQ(channels[c].channel_index, original[c].channel_index);
                EXPECT_EQ(channels[c].atom_index, original[c].atom_index);
                EXPECT_EQ(channels[c].label, original[c].label);
            }
        }
    }
}

TEST(SternheimerABFSStrict2D, IsHermitianPositiveAndComplexLinear)
{
    const Grid grid = skew_grid();
    const QPoint q{0.13, -0.21, 0.0};
    const Complex alpha(0.7, -0.3), beta(-0.4, 0.8);
    std::vector<Channel> rho{density(grid, 0), density(grid, 1), density(grid, 2)};
    for (std::size_t i = 0; i < rho[0].potential_r.size(); ++i)
    {
        rho[2].potential_r[i] = alpha * rho[0].potential_r[i] + beta * rho[1].potential_r[i];
    }
    auto phi = rho;
    ModuleRI::solve_sternheimer_abf_strict2d_coulomb_in_place(phi, grid, q);
    const Complex ab = inner_product(rho[0], phi[1], 1.0);
    const Complex ba = inner_product(rho[1], phi[0], 1.0);
    EXPECT_NEAR(std::abs(ab - std::conj(ba)), 0.0, 2.0e-11);
    for (std::size_t c = 0; c < rho.size(); ++c)
    {
        const Complex norm = inner_product(rho[c], phi[c], 1.0);
        EXPECT_GT(norm.real(), 0.0);
        EXPECT_NEAR(norm.imag(), 0.0, 2.0e-11);
    }
    for (std::size_t i = 0; i < rho[0].potential_r.size(); ++i)
    {
        EXPECT_NEAR(std::abs(phi[2].potential_r[i] - alpha * phi[0].potential_r[i]
                             - beta * phi[1].potential_r[i]), 0.0, 2.0e-12);
    }
}

TEST(SternheimerABFSStrict2D, ImplicitOrthogonalCellMatchesExplicitCell)
{
    Grid implicit_grid = skew_grid();
    implicit_grid.lattice_vectors = {};
    Grid explicit_grid = implicit_grid;
    explicit_grid.lattice_vectors = {{{implicit_grid.nx * implicit_grid.hx, 0.0, 0.0},
                                      {0.0, implicit_grid.ny * implicit_grid.hy, 0.0},
                                      {0.0, 0.0, implicit_grid.nz * implicit_grid.hz}}};
    const QPoint q{0.0, 0.25, 0.0};
    const Channel rho = density(implicit_grid, 0);
    const auto expected = direct_potential(rho, explicit_grid, q);
    std::vector<Channel> phi{rho};
    ModuleRI::solve_sternheimer_abf_strict2d_coulomb_in_place(phi, implicit_grid, q);
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
        EXPECT_NEAR(std::abs(phi[0].potential_r[i] - expected[i]), 0.0, 3.0e-12);
    }
    std::vector<Channel> empty;
    EXPECT_NO_THROW(ModuleRI::solve_sternheimer_abf_strict2d_coulomb_in_place(empty, implicit_grid, q));
    EXPECT_TRUE(empty.empty());
    EXPECT_TRUE(ModuleRI::sternheimer_abf_strict2d_selected_coulomb_integrals(empty, empty, implicit_grid, {}).empty());
}

TEST(SternheimerABFSStrict2D, BoundarySheetHasEqualWeightAndNoPeriodicZWrap)
{
    const Grid grid = skew_grid();
    const QPoint q{0.12, -0.16, 0.0};
    const double Q = wave_number(grid, q[0], q[1]);
    const double dz = 6.5 / grid.nz;
    for (const int source_z: {0, grid.nz - 1})
    {
        Channel rho;
        rho.potential_r.resize(grid.nx * grid.ny * grid.nz, Complex{});
        const Complex amplitude(0.4, -0.9);
        for (int x = 0; x < grid.nx; ++x)
        {
            for (int y = 0; y < grid.ny; ++y)
            {
                rho.potential_r[index(grid, x, y, source_z)]
                    = amplitude * std::exp(Complex(0.0, two_pi * (q[0] * x / grid.nx + q[1] * y / grid.ny)));
            }
        }
        std::vector<Channel> phi{rho};
        ModuleRI::solve_sternheimer_abf_strict2d_coulomb_in_place(phi, grid, q);
        for (int z = 0; z < grid.nz; ++z)
        {
            const Complex expected = amplitude * two_pi * dz / Q * std::exp(-Q * dz * std::abs(z - source_z));
            EXPECT_NEAR(std::abs(phi[0].potential_r[index(grid, 0, 0, z)] - expected), 0.0, 2.0e-12);
        }
        EXPECT_LT(std::abs(phi[0].potential_r[index(grid, 0, 0, grid.nz - 1 - source_z)]),
                  std::abs(phi[0].potential_r[index(grid, 0, 0, source_z)]) * 0.5);
    }
}

TEST(SternheimerABFSStrict2D, SelectedIntegralsUseCellVolumeAndRequestedOrder)
{
    const Grid grid = skew_grid();
    const QPoint q{0.0, 0.25, 0.0};
    const std::vector<Channel> rho{density(grid, 0), density(grid, 1)};
    auto phi = rho;
    ModuleRI::solve_sternheimer_abf_strict2d_coulomb_in_place(phi, grid, q);
    const std::vector<std::pair<std::size_t, std::size_t>> pairs{{1, 0}, {0, 0}, {0, 1}, {1, 0}};
    const auto selected = ModuleRI::sternheimer_abf_strict2d_selected_coulomb_integrals(rho, phi, grid, pairs);
    const double dv = (4.2 * 3.7 - 0.8 * 1.1) * 6.5 / (grid.nx * grid.ny * grid.nz);
    ASSERT_EQ(selected.size(), pairs.size());
    for (std::size_t i = 0; i < pairs.size(); ++i)
    {
        EXPECT_NEAR(std::abs(selected[i] - inner_product(rho[pairs[i].first], phi[pairs[i].second], dv)),
                    0.0, 2.0e-11);
    }
    EXPECT_NEAR(std::abs(selected[0] - std::conj(selected[2])), 0.0, 2.0e-11);
    EXPECT_THROW(ModuleRI::sternheimer_abf_strict2d_selected_coulomb_integrals(rho, phi, grid, {{2, 0}}),
                 std::invalid_argument);
    auto sparse_phi = phi;
    sparse_phi[1].potential_r.clear();
    EXPECT_NO_THROW(ModuleRI::sternheimer_abf_strict2d_selected_coulomb_integrals(rho, sparse_phi, grid, {{1, 0}}));
    EXPECT_THROW(ModuleRI::sternheimer_abf_strict2d_selected_coulomb_integrals(rho, sparse_phi, grid, {{0, 1}}),
                 std::invalid_argument);
}

TEST(SternheimerABFSStrict2D, RejectsGammaQzAndUnsupportedCellsBeforeMutation)
{
    const Grid grid = skew_grid();
    const Channel original = density(grid, 0);
    std::vector<Channel> channels{original};
    for (const QPoint q: {QPoint{0.0, 0.0, 0.0}, QPoint{1.0, -1.0, 0.0}, QPoint{0.1, 0.2, 0.01},
                          QPoint{0.1, 0.2, 1.0}, QPoint{0.1, 0.2, 1.0e-15},
                          QPoint{std::numeric_limits<double>::quiet_NaN(), 0.2, 0.0}})
    {
        EXPECT_THROW(ModuleRI::solve_sternheimer_abf_strict2d_coulomb_in_place(channels, grid, q),
                     std::invalid_argument);
        EXPECT_EQ(channels[0].potential_r, original.potential_r);
    }
    const QPoint q{0.1, 0.2, 0.0};
    for (int defect = 0; defect < 7; ++defect)
    {
        Grid invalid = grid;
        if (defect == 0) invalid.lattice_vectors[2][0] = 0.1;
        if (defect == 1) invalid.lattice_vectors[0][2] = 0.1;
        if (defect == 2) invalid.lattice_vectors[1] = invalid.lattice_vectors[0];
        if (defect == 3) invalid.periodic = false;
        if (defect == 4) invalid.nz = 0;
        if (defect == 5) invalid.lattice_vectors[2][2] = 0.0;
        if (defect == 6) invalid.lattice_vectors[0][0] = std::numeric_limits<double>::infinity();
        EXPECT_THROW(ModuleRI::solve_sternheimer_abf_strict2d_coulomb_in_place(channels, invalid, q),
                     std::invalid_argument);
        EXPECT_EQ(channels[0].potential_r, original.potential_r);
    }
    channels.push_back(original);
    channels[1].potential_r.pop_back();
    EXPECT_THROW(ModuleRI::solve_sternheimer_abf_strict2d_coulomb_in_place(channels, grid, q), std::invalid_argument);
    EXPECT_EQ(channels[0].potential_r, original.potential_r);
    channels[1] = original;
    channels[1].potential_r[0] = Complex(0.0, std::numeric_limits<double>::quiet_NaN());
    EXPECT_THROW(ModuleRI::solve_sternheimer_abf_strict2d_coulomb_in_place(channels, grid, q), std::invalid_argument);
    EXPECT_EQ(channels[0].potential_r, original.potential_r);
    EXPECT_THROW(ModuleRI::sternheimer_abf_strict2d_selected_coulomb_integrals(channels, channels, grid, {{0, 1}}),
                 std::invalid_argument);
}

TEST(SternheimerABFSStrict2D, RejectsFiniteDensityOverflowInForwardRecurrence)
{
    Grid grid;
    grid.nx = grid.ny = 1;
    grid.nz = 3;
    grid.hx = grid.hy = grid.hz = 1.0;
    Channel rho;
    rho.potential_r.assign(3, Complex(0.75 * std::numeric_limits<double>::max(), 0.0));
    rho.max_abs = 17.0;
    std::vector<Channel> channels{rho};
    try
    {
        ModuleRI::solve_sternheimer_abf_strict2d_coulomb_in_place(channels, grid, QPoint{0.001, 0.0, 0.0});
        FAIL() << "Finite density must not silently overflow the z recurrence.";
    }
    catch (const std::overflow_error& error)
    {
        EXPECT_NE(std::string(error.what()).find("forward recurrence"), std::string::npos);
    }
    EXPECT_EQ(channels[0].potential_r, rho.potential_r);
    EXPECT_EQ(channels[0].max_abs, rho.max_abs);
}

TEST(SternheimerABFSStrict2D, RejectsFiniteModeOverflowInInverseFFT)
{
    Grid grid;
    grid.nx = 2;
    grid.ny = grid.nz = 1;
    grid.hx = grid.hy = grid.hz = 1.0;
    Channel rho;
    // q=1/2 gives equal factors of two for both planar FFT modes. Each
    // scaled mode is finite (0.8*max), but their inverse-FFT sum overflows.
    rho.potential_r = {Complex(0.4 * std::numeric_limits<double>::max(), 0.0), Complex{}};
    rho.max_abs = 19.0;
    std::vector<Channel> channels{rho};
    try
    {
        ModuleRI::solve_sternheimer_abf_strict2d_coulomb_in_place(channels, grid, QPoint{0.5, 0.0, 0.0});
        FAIL() << "Finite Fourier coefficients must not silently overflow the inverse FFT.";
    }
    catch (const std::overflow_error& error)
    {
        EXPECT_NE(std::string(error.what()).find("inverse FFT"), std::string::npos);
    }
    EXPECT_EQ(channels[0].potential_r, rho.potential_r);
    EXPECT_EQ(channels[0].max_abs, rho.max_abs);
}

TEST(SternheimerABFSStrict2D, RejectsFiniteFieldOverflowInSelectedIntegrals)
{
    for (int stage = 0; stage < 3; ++stage)
    {
        Grid grid;
        grid.nx = grid.ny = 1;
        grid.nz = stage == 1 ? 2 : 1;
        grid.hx = grid.hy = 1.0;
        grid.hz = stage == 2 ? 4.0 : 1.0;
        Channel rho, phi;
        rho.potential_r.assign(grid.nz, Complex(0.75 * std::numeric_limits<double>::max(), 0.0));
        phi.potential_r.assign(grid.nz, Complex(stage == 0 ? 2.0 : 1.0, 0.0));
        // Independently overflow the product, the running sum, and dV scaling.
        EXPECT_THROW(ModuleRI::sternheimer_abf_strict2d_selected_coulomb_integrals({rho}, {phi}, grid, {{0, 0}}),
                     std::overflow_error);
    }
}

TEST(SternheimerABFSStrict2D, AcceptsLargeFiniteRepresentableDensitiesAndIntegrals)
{
    const Grid grid = skew_grid();
    const QPoint q{0.0, 0.25, 0.0};
    const double scale = 1.0e140;
    const std::vector<Channel> rho{density(grid, 0)};
    auto large_rho = rho;
    for (auto& value: large_rho[0].potential_r) value *= scale;
    auto phi = rho;
    auto large_phi = large_rho;
    ASSERT_NO_THROW(ModuleRI::solve_sternheimer_abf_strict2d_coulomb_in_place(phi, grid, q));
    ASSERT_NO_THROW(ModuleRI::solve_sternheimer_abf_strict2d_coulomb_in_place(large_phi, grid, q));
    for (std::size_t i = 0; i < phi[0].potential_r.size(); ++i)
    {
        EXPECT_NEAR(std::abs(large_phi[0].potential_r[i] / scale - phi[0].potential_r[i]), 0.0, 3.0e-12);
    }
    EXPECT_TRUE(std::isfinite(large_phi[0].max_abs));
    const auto baseline = ModuleRI::sternheimer_abf_strict2d_selected_coulomb_integrals(rho, phi, grid, {{0, 0}});
    const auto large = ModuleRI::sternheimer_abf_strict2d_selected_coulomb_integrals(large_rho, large_phi, grid, {{0, 0}});
    ASSERT_EQ(large.size(), 1u);
    EXPECT_TRUE(std::isfinite(large[0].real()));
    EXPECT_TRUE(std::isfinite(large[0].imag()));
    EXPECT_NEAR(std::abs(large[0] / (scale * scale) - baseline[0]), 0.0, 3.0e-11);
}

#else
TEST(SternheimerABFSStrict2D, RequiresStrict2DCoulombImplementation)
{
    FAIL() << "Missing sternheimer_abfs_strict2d.h: discrete nonzero-q strict2D Coulomb is not implemented.";
}
#endif
