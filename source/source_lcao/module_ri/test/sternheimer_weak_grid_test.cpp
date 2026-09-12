#include "source_lcao/module_ri/sternheimer_weak_grid.h"

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

namespace
{
using H = ModuleRI::SternheimerFDHamiltonian;
using Operator = ModuleRI::SternheimerWeakGridOperator;
using Function = ModuleRI::SternheimerDeltaGridFunction;
using Complex = H::Complex;
using Vector = H::Vector;
constexpr double pi2 = 6.283185307179586476925286766559;

H::Grid Grid(int nx, int ny, int nz)
{
    H::Grid grid{nx, ny, nz, 4.0 / nx, std::sqrt(9.25) / ny, std::sqrt(6.5) / nz, true};
    grid.lattice_vectors = {{{4.0, 0.0, 0.0}, {0.5, 3.0, 0.0}, {0.3, 0.4, 2.5}}};
    grid.kpoint = {0.17, -0.13, 0.21};
    return grid;
}

std::array<double, 3> Position(const H::Grid& grid, int i)
{
    return {{double(i / (grid.ny * grid.nz)) / grid.nx,
             double((i / grid.nz) % grid.ny) / grid.ny, double(i % grid.nz) / grid.nz}};
}

double DV(const H::Grid& grid) { return 30.0 / grid.size(); }

Complex Dot(const Vector& a, const Vector& b)
{
    Complex value{};
    for (std::size_t i = 0; i < a.size(); ++i) value += std::conj(a[i]) * b[i];
    return value;
}

Function PlaneWave(const H::Grid& grid, std::array<int, 3> mode)
{
    const auto dual = ModuleRI::sternheimer_fd_grid_dual_vectors(grid);
    std::array<double, 3> momentum{};
    for (int a = 0; a < 3; ++a)
        for (int j = 0; j < 3; ++j) momentum[a] += pi2 * (mode[j] + grid.kpoint[j]) * dual[j][a];
    Function result;
    result.values.resize(grid.size());
    for (auto& gradient : result.gradients) gradient.resize(grid.size());
    for (int i = 0; i < grid.size(); ++i)
    {
        const auto r = Position(grid, i);
        double phase = 0.0;
        for (int a = 0; a < 3; ++a) phase += (mode[a] + grid.kpoint[a]) * r[a];
        result.values[i] = std::polar(1.0 / std::sqrt(30.0), pi2 * phase);
        for (int a = 0; a < 3; ++a) result.gradients[a][i] = Complex(0.0, momentum[a]) * result.values[i];
    }
    return result;
}

// Explicit Fourier cardinal functions, independent of the FFT transfer code.
std::vector<Function> ExplicitE(const H::Grid& coarse, const H::Grid& fine)
{
    std::vector<Function> result(coarse.size());
    const std::array<int, 3> dims{{coarse.nx, coarse.ny, coarse.nz}};
    for (int j = 0; j < coarse.size(); ++j)
    {
        auto& column = result[j];
        column.values.assign(fine.size(), {});
        for (auto& gradient : column.gradients) gradient.assign(fine.size(), {});
        const auto r = Position(coarse, j);
        for (int m = 0; m < coarse.size(); ++m)
        {
            std::array<int, 3> mode{{m / (coarse.ny * coarse.nz), (m / coarse.nz) % coarse.ny, m % coarse.nz}};
            double phase = 0.0;
            for (int a = 0; a < 3; ++a)
            {
                if (mode[a] > dims[a] / 2) mode[a] -= dims[a];
                phase += (mode[a] + coarse.kpoint[a]) * r[a];
            }
            const auto wave = PlaneWave(fine, mode);
            const Complex weight = std::polar(1.0 / std::sqrt(double(coarse.size())), -pi2 * phase);
            for (int i = 0; i < fine.size(); ++i)
            {
                column.values[i] += weight * wave.values[i];
                for (int a = 0; a < 3; ++a) column.gradients[a][i] += weight * wave.gradients[a][i];
            }
        }
    }
    return result;
}

std::shared_ptr<const H> MakeH(const H::Grid& grid)
{
    std::vector<double> potential(grid.size());
    for (int i = 0; i < grid.size(); ++i)
    {
        const auto r = Position(grid, i);
        potential[i] = 0.9 + 0.2 * std::cos(pi2 * r[0]) + 0.1 * std::sin(pi2 * (r[1] + r[2]));
    }
    using P = ModuleRI::SternheimerFDNonlocalProjector;
    P::ProjectorBlock block;
    block.projectors = {PlaneWave(grid, {0, 0, 0}).values, PlaneWave(grid, {1, 1, 0}).values};
    block.d_matrix = {{Complex(0.4), Complex(0.08, 0.07)}, {Complex(0.08, -0.07), Complex(0.3)}};
    auto projector = std::make_shared<const P>(grid.size(), DV(grid), std::vector<P::ProjectorBlock>{block});
    return std::make_shared<const H>(grid, potential, 0.7, projector, 8);
}

void Near(const Vector& a, const Vector& b, double tolerance = 2e-10)
{
    ASSERT_EQ(a.size(), b.size());
    for (std::size_t i = 0; i < a.size(); ++i) EXPECT_LT(std::abs(a[i] - b[i]), tolerance) << i;
}

TEST(SternheimerWeakGrid, CrossBlocksAndOperatorMatchExplicitFineWeakForm)
{
    const auto coarse = Grid(2, 3, 2), fine = Grid(5, 5, 4);
    const auto h = MakeH(fine);
    Operator op(h, coarse);
    const auto e = ExplicitE(coarse, fine);
    // The mode 6 aliases in values on nx=5, but its analytic gradient does not.
    std::vector<Function> u{PlaneWave(fine, {0, 0, 0}), PlaneWave(fine, {6, 0, 0}), PlaneWave(fine, {2, 1, -1})};
    const auto blocks = op.assemble_blocks(u);
    const auto expected_u = ModuleRI::assemble_delta_sternheimer_grid_matrices(*h, u, DV(fine));
    Near(blocks.state_hamiltonian, expected_u.hamiltonian);
    ASSERT_EQ(blocks.state_count, u.size());
    ASSERT_EQ(blocks.coarse_count, std::size_t(coarse.size()));
    for (int j = 0; j < coarse.size(); ++j)
    {
        Vector unit(coarse.size()), image, lifted;
        unit[j] = 1.0;
        op.lift(unit, lifted);
        Near(lifted, e[j].values);
        op.apply(unit, image);
        Vector local, nonlocal;
        h->apply_local_potential(e[j].values, local);
        h->apply_nonlocal(e[j].values, nonlocal);
        for (std::size_t i = 0; i < u.size(); ++i)
        {
            Complex k = Dot(u[i].values, local) + Dot(u[i].values, nonlocal);
            for (int a = 0; a < 3; ++a) k += h->kinetic_prefactor() * Dot(u[i].gradients[a], e[j].gradients[a]);
            EXPECT_LT(std::abs(blocks.state_coarse_hamiltonian[i + u.size() * j] - DV(fine) * k), 2e-10);
            EXPECT_LT(std::abs(blocks.state_coarse_overlap[i + u.size() * j] - DV(fine) * Dot(u[i].values, e[j].values)), 2e-11);
        }
        for (int i = 0; i < coarse.size(); ++i)
        {
            Complex k = Dot(e[i].values, local) + Dot(e[i].values, nonlocal);
            for (int a = 0; a < 3; ++a) k += h->kinetic_prefactor() * Dot(e[i].gradients[a], e[j].gradients[a]);
            EXPECT_LT(std::abs(image[i] - DV(fine) * k), 2e-10);
        }
    }
}

TEST(SternheimerWeakGrid, FineVertexProjectionIsAdjointOfPhysicalLift)
{
    const auto coarse = Grid(2, 3, 2), fine = Grid(5, 5, 4);
    Operator op(MakeH(fine), coarse);
    Vector x(coarse.size()), g(fine.size()), lifted, projected;
    for (int i = 0; i < coarse.size(); ++i) x[i] = Complex(std::sin(0.4 * i), std::cos(0.3 * i));
    for (int i = 0; i < fine.size(); ++i) g[i] = Complex(std::cos(0.17 * i), std::sin(0.31 * i));
    op.lift(x, lifted);
    op.project(g, projected);
    EXPECT_LT(std::abs(Dot(x, projected) - DV(fine) * Dot(lifted, g)), 2e-11);
    op.project(lifted, projected);
    Near(x, projected);
    Vector expected;
    op.apply(x, expected);
    op.apply(x, x);
    Near(x, expected);
}

TEST(SternheimerWeakGrid, RejectsInvalidMetricAndInputsWithoutChangingOutput)
{
    const auto coarse = Grid(2, 3, 2), fine = Grid(5, 5, 4);
    const auto h = MakeH(fine);
    EXPECT_THROW(Operator(nullptr, coarse), std::invalid_argument);
    EXPECT_THROW(Operator(h, coarse, 1), std::length_error);
    Operator op(h, coarse);
    auto u = std::vector<Function>{PlaneWave(fine, {0, 0, 0})};
    EXPECT_THROW(op.assemble_blocks(u, 1e-8, 1), std::length_error);
    u.push_back(u.front());
    EXPECT_THROW(op.assemble_blocks(u), std::invalid_argument);
    u.pop_back();
    u[0].gradients[0].clear();
    EXPECT_THROW(op.assemble_blocks(u), std::invalid_argument);
    Vector x(coarse.size()), output{Complex(17.0)};
    x[0] = std::numeric_limits<double>::infinity();
    EXPECT_THROW(op.apply(x, output), std::invalid_argument);
    ASSERT_EQ(output.size(), 1u);
    EXPECT_EQ(output[0], Complex(17.0));
}

TEST(SternheimerWeakGrid, RejectsNonHermitianOrNonfiniteNonlocalData)
{
    const auto coarse = Grid(2, 3, 2), fine = Grid(5, 5, 4);
    using P = ModuleRI::SternheimerFDNonlocalProjector;
    P::ProjectorBlock block;
    block.projectors = {PlaneWave(fine, {0, 0, 0}).values};
    block.d_matrix = {{Complex(0.4, 0.2)}};
    const auto make = [&]() {
        auto p = std::make_shared<const P>(fine.size(), DV(fine), std::vector<P::ProjectorBlock>{block});
        return std::make_shared<const H>(fine, std::vector<double>(fine.size(), 0.0), 0.7, p, 8);
    };
    EXPECT_THROW(Operator(make(), coarse), std::invalid_argument);
    block.d_matrix[0][0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(Operator(make(), coarse), std::invalid_argument);
    block.d_matrix[0][0] = 0.4;
    block.projectors[0][0] = std::numeric_limits<double>::infinity();
    EXPECT_THROW(Operator(make(), coarse), std::invalid_argument);
}

using Functions = std::vector<Function>;

void ExpectSameFunctions(const Functions& actual, const Functions& expected)
{
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        EXPECT_EQ(actual[i].values, expected[i].values);
        for (int a = 0; a < 3; ++a) EXPECT_EQ(actual[i].gradients[a], expected[i].gradients[a]);
    }
}

Functions ReferenceStates(const H::Grid& fine)
{
    return {PlaneWave(fine, {0, 0, 0}), PlaneWave(fine, {6, 0, 0}), PlaneWave(fine, {2, 1, -1})};
}

// U=Q*R with known upper R and positive diagonal; the requested Cholesky right
// solve must recover Q, including its supplied (not FFT-reconstructed) gradients.
Functions MixedStates(const Functions& reference)
{
    const std::array<std::array<Complex, 3>, 3> r{{
        {{Complex(1.8), Complex(0.3, 0.2), Complex(-0.15, 0.4)}},
        {{Complex(0.0), Complex(1.4), Complex(-0.2, -0.3)}},
        {{Complex(0.0), Complex(0.0), Complex(0.9)}}}};
    Functions result = reference;
    for (std::size_t j = 0; j < result.size(); ++j)
    {
        std::fill(result[j].values.begin(), result[j].values.end(), Complex(0.0));
        for (auto& gradient : result[j].gradients) std::fill(gradient.begin(), gradient.end(), Complex(0.0));
        for (std::size_t i = 0; i < reference.size(); ++i)
            for (std::size_t k = 0; k < result[j].values.size(); ++k)
            {
                result[j].values[k] += reference[i].values[k] * r[i][j];
                for (int a = 0; a < 3; ++a)
                    result[j].gradients[a][k] += reference[i].gradients[a][k] * r[i][j];
            }
    }
    return result;
}

TEST(SternheimerWeakGrid, InPlaceCholeskyNormalizationPreservesPrefixesAndAnalyticGradients)
{
    const auto fine = Grid(5, 5, 4);
    const Functions reference = ReferenceStates(fine);
    Functions states = MixedStates(reference);
    std::vector<std::array<const Complex*, 4>> addresses;
    for (const auto& state : states)
        addresses.push_back({{state.values.data(), state.gradients[0].data(),
                              state.gradients[1].data(), state.gradients[2].data()}});
    // Far less than a full fine-U copy; forces several tiles and a partial tail.
    const double error = ModuleRI::orthonormalize_sternheimer_weak_states_in_place(states, DV(fine), 1e-12, 1024);
    EXPECT_LT(error, 2e-12);
    ASSERT_EQ(states.size(), reference.size());
    for (std::size_t i = 0; i < states.size(); ++i)
    {
        EXPECT_EQ(states[i].values.data(), addresses[i][0]);
        Near(states[i].values, reference[i].values, 2e-12);
        for (int a = 0; a < 3; ++a)
        {
            EXPECT_EQ(states[i].gradients[a].data(), addresses[i][a + 1]);
            Near(states[i].gradients[a], reference[i].gradients[a], 2e-11);
        }
        for (std::size_t j = 0; j < states.size(); ++j)
            EXPECT_LT(std::abs(DV(fine) * Dot(states[i].values, states[j].values)
                               - Complex(i == j ? 1.0 : 0.0)), 2e-12);
    }
    for (std::size_t prefix = 1; prefix <= states.size(); ++prefix)
    {
        Functions independent(reference.begin(), reference.begin() + prefix);
        const Functions mixed = MixedStates(reference);
        independent.assign(mixed.begin(), mixed.begin() + prefix);
        ModuleRI::orthonormalize_sternheimer_weak_states_in_place(independent, DV(fine));
        for (std::size_t i = 0; i < prefix; ++i)
        {
            Near(independent[i].values, states[i].values, 2e-12);
            for (int a = 0; a < 3; ++a) Near(independent[i].gradients[a], states[i].gradients[a], 2e-11);
        }
    }
}

TEST(SternheimerWeakGrid, NormalizationIsIndependentOfTileBudgetAndHandlesOneState)
{
    const auto fine = Grid(5, 5, 4);
    Functions tiled = MixedStates(ReferenceStates(fine));
    Functions large_tile = tiled;
    ModuleRI::orthonormalize_sternheimer_weak_states_in_place(tiled, DV(fine), 1e-12, 1024);
    ModuleRI::orthonormalize_sternheimer_weak_states_in_place(large_tile, DV(fine));
    for (std::size_t i = 0; i < tiled.size(); ++i)
    {
        Near(tiled[i].values, large_tile[i].values, 2e-12);
        for (int a = 0; a < 3; ++a) Near(tiled[i].gradients[a], large_tile[i].gradients[a], 2e-11);
    }
    Functions single{PlaneWave(fine, {6, 0, 0})};
    const Function original = single[0];
    for (auto& value : single[0].values) value *= Complex(0.0, 2.0);
    for (auto& gradient : single[0].gradients)
        for (auto& value : gradient) value *= Complex(0.0, 2.0);
    EXPECT_LT(ModuleRI::orthonormalize_sternheimer_weak_states_in_place(single, DV(fine)), 2e-12);
    Vector expected = original.values;
    for (auto& value : expected) value *= Complex(0.0, 1.0);
    Near(single[0].values, expected, 2e-12);
}

TEST(SternheimerWeakGrid, NormalizationRejectsSingularAndIllConditionedMetricWithoutMutation)
{
    const auto fine = Grid(5, 5, 4);
    const Functions reference = ReferenceStates(fine);
    Functions states{reference[0], reference[0]};
    Functions original = states;
    EXPECT_THROW(ModuleRI::orthonormalize_sternheimer_weak_states_in_place(states, DV(fine)), std::domain_error);
    ExpectSameFunctions(states, original);
    for (std::size_t k = 0; k < states[1].values.size(); ++k)
    {
        states[1].values[k] += 1e-7 * reference[1].values[k];
        for (int a = 0; a < 3; ++a) states[1].gradients[a][k] += 1e-7 * reference[1].gradients[a][k];
    }
    original = states;
    EXPECT_THROW(ModuleRI::orthonormalize_sternheimer_weak_states_in_place(states, DV(fine)), std::domain_error);
    ExpectSameFunctions(states, original);
    states = reference;
    for (auto& value : states[0].values) value *= 1e8;
    for (auto& gradient : states[0].gradients)
        for (auto& value : gradient) value *= 1e8;
    original = states;
    EXPECT_THROW(ModuleRI::orthonormalize_sternheimer_weak_states_in_place(states, DV(fine)), std::domain_error);
    ExpectSameFunctions(states, original);
}

TEST(SternheimerWeakGrid, NormalizationValidatesArgumentsBeforeMutation)
{
    const auto fine = Grid(5, 5, 4);
    Functions states = MixedStates(ReferenceStates(fine));
    const Functions original = states;
    EXPECT_THROW(ModuleRI::orthonormalize_sternheimer_weak_states_in_place(states, DV(fine), 1e-12, 1),
                 std::length_error);
    for (double dv : {0.0, -0.1, std::numeric_limits<double>::infinity(),
                      std::numeric_limits<double>::quiet_NaN()})
        EXPECT_THROW(ModuleRI::orthonormalize_sternheimer_weak_states_in_place(states, dv), std::invalid_argument);
    for (double threshold : {0.0, -0.1, 1.0, std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN()})
        EXPECT_THROW(ModuleRI::orthonormalize_sternheimer_weak_states_in_place(states, DV(fine), threshold),
                     std::invalid_argument);
    ExpectSameFunctions(states, original);
    states.back().gradients[2].pop_back();
    const Functions malformed = states;
    EXPECT_THROW(ModuleRI::orthonormalize_sternheimer_weak_states_in_place(states, DV(fine)), std::invalid_argument);
    ExpectSameFunctions(states, malformed);
    states = original;
    states.back().gradients[2].back() = std::numeric_limits<double>::infinity();
    const Functions nonfinite = states;
    EXPECT_THROW(ModuleRI::orthonormalize_sternheimer_weak_states_in_place(states, DV(fine)), std::invalid_argument);
    ExpectSameFunctions(states, nonfinite);
    Functions empty;
    EXPECT_THROW(ModuleRI::orthonormalize_sternheimer_weak_states_in_place(empty, DV(fine)), std::invalid_argument);
}

TEST(SternheimerWeakGrid, NormalizationReportsArithmeticOverflow)
{
    Function state;
    state.values = {0.01};
    for (auto& gradient : state.gradients) gradient = {std::numeric_limits<double>::max()};
    Functions states{state};
    EXPECT_THROW(ModuleRI::orthonormalize_sternheimer_weak_states_in_place(states, 1.0), std::overflow_error);
}

#ifdef ABACUS_STERNHEIMER_WEAK_EXACT_APPLY_CACHE
std::shared_ptr<const H> MakeCacheH(const H::Grid& coarse, const H::Grid& fine, bool nonlocal = true)
{
    const std::array<int, 3> c{{coarse.nx, coarse.ny, coarse.nz}};
    const std::array<int, 3> f{{fine.nx, fine.ny, fine.nz}};
    std::vector<double> potential(fine.size(), 0.7);
    for (int i = 0; i < fine.size(); ++i)
    {
        const auto r = Position(fine, i);
        // Required difference-window endpoints, just-outside modes and fine Nyquist.
        for (int a = 0; a < 3; ++a)
            potential[i] += 0.17 * std::sin(pi2 * (c[a] - 1) * r[a])
                            + 0.23 * std::cos(pi2 * c[a] * r[a])
                            + 0.13 * std::cos(pi2 * (f[a] / 2) * r[a]);
        potential[i] += 0.19 * std::sin(pi2 * ((c[0] - 1) * r[0] - (c[1] - 1) * r[1]
                                              + (c[2] - 1) * r[2]));
    }
    using P = ModuleRI::SternheimerFDNonlocalProjector;
    std::shared_ptr<const P> projector;
    if (nonlocal)
    {
        P::ProjectorBlock first, second;
        first.projectors.assign(2, Vector(fine.size()));
        second.projectors.assign(1, Vector(fine.size()));
        for (int i = 0; i < fine.size(); ++i)
        {
            first.projectors[0][i] = Complex(std::sin(0.31 * i), std::cos(0.57 * i)) / std::sqrt(30.0);
            first.projectors[1][i] = Complex(std::cos(0.71 * i), std::sin(0.29 * i)) / std::sqrt(30.0);
            second.projectors[0][i] = i % 3 ? Complex(0.0) : Complex(0.11, -0.09);
        }
        first.d_matrix = {{Complex(0.4), Complex(0.08, 0.27)}, {Complex(0.08, -0.27), Complex(-0.3)}};
        second.d_matrix = {{Complex(-0.17)}};
        projector = std::make_shared<const P>(fine.size(), DV(fine), std::vector<P::ProjectorBlock>{first, second});
    }
    return std::make_shared<const H>(fine, potential, 0.7, projector, 8);
}

void CheckCachedActions(const H::Grid& coarse, const H::Grid& fine, bool nonlocal, bool local)
{
    const auto h = MakeCacheH(coarse, fine);
    Operator original(h, coarse), cached(h, coarse);
    EXPECT_FALSE(cached.has_exact_nonlocal_cache());
    EXPECT_FALSE(cached.has_exact_local_cache());
    EXPECT_EQ(cached.exact_cache_storage_bytes(), 0u);
    cached.enable_exact_apply_cache(nonlocal, local);
    EXPECT_EQ(cached.has_exact_nonlocal_cache(), nonlocal);
    EXPECT_EQ(cached.has_exact_local_cache(), local);
    Vector x(coarse.size()), y(coarse.size()), expected, actual, hy;
    for (int i = 0; i < coarse.size(); ++i)
    {
        x[i] = Complex(std::sin(0.37 * i + 0.2), std::cos(0.53 * i));
        y[i] = Complex(std::cos(0.61 * i), std::sin(0.19 * i + 0.3));
    }
    original.apply(x, expected);
    cached.apply(x, actual);
    Near(actual, expected);
    cached.apply(y, hy);
    EXPECT_LT(std::abs(Dot(y, actual) - Dot(hy, x)), 2e-9);
    cached.apply(x, x);
    Near(x, expected);
    // Visit every mode, including the asymmetric positive coarse Nyquist.
    const std::array<int, 3> dims{{coarse.nx, coarse.ny, coarse.nz}};
    for (int i = 0; i < coarse.size(); ++i)
    {
        std::array<int, 3> mode{{i / (coarse.ny * coarse.nz), (i / coarse.nz) % coarse.ny, i % coarse.nz}};
        for (int a = 0; a < 3; ++a) if (mode[a] > dims[a] / 2) mode[a] -= dims[a];
        x = PlaneWave(coarse, mode).values;
        original.apply(x, expected);
        cached.apply(x, actual);
        Near(actual, expected);
    }
}

TEST(SternheimerWeakGrid, ExactNonlocalCacheMatchesFineApplyAtSkewNonzeroK)
{
    CheckCachedActions(Grid(2, 3, 2), Grid(7, 7, 6), true, false);
    CheckCachedActions(Grid(3, 2, 3), Grid(8, 5, 8), true, false);
}

TEST(SternheimerWeakGrid, ExactLocalCacheRetainsBoundaryDifferencesAndRejectsAliasedHighModes)
{
    CheckCachedActions(Grid(2, 3, 2), Grid(7, 9, 8), false, true);
    CheckCachedActions(Grid(3, 2, 3), Grid(9, 7, 9), false, true);
}

TEST(SternheimerWeakGrid, ExactCombinedCacheHandlesFineLimitedAxesAndIdentityTransfer)
{
    CheckCachedActions(Grid(3, 4, 2), Grid(4, 7, 3), true, true);
    CheckCachedActions(Grid(3, 2, 2), Grid(3, 2, 2), true, true);
    CheckCachedActions(Grid(1, 2, 1), Grid(1, 5, 4), true, true);
}

TEST(SternheimerWeakGrid, ExactCacheSupportsLargeExplicitBlochLabelsAnd144Axis)
{
    auto coarse = Grid(3, 2, 2), fine = Grid(144, 5, 6);
    coarse.kpoint = fine.kpoint = {3.17, -2.13, 1.21};
    CheckCachedActions(coarse, fine, true, true);
}

TEST(SternheimerWeakGrid, ExactCacheLeavesFineBlocksLiftAndProjectUnchanged)
{
    const auto coarse = Grid(2, 3, 2), fine = Grid(5, 7, 6);
    Operator op(MakeCacheH(coarse, fine), coarse);
    const Functions states{PlaneWave(fine, {0, 0, 0}), PlaneWave(fine, {6, 0, 0})};
    const auto before = op.assemble_blocks(states);
    Vector x(coarse.size(), Complex(0.2, -0.1)), lifted, projected, after_lift, after_project;
    op.lift(x, lifted);
    op.project(lifted, projected);
    op.enable_exact_apply_cache(true, true);
    const auto after = op.assemble_blocks(states);
    Near(before.state_hamiltonian, after.state_hamiltonian);
    Near(before.state_coarse_overlap, after.state_coarse_overlap);
    Near(before.state_coarse_hamiltonian, after.state_coarse_hamiltonian);
    op.lift(x, after_lift);
    op.project(lifted, after_project);
    Near(lifted, after_lift);
    Near(projected, after_project);
}

TEST(SternheimerWeakGrid, ExactCacheBudgetFailuresAreTransactionalAndCachesCanBeDisabled)
{
    const auto coarse = Grid(2, 3, 2), fine = Grid(7, 9, 8);
    const auto h = MakeCacheH(coarse, fine);
    Operator op(h, coarse);
    const auto budget = Operator::exact_cache_workspace_bytes_required(*h, coarse, true, true);
    EXPECT_GT(budget, 0u);
    EXPECT_THROW(op.enable_exact_apply_cache(true, true, budget - 1), std::length_error);
    EXPECT_EQ(op.exact_cache_storage_bytes(), 0u);
    op.enable_exact_apply_cache(true, true, budget);
    EXPECT_GT(op.exact_cache_storage_bytes(), 0u);
    EXPECT_LE(op.exact_cache_storage_bytes(), budget);
    Vector x(coarse.size(), Complex(0.2, 0.7)), expected, actual;
    op.apply(x, expected);
    EXPECT_THROW(op.enable_exact_apply_cache(true, false, 0), std::length_error);
    EXPECT_TRUE(op.has_exact_local_cache());
    EXPECT_TRUE(op.has_exact_nonlocal_cache());
    op.apply(x, actual);
    Near(actual, expected);
    Vector sentinel{Complex(17.0, -2.0)};
    EXPECT_THROW(op.apply(Vector(1), sentinel), std::invalid_argument);
    x[0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(op.apply(x, sentinel), std::invalid_argument);
    ASSERT_EQ(sentinel.size(), 1u);
    EXPECT_EQ(sentinel[0], Complex(17.0, -2.0));
    op.enable_exact_apply_cache(false, false, 0);
    EXPECT_EQ(op.exact_cache_storage_bytes(), 0u);
    EXPECT_FALSE(op.has_exact_nonlocal_cache());
    EXPECT_FALSE(op.has_exact_local_cache());
}

TEST(SternheimerWeakGrid, ExactCacheHandlesAbsentProjectorsAndValidatesGridBeforeAllocation)
{
    const auto coarse = Grid(3, 2, 2), fine = Grid(7, 5, 6);
    const auto h = MakeCacheH(coarse, fine, false);
    Operator original(h, coarse), cached(h, coarse);
    cached.enable_exact_apply_cache(true, true);
    Vector x(coarse.size(), Complex(0.7, 0.3)), expected, actual;
    original.apply(x, expected);
    cached.apply(x, actual);
    Near(expected, actual);
    auto invalid = coarse;
    invalid.nx = 0;
    EXPECT_THROW(Operator::exact_cache_workspace_bytes_required(*h, invalid, true, true), std::invalid_argument);
    invalid = coarse;
    invalid.kpoint[0] += 1;
    EXPECT_THROW(Operator::exact_cache_workspace_bytes_required(*h, invalid, true, true), std::invalid_argument);
    EXPECT_EQ(Operator::exact_cache_workspace_bytes_required(*h, coarse, false, false), 0u);
}
#else
TEST(SternheimerWeakGrid, ExactApplyCacheRequiresImplementation)
{
    FAIL() << "Missing opt-in exact apply cache: expected test-first RED before implementation";
}
#endif
} // namespace
