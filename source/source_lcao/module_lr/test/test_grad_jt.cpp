#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "../grad_degen.h"

// ----------------------------- the Jahn-Teller direction search -----------------------------

namespace
{
    /// A fixed, deliberately unsymmetric G tensor: `ncoord` blocks of d x d, symmetric in (k, l).
    std::vector<double> sample_gflat(const int ncoord, const int d, const unsigned seed)
    {
        std::vector<double> g(static_cast<size_t>(ncoord) * d * d, 0.0);
        unsigned x = seed;
        const auto next = [&x]() {
            x = x * 1664525u + 1013904223u;   // deterministic; a test must not depend on rand()
            return static_cast<double>(static_cast<int>(x >> 8) % 2000 - 1000) / 500.0;
        };
        for (int a = 0; a < ncoord; ++a)
        {
            double* const b = g.data() + static_cast<size_t>(a) * d * d;
            for (int k = 0; k < d; ++k)
            {
                for (int l = k; l < d; ++l)
                {
                    const double v = next();
                    b[k * d + l] = v;
                    b[l * d + k] = v;
                }
            }
        }
        return g;
    }

    double branch_grad_norm(const std::vector<double>& g, const int ncoord, const int d,
        const std::vector<double>& v)
    {
        double s = 0.0;
        for (int a = 0; a < ncoord; ++a)
        {
            const double* const b = g.data() + static_cast<size_t>(a) * d * d;
            double q = 0.0;
            for (int k = 0; k < d; ++k)
            {
                for (int l = 0; l < d; ++l) { q += v[k] * b[k * d + l] * v[l]; }
            }
            s += q * q;
        }
        return std::sqrt(s);
    }

    /// lambda_min of sum_a u_a G^(a), by brute force over the 2x2 / 3x3 block.
    double lambda_min_along(const std::vector<double>& g, const int ncoord, const int d,
        const std::vector<double>& u)
    {
        std::vector<double> m(static_cast<size_t>(d) * d, 0.0);
        for (int a = 0; a < ncoord; ++a)
        {
            const double* const b = g.data() + static_cast<size_t>(a) * d * d;
            for (int i = 0; i < d * d; ++i) { m[i] += u[a] * b[i]; }
        }
        // smallest Rayleigh quotient, sampled densely; d is 2 or 3 in these tests
        double lo = 1e300;
        const int n = 2000;
        if (d == 2)
        {
            for (int i = 0; i <= n; ++i)
            {
                const double t = M_PI * i / n;
                const double v[2] = { std::cos(t), std::sin(t) };
                const double r = v[0] * v[0] * m[0] + 2.0 * v[0] * v[1] * m[1] + v[1] * v[1] * m[3];
                lo = std::min(lo, r);
            }
        }
        else
        {
            for (int i = 0; i <= 200; ++i)
            {
                for (int j = 0; j <= 400; ++j)
                {
                    const double th = M_PI * i / 200;
                    const double ph = 2.0 * M_PI * j / 400;
                    const double v[3] = { std::sin(th) * std::cos(ph), std::sin(th) * std::sin(ph),
                                          std::cos(th) };
                    double r = 0.0;
                    for (int p = 0; p < 3; ++p)
                    {
                        for (int q = 0; q < 3; ++q) { r += v[p] * m[p * 3 + q] * v[q]; }
                    }
                    lo = std::min(lo, r);
                }
            }
        }
        return lo;
    }
}

/// The identity the whole search rests on:
///     min_{|u|=1} lambda_min(M(u)) = -max_{|v|=1} |q(v)|.
/// Checked against a brute-force scan of the left-hand side over directions built from the
/// returned optimum, so a wrong reduction cannot pass.
TEST(JTDirection, ReductionIdentityHolds)
{
    for (const int d : { 2, 3 })
    {
        const int ncoord = 6;
        const std::vector<double> g = sample_gflat(ncoord, d, 7u + d);
        const LR::JTDirection jt = LR::find_jt_direction(g, ncoord, d);
        ASSERT_EQ(static_cast<int>(jt.displacement.size()), ncoord) << "d=" << d;
        ASSERT_EQ(static_cast<int>(jt.mixing.size()), d);
        EXPECT_GT(jt.slope, 0.0);

        // |q(v*)| must equal the reported slope
        EXPECT_NEAR(branch_grad_norm(g, ncoord, d, jt.mixing), jt.slope, 1e-9 * jt.slope);
        // the displacement is the normalized branch gradient
        double n = 0.0;
        for (int a = 0; a < ncoord; ++a) { n += jt.displacement[a] * jt.displacement[a]; }
        EXPECT_NEAR(std::sqrt(n), 1.0, 1e-12);
        // and along MINUS it the smallest eigenvalue is -slope: the two sides of the identity
        std::vector<double> u(ncoord);
        for (int a = 0; a < ncoord; ++a) { u[a] = -jt.displacement[a]; }
        EXPECT_NEAR(lambda_min_along(g, ncoord, d, u), -jt.slope, 1e-4 * jt.slope) << "d=" << d;
    }
}

/// Optimality: no other unit mixing may give a larger branch-gradient norm. Brute-forced over the
/// subspace sphere, which is the independent check that the alternating iteration converged to the
/// global optimum and not merely to a stationary point.
TEST(JTDirection, IsGlobalOverTheSubspace)
{
    const int ncoord = 6;
    {
        const int d = 2;
        const std::vector<double> g = sample_gflat(ncoord, d, 9u);
        const LR::JTDirection jt = LR::find_jt_direction(g, ncoord, d);
        double best = 0.0;
        for (int i = 0; i <= 4000; ++i)
        {
            const double t = M_PI * i / 4000;
            best = std::max(best, branch_grad_norm(g, ncoord, d, { std::cos(t), std::sin(t) }));
        }
        EXPECT_NEAR(jt.slope, best, 1e-6 * best);
    }
    {
        const int d = 3;
        const std::vector<double> g = sample_gflat(ncoord, d, 11u);
        const LR::JTDirection jt = LR::find_jt_direction(g, ncoord, d);
        double best = 0.0;
        for (int i = 0; i <= 300; ++i)
        {
            for (int j = 0; j <= 600; ++j)
            {
                const double th = M_PI * i / 300;
                const double ph = 2.0 * M_PI * j / 600;
                best = std::max(best, branch_grad_norm(g, ncoord, d,
                    { std::sin(th) * std::cos(ph), std::sin(th) * std::sin(ph), std::cos(th) }));
            }
        }
        EXPECT_NEAR(jt.slope, best, 1e-4 * best);
    }
}

/// A multiplet whose every G block is a multiple of the identity has no Jahn-Teller direction to
/// find: all branches share one gradient, so |q(v)| is the same for every v and nothing splits.
/// The search must still return that common direction rather than something arbitrary.
TEST(JTDirection, DegenerateCaseGivesTheCommonGradient)
{
    const int ncoord = 6;
    const int d = 2;
    std::vector<double> g(static_cast<size_t>(ncoord) * d * d, 0.0);
    const double comm[6] = { 0.3, -0.7, 1.1, 0.0, 0.5, -0.2 };
    for (int a = 0; a < ncoord; ++a)
    {
        g[static_cast<size_t>(a) * 4 + 0] = comm[a];
        g[static_cast<size_t>(a) * 4 + 3] = comm[a];
    }
    const LR::JTDirection jt = LR::find_jt_direction(g, ncoord, d);
    double n = 0.0;
    for (int a = 0; a < ncoord; ++a) { n += comm[a] * comm[a]; }
    n = std::sqrt(n);
    EXPECT_NEAR(jt.slope, n, 1e-10);
    for (int a = 0; a < ncoord; ++a) { EXPECT_NEAR(jt.displacement[a], comm[a] / n, 1e-10); }
    // every start reaches the same value, since the objective is constant on the sphere
    EXPECT_GE(jt.restarts_agreeing, 2);
}

/// Homogeneity: scaling G scales the slope and leaves the direction alone.
TEST(JTDirection, ScalesLinearly)
{
    const int ncoord = 6;
    const int d = 3;
    const std::vector<double> g = sample_gflat(ncoord, d, 13u);
    std::vector<double> g3 = g;
    for (size_t i = 0; i < g3.size(); ++i) { g3[i] *= 3.0; }
    const LR::JTDirection a = LR::find_jt_direction(g, ncoord, d);
    const LR::JTDirection b = LR::find_jt_direction(g3, ncoord, d);
    EXPECT_NEAR(b.slope, 3.0 * a.slope, 1e-8 * a.slope);
    for (int i = 0; i < ncoord; ++i) { EXPECT_NEAR(b.displacement[i], a.displacement[i], 1e-8); }
}

/// The symmetric / Jahn-Teller split must reconstruct the branch gradient, and the symmetric part
/// must be the multiplet average -- independent of which branch was chosen.
TEST(JTDirection, SymmetricSplitReconstructs)
{
    const int ncoord = 6;
    const int d = 3;
    const std::vector<double> g = sample_gflat(ncoord, d, 17u);
    const LR::JTDirection jt = LR::find_jt_direction(g, ncoord, d);
    std::vector<double> jt_part;
    const std::vector<double> sym = LR::split_symmetric_part(g, ncoord, d, jt.mixing, jt_part);
    for (int a = 0; a < ncoord; ++a)
    {
        const double* const b = g.data() + static_cast<size_t>(a) * d * d;
        double q = 0.0;
        for (int k = 0; k < d; ++k)
        {
            for (int l = 0; l < d; ++l) { q += jt.mixing[k] * b[k * d + l] * jt.mixing[l]; }
        }
        EXPECT_NEAR(sym[a] + jt_part[a], q, 1e-12);
        double tr = 0.0;
        for (int k = 0; k < d; ++k) { tr += b[k * d + k]; }
        EXPECT_NEAR(sym[a], tr / d, 1e-12);
    }
    // a different mixing keeps the same symmetric part
    std::vector<double> other(d, 0.0);
    other[0] = 1.0;
    std::vector<double> jt2;
    const std::vector<double> sym2 = LR::split_symmetric_part(g, ncoord, d, other, jt2);
    for (int a = 0; a < ncoord; ++a) { EXPECT_NEAR(sym2[a], sym[a], 1e-14); }
}

TEST(JTDirection, ZeroGradientReturnsAValidStationaryBranch)
{
    for (const int d : {2, 3})
    {
        const int ncoord = 3;
        const size_t size = static_cast<size_t>(ncoord) * d * d;
        const std::vector<double> zero(size, 0.0);
        const LR::JTDirection jt = LR::find_jt_direction(zero, ncoord, d);
        ASSERT_EQ(jt.mixing.size(), d);
        ASSERT_EQ(jt.displacement.size(), ncoord);
        EXPECT_DOUBLE_EQ(jt.slope, 0.0);
        EXPECT_DOUBLE_EQ(jt.mixing[0], 1.0);
        double norm_squared = 0.0;
        for (const double value : jt.mixing) { norm_squared += value * value; }
        EXPECT_DOUBLE_EQ(norm_squared, 1.0);
        std::vector<double> jt_part;
        const std::vector<double> sym = LR::split_symmetric_part(zero, ncoord, d, jt.mixing, jt_part);
        EXPECT_EQ(sym, (std::vector<double>(ncoord, 0.0)));
        EXPECT_EQ(jt_part, sym);
        EXPECT_EQ(jt.displacement, sym);
    }
}
