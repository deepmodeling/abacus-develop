#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "source_base/matrix.h"
#include "../grad_matrix_degenerate.h"

/// Tests for the degenerate-subspace gradient matrix algebra; `../grad_matrix_degenerate.h`
/// states the identity being exercised.
///
/// The interesting test is not the arithmetic of `assemble_grad_matrix` but whether the route it
/// implements really recovers the bilinear form: `QuadraticForm` below plays the part of the
/// analytic-gradient pipeline (a genuine quadratic form in X), and the tests check that feeding it
/// the normalized combinations reproduces $G_{kl}=X_k^\top B X_l$ and that the result is covariant
/// under a rotation of the basis of the degenerate subspace -- the two self-checks the document
/// lists in section 5.4.5 as the ones needing no finite differences.

namespace
{
    constexpr int nov = 5;   ///< ambient particle-hole space, stands in for nocc*nvirt

    /// $\mathcal F[X]=X^\top B X$ with B symmetric: the same shape the real force map has on a
    /// multiplet. Scalar-valued, which is enough -- the real one is just 3N independent copies.
    struct QuadraticForm
    {
        std::vector<double> b;   ///< nov*nov, symmetric

        QuadraticForm() : b(nov * nov, 0.0)
        {
            // arbitrary but fixed and definitely not diagonal, so the off-diagonal elements of G
            // are nonzero and the test can fail
            const double raw[nov][nov] = { { 1.3, -0.7, 0.4, 0.9, -0.2 },
                                           { 0.0, 2.1, -1.1, 0.3, 0.6 },
                                           { 0.0, 0.0, -0.5, 0.8, 1.4 },
                                           { 0.0, 0.0, 0.0, 0.7, -0.9 },
                                           { 0.0, 0.0, 0.0, 0.0, 1.9 } };
            for (int i = 0; i < nov; ++i)
            {
                for (int j = i; j < nov; ++j)
                {
                    b[i * nov + j] = raw[i][j];
                    b[j * nov + i] = raw[i][j];
                }
            }
        }

        /// the bilinear form B behind it, i.e. the exact answer
        double bilinear(const std::vector<double>& x, const std::vector<double>& y) const
        {
            double s = 0.0;
            for (int i = 0; i < nov; ++i)
            {
                for (int j = 0; j < nov; ++j) { s += x[i] * b[i * nov + j] * y[j]; }
            }
            return s;
        }

        /// what "the code" returns for one state
        double operator()(const std::vector<double>& x) const { return bilinear(x, x); }
    };

    /// three orthonormal vectors spanning the stand-in degenerate subspace
    std::vector<std::vector<double>> subspace_basis()
    {
        std::vector<std::vector<double>> x;
        x.push_back({ 1.0, 0.0, 0.0, 0.0, 0.0 });
        x.push_back({ 0.0, 0.6, 0.0, 0.8, 0.0 });
        x.push_back({ 0.0, 0.8, 0.0, -0.6, 0.0 });
        return x;
    }

    /// run route A2 over a subspace basis and return G
    std::vector<std::vector<double>> route_a2(const QuadraticForm& f,
        const std::vector<std::vector<double>>& x)
    {
        const int d = static_cast<int>(x.size());
        const std::vector<std::pair<int, int>> pairs = LR::degenerate_pairs(d);
        std::vector<double> diag(d);
        for (int k = 0; k < d; ++k) { diag[k] = f(x[k]); }
        std::vector<double> plus(pairs.size());
        for (size_t ip = 0; ip < pairs.size(); ++ip)
        {
            std::vector<double> xp(nov);
            LR::combine_normalized(x[pairs[ip].first].data(), x[pairs[ip].second].data(), nov,
                xp.data());
            plus[ip] = f(xp);
        }
        return LR::assemble_grad_matrix(diag, plus, pairs);
    }
}

/// The whole point: route A2 reproduces the bilinear form exactly, off-diagonal included.
TEST(GradMatrixDegenerate, PolarizationRecoversBilinearForm)
{
    const QuadraticForm f;
    const std::vector<std::vector<double>> x = subspace_basis();
    const std::vector<std::vector<double>> g = route_a2(f, x);

    const int d = static_cast<int>(x.size());
    bool any_offdiag = false;
    for (int k = 0; k < d; ++k)
    {
        for (int l = 0; l < d; ++l)
        {
            EXPECT_NEAR(g[k][l], f.bilinear(x[k], x[l]), 1e-12) << "k=" << k << " l=" << l;
            if (k != l && std::abs(g[k][l]) > 1e-6) { any_offdiag = true; }
        }
    }
    // guard against a degenerate test case that would pass with G assembled as zero off-diagonal
    EXPECT_TRUE(any_offdiag);
}

/// G is symmetric by construction; assert it so a future refactor cannot lose it silently.
TEST(GradMatrixDegenerate, IsSymmetric)
{
    const QuadraticForm f;
    const std::vector<std::vector<double>> g = route_a2(f, subspace_basis());
    for (size_t k = 0; k < g.size(); ++k)
    {
        for (size_t l = 0; l < g.size(); ++l) { EXPECT_DOUBLE_EQ(g[k][l], g[l][k]); }
    }
}

/// Section 5.4.5(ii), the strongest check: rotating the basis of the degenerate subspace must
/// rotate G, $G'=U^\top G U$. This is what would fail if the pipeline were not a pure quadratic
/// form, and it needs no finite differences.
TEST(GradMatrixDegenerate, RotationCovariance)
{
    const QuadraticForm f;
    const std::vector<std::vector<double>> x = subspace_basis();
    const std::vector<std::vector<double>> g = route_a2(f, x);
    const int d = static_cast<int>(x.size());

    // an orthogonal U mixing all three members (rotation by theta in (0,1), then by phi in (1,2))
    const double ct = std::cos(0.7);
    const double st = std::sin(0.7);
    const double cp = std::cos(0.4);
    const double sp = std::sin(0.4);
    std::vector<std::vector<double>> u(d, std::vector<double>(d, 0.0));
    u[0][0] = ct;
    u[0][1] = -st * cp;
    u[0][2] = st * sp;
    u[1][0] = st;
    u[1][1] = ct * cp;
    u[1][2] = -ct * sp;
    u[2][0] = 0.0;
    u[2][1] = sp;
    u[2][2] = cp;

    // rotated basis $X'_k=\sum_l U_{lk}X_l$
    std::vector<std::vector<double>> xr(d, std::vector<double>(nov, 0.0));
    for (int k = 0; k < d; ++k)
    {
        for (int l = 0; l < d; ++l)
        {
            for (int i = 0; i < nov; ++i) { xr[k][i] += u[l][k] * x[l][i]; }
        }
    }
    const std::vector<std::vector<double>> gr = route_a2(f, xr);

    for (int k = 0; k < d; ++k)
    {
        for (int l = 0; l < d; ++l)
        {
            double expect = 0.0;
            for (int p = 0; p < d; ++p)
            {
                for (int q = 0; q < d; ++q) { expect += u[p][k] * g[p][q] * u[q][l]; }
            }
            EXPECT_NEAR(gr[k][l], expect, 1e-12) << "k=" << k << " l=" << l;
        }
    }
}

/// The trace is basis-invariant while the individual diagonal elements are not -- the reason the
/// validation tables compare group means (section 1.1).
TEST(GradMatrixDegenerate, TraceIsBasisInvariantButDiagonalIsNot)
{
    const QuadraticForm f;
    const std::vector<std::vector<double>> x = subspace_basis();
    // swap-and-mix the last two members only
    std::vector<std::vector<double>> xr = x;
    const double r = 1.0 / std::sqrt(2.0);
    for (int i = 0; i < nov; ++i)
    {
        xr[1][i] = r * (x[1][i] + x[2][i]);
        xr[2][i] = r * (x[1][i] - x[2][i]);
    }
    const std::vector<std::vector<double>> g = route_a2(f, x);
    const std::vector<std::vector<double>> gr = route_a2(f, xr);

    double tr = 0.0;
    double trr = 0.0;
    for (size_t k = 0; k < g.size(); ++k)
    {
        tr += g[k][k];
        trr += gr[k][k];
    }
    EXPECT_NEAR(tr, trr, 1e-12);
    // and the mixing really did change the individual entries
    EXPECT_GT(std::abs(g[1][1] - gr[1][1]), 1e-6);
}

TEST(GradMatrixDegenerate, CombineNormalizedKeepsNorm)
{
    const std::vector<double> a = { 1.0, 0.0, 0.0, 0.0, 0.0 };
    const std::vector<double> b = { 0.0, 1.0, 0.0, 0.0, 0.0 };
    std::vector<double> p(nov);
    LR::combine_normalized(a.data(), b.data(), nov, p.data());
    double n = 0.0;
    for (int i = 0; i < nov; ++i) { n += p[i] * p[i]; }
    EXPECT_NEAR(n, 1.0, 1e-14);
    EXPECT_NEAR(p[0], 1.0 / std::sqrt(2.0), 1e-14);
    EXPECT_NEAR(p[1], 1.0 / std::sqrt(2.0), 1e-14);
}

TEST(GradMatrixDegenerate, PairOrder)
{
    const std::vector<std::pair<int, int>> p = LR::degenerate_pairs(3);
    ASSERT_EQ(p.size(), 3u);
    EXPECT_EQ(p[0], std::make_pair(0, 1));
    EXPECT_EQ(p[1], std::make_pair(0, 2));
    EXPECT_EQ(p[2], std::make_pair(1, 2));
    EXPECT_TRUE(LR::degenerate_pairs(1).empty());
    EXPECT_TRUE(LR::degenerate_pairs(0).empty());
    // d + d(d-1)/2 = d(d+1)/2 evaluations in total, the number of independent components
    for (int d = 1; d < 8; ++d)
    {
        EXPECT_EQ(static_cast<int>(LR::degenerate_pairs(d).size()) + d, d * (d + 1) / 2);
    }
}

TEST(GradMatrixDegenerate, GroupingCollectsMultiplets)
{
    // a triplet, then a singlet, then a doublet
    const std::vector<double> omega = { 0.5000000, 0.5000001, 0.5000002, 0.9, 1.30, 1.3000005 };
    const std::vector<std::vector<int>> g = LR::group_degenerate_states(omega, 2e-3);
    ASSERT_EQ(g.size(), 3u);
    EXPECT_EQ(g[0], std::vector<int>({ 0, 1, 2 }));
    EXPECT_EQ(g[1], std::vector<int>({ 3 }));
    EXPECT_EQ(g[2], std::vector<int>({ 4, 5 }));
}

TEST(GradMatrixDegenerate, GroupingSortsAndHandlesEdgeCases)
{
    // unsorted input must still group correctly, and the returned indices are into `omega`
    const std::vector<double> omega = { 1.3, 0.5, 1.3000005, 0.5000001 };
    const std::vector<std::vector<int>> g = LR::group_degenerate_states(omega, 2e-3);
    ASSERT_EQ(g.size(), 2u);
    EXPECT_EQ(g[0], std::vector<int>({ 1, 3 }));
    EXPECT_EQ(g[1], std::vector<int>({ 0, 2 }));

    // thr <= 0 disables grouping entirely (the default: current per-state behaviour)
    const std::vector<std::vector<int>> none = LR::group_degenerate_states(omega, 0.0);
    EXPECT_EQ(none.size(), omega.size());
    for (const std::vector<int>& grp : none) { EXPECT_EQ(grp.size(), 1u); }

    EXPECT_TRUE(LR::group_degenerate_states(std::vector<double>(), 2e-3).empty());
}

/// Anchoring to the group's first member, not the predecessor: a ladder of steps each below `thr`
/// must NOT chain into one group whose total spread exceeds it.
TEST(GradMatrixDegenerate, GroupingDoesNotChain)
{
    const double thr = 1e-3;
    std::vector<double> omega;
    for (int i = 0; i < 6; ++i) { omega.push_back(0.5 + i * 0.9e-3); }
    const std::vector<std::vector<int>> g = LR::group_degenerate_states(omega, thr);
    for (const std::vector<int>& grp : g)
    {
        const double lo = omega[grp.front()];
        const double hi = omega[grp.back()];
        EXPECT_LT(hi - lo, thr);
    }
    EXPECT_GT(g.size(), 1u);
}

/// The template has to work on the type the driver actually uses.
TEST(GradMatrixDegenerate, AssemblesModuleBaseMatrix)
{
    constexpr int nat = 2;
    std::vector<ModuleBase::matrix> diag(2, ModuleBase::matrix(nat, 3));
    diag[0](0, 0) = 1.0;
    diag[1](0, 0) = 3.0;
    std::vector<ModuleBase::matrix> plus(1, ModuleBase::matrix(nat, 3));
    plus[0](0, 0) = 5.0;   // -> off-diagonal 5 - (1+3)/2 = 3
    const std::vector<std::vector<ModuleBase::matrix>> g
        = LR::assemble_grad_matrix(diag, plus, LR::degenerate_pairs(2));
    ASSERT_EQ(g.size(), 2u);
    EXPECT_DOUBLE_EQ(g[0][0](0, 0), 1.0);
    EXPECT_DOUBLE_EQ(g[1][1](0, 0), 3.0);
    EXPECT_DOUBLE_EQ(g[0][1](0, 0), 3.0);
    EXPECT_DOUBLE_EQ(g[1][0](0, 0), 3.0);
}
