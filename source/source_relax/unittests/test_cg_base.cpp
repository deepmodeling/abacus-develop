#include "../cg_base.h"

#include <gtest/gtest.h>

#include <cmath>

// Unit tests for CG_Base. normalize() is protected, so it is exercised
// through a test subclass that re-exposes it (instead of a
// #define protected public access hack).
class TestableCG : public CG_Base
{
  public:
    using CG_Base::normalize;
};

TEST(CGBase, SetupCgGradRestartOnCounter)
{
    // ncggrad % 10000 == 0 triggers the restart branch: cg_grad = grad.
    const int dim = 3;
    double grad[dim] = {1.0, 2.0, 3.0};
    double grad0[dim] = {9.0, 9.0, 9.0};
    double cg_grad[dim] = {0.0, 0.0, 0.0};
    double cg_grad0[dim] = {7.0, 8.0, 9.0};
    int ncggrad = 0;
    int flag = 0;

    TestableCG cg;
    cg.setup_cg_grad(dim, grad, grad0, cg_grad, cg_grad0, ncggrad, flag);

    for (int i = 0; i < dim; ++i)
    {
        EXPECT_DOUBLE_EQ(cg_grad[i], grad[i]);
    }
}

TEST(CGBase, SetupCgGradRestartOnFlagTwo)
{
    // flag == 2 forces the restart branch regardless of ncggrad.
    const int dim = 2;
    double grad[dim] = {4.0, 5.0};
    double grad0[dim] = {1.0, 1.0};
    double cg_grad[dim] = {0.0, 0.0};
    double cg_grad0[dim] = {2.0, 2.0};
    int ncggrad = 5; // not a restart multiple
    int flag = 2;

    TestableCG cg;
    cg.setup_cg_grad(dim, grad, grad0, cg_grad, cg_grad0, ncggrad, flag);

    for (int i = 0; i < dim; ++i)
    {
        EXPECT_DOUBLE_EQ(cg_grad[i], grad[i]);
    }
}

TEST(CGBase, SetupCgGradFRBranch)
{
    // gamma1 = gg/gp_gp < 0.5 selects the Fletcher-Reeves coefficient.
    // grad0 = {1,0}: gp_gp = 1. grad = {0.5,0}: gg = 0.25 -> gamma1 = 0.25.
    const int dim = 2;
    double grad[dim] = {0.5, 0.0};
    double grad0[dim] = {1.0, 0.0};
    double cg_grad[dim] = {0.0, 0.0};
    double cg_grad0[dim] = {2.0, 0.0};
    int ncggrad = 1;
    int flag = 0;

    TestableCG cg;
    cg.setup_cg_grad(dim, grad, grad0, cg_grad, cg_grad0, ncggrad, flag);

    const double gamma = 0.25;
    EXPECT_DOUBLE_EQ(cg_grad[0], grad[0] + gamma * cg_grad0[0]);
    EXPECT_DOUBLE_EQ(cg_grad[1], grad[1] + gamma * cg_grad0[1]);
}

TEST(CGBase, SetupCgGradPRPBranch)
{
    // gamma1 = gg/gp_gp >= 0.5 selects the Polak-Ribiere coefficient
    // gamma2 = (gg - g_gp)/gp_gp.
    // grad0 = {1,0}: gp_gp = 1. grad = {1,0}: gg = 1, g_gp = 1
    //   -> gamma1 = 1 (>= 0.5), gamma2 = 0.
    const int dim = 2;
    double grad[dim] = {1.0, 0.0};
    double grad0[dim] = {1.0, 0.0};
    double cg_grad[dim] = {0.0, 0.0};
    double cg_grad0[dim] = {3.0, 0.0};
    int ncggrad = 1;
    int flag = 0;

    TestableCG cg;
    cg.setup_cg_grad(dim, grad, grad0, cg_grad, cg_grad0, ncggrad, flag);

    const double gamma = 0.0;
    EXPECT_DOUBLE_EQ(cg_grad[0], grad[0] + gamma * cg_grad0[0]);
    EXPECT_DOUBLE_EQ(cg_grad[1], grad[1] + gamma * cg_grad0[1]);
}

TEST(CGBase, FCalProjectsNormalized)
{
    // f_value = (g0 . g1) / |g0|.
    const int dim = 2;
    double g0[dim] = {3.0, 4.0}; // |g0| = 5
    double g1[dim] = {1.0, 0.0};
    double f_value = 0.0;

    TestableCG cg;
    cg.f_cal(dim, g0, g1, f_value);

    EXPECT_DOUBLE_EQ(f_value, 0.6); // 3/5
}

TEST(CGBase, SetupMoveNegativeScaled)
{
    const int dim = 3;
    double move[dim] = {0.0, 0.0, 0.0};
    double cg_gradn[dim] = {1.0, 0.0, -2.0};
    double trust_radius = 0.5;

    TestableCG cg;
    cg.setup_move(dim, move, cg_gradn, trust_radius);

    EXPECT_DOUBLE_EQ(move[0], -0.5);
    EXPECT_DOUBLE_EQ(move[1], 0.0);
    EXPECT_DOUBLE_EQ(move[2], 1.0);
}

TEST(CGBase, NormalizeUnitVector)
{
    const int dim = 2;
    double cg_grad[dim] = {3.0, 4.0};
    double cg_gradn[dim] = {0.0, 0.0};

    TestableCG cg;
    cg.normalize(dim, cg_gradn, cg_grad);

    EXPECT_DOUBLE_EQ(cg_gradn[0], 0.6);
    EXPECT_DOUBLE_EQ(cg_gradn[1], 0.8);
    const double norm = std::sqrt(cg_gradn[0] * cg_gradn[0] + cg_gradn[1] * cg_gradn[1]);
    EXPECT_DOUBLE_EQ(norm, 1.0);
}

TEST(CGBase, NormalizeZeroVectorSafe)
{
    // A zero input must not divide by zero; cg_gradn is left untouched.
    const int dim = 2;
    double cg_grad[dim] = {0.0, 0.0};
    double cg_gradn[dim] = {1.0, 1.0};

    TestableCG cg;
    cg.normalize(dim, cg_gradn, cg_grad);

    EXPECT_DOUBLE_EQ(cg_gradn[0], 1.0);
    EXPECT_DOUBLE_EQ(cg_gradn[1], 1.0);
}

TEST(CGBase, BrentSameSignBranch)
{
    // fa*fb > 0: linear extrapolation, dmove = (xc*fa - xa*fc)/(fa - fc),
    // best_x = dmove - xpt, then xpt = xc = dmove, xb/fb shift to xc/fc.
    TestableCG cg;
    double fa = 1.0;
    double fb = 2.0;
    double fc = 0.5;
    double xa = 0.0;
    double xb = 1.0;
    double xc = 0.5;
    double best_x = 0.0;
    double xpt = 0.0;

    cg.Brent(fa, fb, fc, xa, xb, xc, best_x, xpt);

    const double dmove = (0.5 * 1.0 - 0.0 * 0.5) / (1.0 - 0.5); // = 1.0
    EXPECT_DOUBLE_EQ(best_x, dmove);      // xpt was 0
    EXPECT_DOUBLE_EQ(xpt, dmove);
    EXPECT_DOUBLE_EQ(xc, dmove);
    EXPECT_DOUBLE_EQ(xb, 0.5);            // old xc
    EXPECT_DOUBLE_EQ(fb, 0.5);            // old fc
}

TEST(CGBase, BrentOppositeSignQuadratic)
{
    // fa*fb < 0 fits a quadratic through (xa,fa),(xb,fb),(xc,fc) and picks the
    // stationary point with the lower integrated value. Only sanity-check the
    // invariant: xpt and xc move to the returned displacement + old xpt.
    TestableCG cg;
    double fa = 1.0;
    double fb = -1.0;
    double fc = 0.5;
    double xa = 0.0;
    double xb = 1.0;
    double xc = 0.5;
    double best_x = 0.0;
    double xpt = 0.0;

    cg.Brent(fa, fb, fc, xa, xb, xc, best_x, xpt);

    EXPECT_DOUBLE_EQ(xpt, best_x);  // xpt = dmove = best_x + old xpt(0)
    EXPECT_DOUBLE_EQ(xc, xpt);
}

TEST(CGBase, ThirdOrderFallbackDmoveh)
{
    // When |k3/k1| < 0.01 the cubic term is negligible and the linear
    // estimate dmoveh = x*fb/(fa - fb) is used. Choose fa ~ -fb so the cubic
    // coefficient stays tiny.
    TestableCG cg;
    double best_x = 0.0;
    const double x = 1.0;
    const double fa = -1.0;
    const double fb = 1.0;
    // e0, e1 picked so that e1 - e0 = (fa+fb)*x/2 = 0 makes k3 = 0 exactly.
    const double e0 = 0.0;
    const double e1 = 0.0;

    cg.third_order(e0, e1, fa, fb, x, best_x);

    const double dmoveh = x * fb / (fa - fb); // = -0.5
    EXPECT_DOUBLE_EQ(best_x, dmoveh);
}

TEST(CGBase, ThirdOrderCubicFit)
{
    // Choose inputs that take the genuine cubic-fit branch: k3 != 0, k3 > 0,
    // a non-negative discriminant, and |k3/k1| >= 0.01 so neither fallback
    // fires. With x=1, fa=-1, fb=1, e1-e0=-0.5:
    //   k3 = 3*((fb+fa)*x - 2*(e1-e0)) / x^3 = 3
    //   k2 = (fb-fa)/x - k3*x = -1
    //   k1 = fa = -1
    TestableCG cg;
    double best_x = 0.0;
    const double x = 1.0;
    const double fa = -1.0;
    const double fb = 1.0;
    const double e0 = 0.0;
    const double e1 = -0.5;

    cg.third_order(e0, e1, fa, fb, x, best_x);

    // best_x must be one of the two stationary points minus x, whichever gives
    // the lower cubic energy.
    const double k3 = 3.0;
    const double k2 = -1.0;
    const double k1 = -1.0;
    const double disc = std::sqrt(1.0 - 4.0 * k1 * k3 / (k2 * k2)); // sqrt(13)
    const double dmove1 = -k2 * (1.0 - disc) / (2.0 * k3);
    const double dmove2 = -k2 * (1.0 + disc) / (2.0 * k3);
    const double ecal1 = k3 * dmove1 * dmove1 * dmove1 / 3.0 + k2 * dmove1 * dmove1 / 2.0 + k1 * dmove1;
    const double ecal2 = k3 * dmove2 * dmove2 * dmove2 / 3.0 + k2 * dmove2 * dmove2 / 2.0 + k1 * dmove2;
    const double expected = (ecal2 > ecal1) ? (dmove1 - x) : (dmove2 - x);

    EXPECT_TRUE(std::isfinite(best_x));
    EXPECT_DOUBLE_EQ(best_x, expected);
}

// --- Cases merged from test_ions_move_cg.cpp / test_lattice_change_cg.cpp ---
// These exercise the CG_Base line-search helpers directly (they were
// previously tested redundantly through the Ions_Move_CG and
// Lattice_Change_CG subclasses). References are pre-verified numerical
// outputs of the algorithm.

TEST(CGBase, BrentSameSignMovesToXcCase)
{
    // fa, fb, fc all > 0 with fa*fb > 0: pure linear extrapolation moves the
    // bracket to xc.
    TestableCG cg;
    double fa = 2.0;
    double fb = 1.0;
    double fc = 1.0;
    double xa = -3.0;
    double xb = 2.0;
    double xc = 1.0;
    double best_x = 0.0;
    double xpt = 0.0;

    cg.Brent(fa, fb, fc, xa, xb, xc, best_x, xpt);

    EXPECT_DOUBLE_EQ(fa, 2.0);
    EXPECT_DOUBLE_EQ(xb, 1.0);
    EXPECT_DOUBLE_EQ(xc, 4.0);
    EXPECT_DOUBLE_EQ(best_x, 4.0);
    EXPECT_DOUBLE_EQ(xpt, 4.0);
}

TEST(CGBase, BrentQuadraticCaseA)
{
    TestableCG cg;
    double fa = -2.0;
    double fb = 3.0;
    double fc = -4.0;
    double xa = 1.0;
    double xb = 2.0;
    double xc = 3.0;
    double best_x = 0.0;
    double xpt = 0.0;

    cg.Brent(fa, fb, fc, xa, xb, xc, best_x, xpt);

    EXPECT_DOUBLE_EQ(fa, -4.0);
    EXPECT_DOUBLE_EQ(xa, 3.0);
    EXPECT_DOUBLE_EQ(xb, 2.0);
    EXPECT_NEAR(xc, 1.2046663545568725, 1e-12);
    EXPECT_NEAR(best_x, 1.2046663545568725, 1e-12);
    EXPECT_NEAR(xpt, 1.2046663545568725, 1e-12);
}

TEST(CGBase, BrentQuadraticCaseB)
{
    TestableCG cg;
    double fa = 1.0;
    double fb = -3.0;
    double fc = -4.0;
    double xa = 3.0;
    double xb = 2.0;
    double xc = 1.0;
    double best_x = 0.0;
    double xpt = 0.0;

    cg.Brent(fa, fb, fc, xa, xb, xc, best_x, xpt);

    EXPECT_DOUBLE_EQ(fa, 1.0);
    EXPECT_DOUBLE_EQ(xa, 3.0);
    EXPECT_DOUBLE_EQ(xb, 1.0);
    EXPECT_NEAR(xc, 2.8081429669660172, 1e-12);
    EXPECT_NEAR(best_x, 2.8081429669660172, 1e-12);
    EXPECT_NEAR(xpt, 2.8081429669660172, 1e-12);
}

TEST(CGBase, BrentQuadraticCaseC)
{
    TestableCG cg;
    double fa = 2.0;
    double fb = -3.0;
    double fc = 4.0;
    double xa = 0.0;
    double xb = 2.0;
    double xc = 1.0;
    double best_x = 0.0;
    double xpt = 0.0;

    cg.Brent(fa, fb, fc, xa, xb, xc, best_x, xpt);

    EXPECT_DOUBLE_EQ(fa, 4.0);
    EXPECT_DOUBLE_EQ(xa, 1.0);
    EXPECT_DOUBLE_EQ(xb, 2.0);
    EXPECT_DOUBLE_EQ(xc, 2.0);
    EXPECT_DOUBLE_EQ(best_x, 2.0);
    EXPECT_DOUBLE_EQ(xpt, 2.0);
}

TEST(CGBase, ThirdOrderLinearFallbackPositiveFa)
{
    // |k3/k1| small: linear estimate dmoveh = x*fb/(fa-fb) is returned.
    TestableCG cg;
    double best_x = -1.0;

    cg.third_order(1.0, 1.0, 10.0, -9.99, 1.0, best_x);

    EXPECT_DOUBLE_EQ(best_x, 1.0 * -9.99 / (10.0 - -9.99));
}

TEST(CGBase, ThirdOrderLinearFallbackNegativeFa)
{
    TestableCG cg;
    double best_x = -1.0;

    cg.third_order(1.0, 1.0, -10.0, 9.9, 1.0, best_x);

    EXPECT_DOUBLE_EQ(best_x, 1.0 * 9.9 / (-10.0 - 9.9));
}

TEST(CGBase, ThirdOrderLinearFallbackMixedSign)
{
    TestableCG cg;
    double best_x = -1.0;

    cg.third_order(1.0, 1.0, 10.0, -10.1, 1.0, best_x);

    EXPECT_DOUBLE_EQ(best_x, 1.0 * -10.1 / (10.0 - -10.1));
}

TEST(CGBase, FCalUniformNine)
{
    // g0 = g1 = all-ones in 9 dims: f_value = (9)/sqrt(9) = 3.
    TestableCG cg;
    const int dim = 9;
    double g0[dim] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double g1[dim] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double f_value = 0.0;

    cg.f_cal(dim, g0, g1, f_value);

    EXPECT_DOUBLE_EQ(f_value, 3.0);
}

TEST(CGBase, SetupMoveUniformNine)
{
    TestableCG cg;
    const int dim = 9;
    double trust_radius = 1.0;
    double cg_gradn[dim] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double move[dim] = {0.0};

    cg.setup_move(dim, move, cg_gradn, trust_radius);

    for (int i = 0; i < dim; ++i)
    {
        EXPECT_DOUBLE_EQ(move[i], -1.0);
    }
}
