#include "../sto_func.h"

#include <gtest/gtest.h>

#include <cmath>

/************************************************
 *  unit test of sto_func.cpp  (occupation / entropy kernels of SDFT)
 *
 *  Goal: turn the "Ne vs TS Chebyshev precision asymmetry" failure into a
 *  self-checking regression test. Two known failure modes:
 *   (1) computing `1 - f` by subtraction loses accuracy when f -> 1 or f -> 0,
 *       which destroys the relative accuracy of the entropy density g;
 *   (2) an artificial cutoff such as `if (x > 36) return 0` introduces a
 *       discontinuity in the expanded function, which lifts the tail of the
 *       Chebyshev coefficients by several orders of magnitude and is then
 *       amplified by |x| > 1 into an O(1) spurious entropy.
 *
 *  Only mathematical/physical properties are asserted here (no implementation
 *  detail), so any correct implementation passes.
 *
 *  Sign convention of sto_func.cpp (must be respected, otherwise TS flips sign):
 *      fdlnfd(e)           =  f ln f + (1-f) ln(1-f)   <= 0   (i.e. -g)
 *      nfdlnfd(rawe)       =  the same negative value         (i.e. -g)
 *      n_root_fdlnfd(rawe) =  sqrt(g) >= 0                    (expanded by method 2)
 *  with g = -f ln f - (1-f) ln(1-f) >= 0 the true entropy density.
 *  The "n*" variants take the normalized Chebyshev coordinate rawe in [-1,1]
 *  and restore the energy through [Emin, Emax].
 ***********************************************/

namespace
{
const double kT = 25.57754484; // 348 eV, the production temperature [Ry]
const double mu = 0.3;         // arbitrary chemical potential inside the spectrum [Ry]

// Independent closed form (equivalent, but free of cancellation):
//   g(x) = ln(2cosh(x/2)) - (x/2)tanh(x/2) = |x|*y/(1+y) + ln(1+y), y = exp(-|x|)
// It follows a completely different path than "compute f, then take log" in
// ABACUS and agrees with a 50-digit mpmath reference to ~1e-16 relative
// for |x| <= 60.
double g_reference(double x)
{
    const double y = std::exp(-std::fabs(x));
    return std::fabs(x) * y / (1.0 + y) + std::log1p(y);
}
} // namespace

/**
 * - Tested Functions:
 *   - Sto_Func::fd / root_fd            (occupation f and sqrt(f), input = physical energy)
 *   - Sto_Func::nroot_fd                (sqrt(f), input = normalized coordinate)
 *   - Sto_Func::fdlnfd                  (f ln f + (1-f) ln(1-f) = -g, input = physical energy)
 *   - Sto_Func::n_root_fdlnfd           (sqrt(g), input = normalized coordinate)
 */
class TestStoFunc : public ::testing::Test
{
  protected:
    Sto_Func<double> sto;
    double emin = 0.0;
    double emax = 0.0;
    double xc = 0.0; // centre of [xmin,xmax] in units of x=(E-mu)/kT
    double xh = 1.0; // half width

    void SetUp() override
    {
        sto.tem = kT;
        sto.mu = mu;
        set_x_range(-10.0, 50.0);
    }

    // Chebyshev interval used by the "n*" variants: x in [xmin,xmax] <=> rawe in [-1,1]
    void set_x_range(double xmin, double xmax)
    {
        xc = 0.5 * (xmin + xmax);
        xh = 0.5 * (xmax - xmin);
        emin = mu + xmin * kT;
        emax = mu + xmax * kT;
        sto.set_E_range(&emin, &emax);
    }

    double to_rawe(double x) const
    {
        return (x - xc) / xh;
    }
};

// 1) identities and independent reference: f(x)+f(-x)=1, g symmetric, sqrt(g)^2=g, f(0)=1/2, g(0)=ln2
TEST_F(TestStoFunc, identities_and_independent_reference)
{
    EXPECT_NEAR(sto.fd(mu), 0.5, 1e-16);
    EXPECT_NEAR(-sto.fdlnfd(mu), std::log(2.0), 1e-14);
    EXPECT_NEAR(sto.n_root_fdlnfd(to_rawe(0.0)), std::sqrt(std::log(2.0)), 1e-13);

    set_x_range(-20.0, 20.0);
    for (double x = -20.0; x <= 20.0; x += 0.5)
    {
        const double e = mu + x * kT;
        const double e_neg = mu - x * kT;
        const double g = -sto.fdlnfd(e);

        EXPECT_NEAR(sto.fd(e) + sto.fd(e_neg), 1.0, 1e-14);          // particle-hole symmetry
        EXPECT_NEAR(g, -sto.fdlnfd(e_neg), 1e-13);                   // g is even in x
        EXPECT_GE(g, 0.0);                                           // entropy is non-negative
        EXPECT_LE(sto.fdlnfd(e), 0.0);                               // sign convention
        EXPECT_NEAR(g, g_reference(x), 1e-13 * (1.0 + std::fabs(g_reference(x))));
        EXPECT_NEAR(sto.n_root_fdlnfd(to_rawe(x)) * sto.n_root_fdlnfd(to_rawe(x)), g, 1e-13);
        EXPECT_GE(sto.fd(e), 0.0);
        EXPECT_LE(sto.fd(e), 1.0);
    }
}

// 2) key regression: the expanded function must not contain an artificial jump.
//    The true value is still 1e-13 ... 9e-8 for x in (36,45); `if (x>36) return 0`
//    flattens it, leaving a step of sqrt(g(36)) = 9.26e-8 whose Chebyshev
//    coefficient tail is ~2e-9, whereas a correct implementation stays ~1e-15.
TEST_F(TestStoFunc, no_artificial_cutoff_jump)
{
    // (a) for x in (36,45) g must still be positive and match the reference
    for (double x = 36.5; x <= 45.0; x += 0.5)
    {
        const double e = mu + x * kT;
        const double g = -sto.fdlnfd(e);
        const double g_ref = g_reference(x);
        EXPECT_GT(g, 0.0);                     // the original returns 0 here -> FAIL
        EXPECT_NEAR(g, g_ref, 1e-12 * g_ref + 1e-300);
        EXPECT_GT(sto.n_root_fdlnfd(to_rawe(x)), 0.0);   // same for sqrt(g)
    }

    // (b) sqrt(g) is continuous on [33,45]: with a 1e-2 step the difference
    //     between neighbours must stay far below the 9.3e-8 artificial step.
    //     The steepest legitimate slope in this range is at x=33,
    //     |d sqrt(g)/dx| ~ 0.5*sqrt(g(33)) ~ 2.5e-7, i.e. ~2.5e-9 per 1e-2 step.
    double max_step = 0.0;
    double prev = sto.n_root_fdlnfd(to_rawe(33.0));
    for (int i = 1; i <= 1200; ++i)
    {
        const double x = 33.0 + 0.01 * i;
        const double cur = sto.n_root_fdlnfd(to_rawe(x));
        max_step = std::max(max_step, std::fabs(cur - prev));
        EXPECT_LE(cur, prev + 1e-15);          // sqrt(g) is non-increasing for x>0
        prev = cur;
    }
    EXPECT_LE(max_step, 1e-8);                 // artificial step 9.3e-8 -> FAIL
}

// 3) tail asymptotics and extremes: no NaN/Inf, and no crude truncation of a meaningful tail
TEST_F(TestStoFunc, tail_asymptotics_and_extremes)
{
    set_x_range(-800.0, 800.0);
    for (double x: {50.0, 80.0, 120.0, 200.0})
    {
        const double g = -sto.fdlnfd(mu + x * kT);
        const double asym = std::exp(-x) * (x + 1.0);
        EXPECT_GT(g, 0.0);                     // the original is 0 for x>36 -> FAIL
        EXPECT_NEAR(g / asym, 1.0, 1e-9);
    }
    for (double x: {-50.0, -80.0, -120.0, -200.0})
    {
        const double g = -sto.fdlnfd(mu + x * kT);
        const double asym = std::exp(x) * (-x + 1.0);
        EXPECT_GT(g, 0.0);
        EXPECT_NEAR(g / asym, 1.0, 1e-9);
    }
    for (double x: {-700.0, 700.0, -745.0, 745.0, 1000.0, -1000.0})
    {
        const double e = mu + x * kT;
        EXPECT_TRUE(std::isfinite(sto.fd(e)));
        EXPECT_TRUE(std::isfinite(sto.fdlnfd(e)));
        EXPECT_TRUE(std::isfinite(sto.n_root_fdlnfd(to_rawe(x))));
        EXPECT_GE(sto.fd(e), 0.0);
        EXPECT_LE(sto.fd(e), 1.0);
        EXPECT_GE(-sto.fdlnfd(e), 0.0);
    }
}

// 4) monotonicity and sqrt(f) consistency (method 2 uses sqrt(f) for Ne)
TEST_F(TestStoFunc, monotonicity_and_sqrt_f)
{
    set_x_range(-50.0, 50.0);
    double prev = sto.fd(mu - 50.0 * kT);
    for (int i = -5000; i <= 5000; ++i)
    {
        const double x = 0.01 * i;
        const double e = mu + x * kT;
        const double cur = sto.fd(e);
        EXPECT_LE(cur, prev + 1e-16);          // f is non-increasing
        prev = cur;
        const double sf = sto.root_fd(e);
        EXPECT_GE(sf, 0.0);
        EXPECT_NEAR(sf * sf, cur, 1e-14);
        const double nsf = sto.nroot_fd(to_rawe(x));
        EXPECT_NEAR(nsf * nsf, cur, 1e-13);
    }
}
