#include <gtest/gtest.h>
#include "../window_response_math.h"
#include <cmath>
#include <stdexcept>

TEST(WindowResponse, MovingGeneralizedEigenvectorMatchesFiniteDifference)
{
    const double eps[] = {0.7, 2.3};
    const double X[] = {1.0, 0.0};
    const double AX[] = {1.1, 0.31};
    const auto W = LR::window_multipliers(1, 2, 1, eps, X, AX);
    // Generalized F c = S c eps, with off-diagonal F'=0.23, S'=0.17.
    // The truncated excitation includes a fixed AO kernel K.
    const auto energy = [](const double t)
    {
        const double b = 0.23 * t;
        const double s = 0.17 * t;
        const double qa = 1.0 - s * s;
        const double qb = -(0.7 + 2.3) + 2.0 * b * s;
        const double qc = 0.7 * 2.3 - b * b;
        const double eig = (-qb - std::sqrt(qb * qb - 4.0 * qa * qc)) / (2.0 * qa);
        const double ratio = -(b - eig * s) / (2.3 - eig);
        const double norm = 1.0 + 2.0 * s * ratio + ratio * ratio;
        return eig + (0.4 + 2.0 * 0.31 * ratio + 0.2 * ratio * ratio) / norm;
    };
    const double step = 1e-5;
    const double fd = (energy(step) - energy(-step)) / (2.0 * step);
    EXPECT_NEAR(W[1] * (0.23 - eps[0] * 0.17), fd, 1e-10);
    EXPECT_EQ(W[0], 0.0);
    EXPECT_EQ(W[2], 0.0);
    EXPECT_EQ(W[3], 0.0);
}

TEST(WindowResponse, FullWindowVanishesAndDegenerateBoundaryIsRejected)
{
    const double eps[] = {1.0, 1.0};
    const double X[] = {1.0, 0.0};
    const double AX[] = {1.0, 0.3};
    EXPECT_THROW(LR::window_multipliers(1, 2, 1, eps, X, AX), std::runtime_error);
    const auto W = LR::window_multipliers(1, 2, 2, eps, X, AX);
    for (const double value : W)
    {
        EXPECT_EQ(value, 0.0);
    }
}

TEST(WindowResponse, QuadraticAndIndependentOfAmplitudePhase)
{
    const double eps[] = {1.0, 2.0, 3.0};
    const double X[] = {0.3, 0.4, 0.0, 0.5, 0.2, 0.0};
    const double AX[] = {0.1, 0.2, 0.6, 0.3, 0.4, -0.7};
    double scaled_x[6];
    double scaled_ax[6];
    for (int i = 0; i < 6; ++i)
    {
        scaled_x[i] = -2.0 * X[i];
        scaled_ax[i] = -2.0 * AX[i];
    }
    const auto W = LR::window_multipliers(2, 3, 2, eps, X, AX);
    const auto scaled_w = LR::window_multipliers(2, 3, 2, eps, scaled_x, scaled_ax);
    for (int i = 0; i < 9; ++i)
    {
        EXPECT_NEAR(scaled_w[i], 4.0 * W[i], 1e-15);
    }
}


TEST(WindowResponse, CoupledSpinEnergyUsesJointNormalization)
{
    const double up_eps[] = {0.7, 2.3};
    const double down_eps[] = {0.9, 3.1};
    // The retained 2x2 Casida matrix has eigenvalue 1 and vector (0.6, 0.8).
    const double up_x[] = {0.6, 0.0};
    const double down_x[] = {0.8, 0.0};
    const double up_ax[] = {0.6, 0.31 * 0.6 - 0.17 * 0.8};
    const double down_ax[] = {0.8, 0.23 * 0.6 - 0.29 * 0.8};
    const auto up_w = LR::window_multipliers(1, 2, 1, up_eps, up_x, up_ax);
    const auto down_w = LR::window_multipliers(1, 2, 1, down_eps, down_x, down_ax);
    const auto orbital = [](const double t, const double e0, const double e1,
                            const double fprime, const double sprime)
    {
        const double b = fprime * t;
        const double s = sprime * t;
        const double qa = 1.0 - s * s;
        const double qb = -(e0 + e1) + 2.0 * b * s;
        const double qc = e0 * e1 - b * b;
        const double eig = (-qb - std::sqrt(qb * qb - 4.0 * qa * qc)) / (2.0 * qa);
        const double ratio = -(b - eig * s) / (e1 - eig);
        const double norm = std::sqrt(1.0 + 2.0 * s * ratio + ratio * ratio);
        return std::vector<double>{eig, 1.0 / norm, ratio / norm};
    };
    const auto energy = [&orbital](const double t)
    {
        const auto u = orbital(t, 0.7, 2.3, 0.23, 0.17);
        const auto d = orbital(t, 0.9, 3.1, -0.19, 0.11);
        const double auu = u[0] + 1.1 * u[1] * u[1] + 0.62 * u[1] * u[2] + 0.2 * u[2] * u[2];
        const double add = d[0] + 0.55 * d[1] * d[1] - 0.58 * d[1] * d[2] + 0.4 * d[2] * d[2];
        const double aud = -0.6 * u[1] * d[1] + 0.23 * u[1] * d[2]
            - 0.17 * u[2] * d[1] + 0.1 * u[2] * d[2];
        const double delta = auu - add;
        return 0.5 * (auu + add - std::sqrt(delta * delta + 4.0 * aud * aud));
    };
    const double analytic = up_w[1] * (0.23 - up_eps[0] * 0.17)
        + down_w[1] * (-0.19 - down_eps[0] * 0.11);
    const double step = 1e-5;
    const double fd = (energy(step) - energy(-step)) / (2.0 * step);
    EXPECT_NEAR(energy(0.0), 1.0, 1e-14);
    EXPECT_NEAR(analytic, fd, 1e-10);
}
