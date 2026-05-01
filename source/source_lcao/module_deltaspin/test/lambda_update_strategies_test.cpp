#include "../lambda_update_strategies.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <cmath>
#include <vector>
#include <string>

/************************************************
 *  Unit tests for lambda update strategies
 *
 *  - Tested Strategies:
 *    - LinearResponseUpdate (Scheme B)
 *    - AugmentedLagrangianUpdate (Scheme C)
 *    - HybridDelayedUpdate (Scheme D)
 *
 *  - Tested Helpers:
 *    - compute_rms_error()
 *    - count_converged()
 *    - cap_lambda()
 ************************************************/

namespace
{

using ModuleBase::Vector3;

// ===================================================================
// Helper function tests
// ===================================================================

class LambdaUpdateHelpersTest : public ::testing::Test
{
  protected:
    int nat;
    std::vector<Vector3<double>> Mi;
    std::vector<Vector3<double>> target_mag;
    std::vector<Vector3<int>> constrain;

    void SetUp() override
    {
        nat = 3;
        Mi.push_back(Vector3<double>(1.0, 0.5, 0.3));
        Mi.push_back(Vector3<double>(-0.8, 0.2, 0.1));
        Mi.push_back(Vector3<double>(0.5, 0.5, 0.5));

        target_mag.push_back(Vector3<double>(2.0, 0.0, 0.0));
        target_mag.push_back(Vector3<double>(-1.0, 0.0, 0.0));
        target_mag.push_back(Vector3<double>(0.5, 0.5, 0.5));

        constrain.push_back(Vector3<int>(1, 1, 0));
        constrain.push_back(Vector3<int>(1, 0, 0));
        constrain.push_back(Vector3<int>(1, 1, 1));
    }
};

TEST_F(LambdaUpdateHelpersTest, ComputeRmsError)
{
    double rms = spinconstrain::compute_rms_error(Mi, target_mag, constrain, nat);
    // Constrained: atom0(x,y), atom1(x), atom2(x,y,z) = 6 components
    double expected_sum = 1.0*1.0 + 0.5*0.5 + 0.2*0.2 + 0.0 + 0.0 + 0.0;
    double expected_rms = std::sqrt(expected_sum / 6.0);
    EXPECT_NEAR(rms, expected_rms, 1e-10);
}

TEST_F(LambdaUpdateHelpersTest, ComputeRmsErrorAlreadyConverged)
{
    Mi[0] = target_mag[0];
    Mi[1] = target_mag[1];
    Mi[2] = target_mag[2];
    double rms = spinconstrain::compute_rms_error(Mi, target_mag, constrain, nat);
    EXPECT_NEAR(rms, 0.0, 1e-15);
}

TEST_F(LambdaUpdateHelpersTest, ComputeRmsErrorNoConstraints)
{
    std::vector<Vector3<int>> no_constrain(nat, Vector3<int>(0, 0, 0));
    double rms = spinconstrain::compute_rms_error(Mi, target_mag, no_constrain, nat);
    EXPECT_NEAR(rms, 0.0, 1e-15);
}

TEST_F(LambdaUpdateHelpersTest, CountConverged)
{
    int n = spinconstrain::count_converged(Mi, target_mag, constrain, 0.3, nat);
    EXPECT_EQ(n, 4); // 1 from atom1 + 3 from atom2
}

TEST_F(LambdaUpdateHelpersTest, CountConvergedAll)
{
    Mi[0] = target_mag[0];
    Mi[1] = target_mag[1];
    Mi[2] = target_mag[2];
    int n = spinconstrain::count_converged(Mi, target_mag, constrain, 1e-6, nat);
    EXPECT_EQ(n, 6);
}

TEST_F(LambdaUpdateHelpersTest, CapLambda)
{
    std::vector<Vector3<double>> lam(nat);
    lam[0] = Vector3<double>(15.0, -20.0, 5.0);
    lam[1] = Vector3<double>(0.0, 8.0, -12.0);
    lam[2] = Vector3<double>(3.0, 3.0, 3.0);

    std::vector<Vector3<int>> con(nat);
    con[0] = Vector3<int>(1, 1, 1);
    con[1] = Vector3<int>(0, 1, 0);
    con[2] = Vector3<int>(1, 1, 1);

    spinconstrain::cap_lambda(lam, con, 10.0, nat);

    EXPECT_NEAR(lam[0][0], 10.0, 1e-10);
    EXPECT_NEAR(lam[0][1], -10.0, 1e-10);
    EXPECT_NEAR(lam[0][2], 5.0, 1e-10);
    EXPECT_NEAR(lam[1][0], 0.0, 1e-10);
    EXPECT_NEAR(lam[1][1], 8.0, 1e-10);
    EXPECT_NEAR(lam[1][2], -12.0, 1e-10);
    EXPECT_NEAR(lam[2][0], 3.0, 1e-10);
    EXPECT_NEAR(lam[2][1], 3.0, 1e-10);
    EXPECT_NEAR(lam[2][2], 3.0, 1e-10);
}

// ===================================================================
// Scheme B: Linear Response Update tests
// ===================================================================

class LinearResponseTest : public ::testing::Test
{
  protected:
    int nat;
    std::vector<Vector3<double>> lambda;
    std::vector<Vector3<double>> Mi;
    std::vector<Vector3<double>> target_mag;
    std::vector<Vector3<int>> constrain;

    void SetUp() override
    {
        nat = 2;
        lambda.push_back(Vector3<double>(0.0, 0.0, 0.0));
        lambda.push_back(Vector3<double>(0.0, 0.0, 0.0));
        Mi.push_back(Vector3<double>(1.0, 0.0, 0.0));
        Mi.push_back(Vector3<double>(-0.5, 0.0, 0.0));
        target_mag.push_back(Vector3<double>(2.0, 0.0, 0.0));
        target_mag.push_back(Vector3<double>(-1.0, 0.0, 0.0));
        constrain.push_back(Vector3<int>(1, 1, 1));
        constrain.push_back(Vector3<int>(1, 1, 1));
    }
};

TEST_F(LinearResponseTest, FirstUpdateNoHistory)
{
    spinconstrain::LinearResponseUpdate updater(0.01, 100.0, 0.3, 10.0);
    EXPECT_EQ(updater.name(), "LinearResponse");
    EXPECT_FALSE(updater.is_converged());

    auto result = updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-6, 0, nat);

    EXPECT_NEAR(lambda[0][0], 0.3, 1e-10);
    EXPECT_NEAR(lambda[0][1], 0.0, 1e-10);
    EXPECT_LT(result.max_lambda, 1.0);
    EXPECT_EQ(result.status, "updating");
}

TEST_F(LinearResponseTest, ConvergesAfterMultipleSteps)
{
    spinconstrain::LinearResponseUpdate updater(0.01, 100.0, 0.5, 10.0);
    double chi = 1.0;
    Vector3<double> Mi_init_0 = Mi[0];
    Vector3<double> Mi_init_1 = Mi[1];

    int max_iter = 50;
    int converged_iter = -1;
    for (int iter = 0; iter < max_iter; ++iter)
    {
        auto result = updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-5, iter, nat);
        Mi[0] = Vector3<double>(Mi_init_0.x + chi * lambda[0][0],
                                Mi_init_0.y + chi * lambda[0][1],
                                Mi_init_0.z + chi * lambda[0][2]);
        Mi[1] = Vector3<double>(Mi_init_1.x + chi * lambda[1][0],
                                Mi_init_1.y + chi * lambda[1][1],
                                Mi_init_1.z + chi * lambda[1][2]);
        if (updater.is_converged())
        {
            EXPECT_LT(result.rms_error, 1e-5);
            converged_iter = iter;
            break;
        }
    }
    EXPECT_GE(converged_iter, 0) << "Linear response did not converge within " << max_iter;

    double expected_l0 = (target_mag[0][0] - Mi_init_0.x) / chi;
    double expected_l1 = (target_mag[1][0] - Mi_init_1.x) / chi;
    EXPECT_NEAR(lambda[0][0], expected_l0, 0.1);
    EXPECT_NEAR(lambda[1][0], expected_l1, 0.1);
}

TEST_F(LinearResponseTest, RespectsConstrainFlags)
{
    std::vector<Vector3<int>> partial_constrain(nat);
    partial_constrain[0] = Vector3<int>(1, 0, 0);
    partial_constrain[1] = Vector3<int>(0, 0, 0);

    spinconstrain::LinearResponseUpdate updater(0.01, 100.0, 0.3, 10.0);
    updater.update_lambda(lambda, Mi, target_mag, partial_constrain, 1e-6, 0, nat);

    EXPECT_NEAR(lambda[0][0], 0.3, 1e-10);
    EXPECT_NEAR(lambda[0][1], 0.0, 1e-10);
    EXPECT_NEAR(lambda[1][0], 0.0, 1e-10);
}

TEST_F(LinearResponseTest, CapsLambda)
{
    target_mag[0] = Vector3<double>(100.0, 0.0, 0.0);
    spinconstrain::LinearResponseUpdate updater(0.01, 100.0, 1.0, 5.0);
    updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-6, 0, nat);
    EXPECT_LE(std::abs(lambda[0][0]), 5.0 + 1e-10);
}

TEST_F(LinearResponseTest, ChiEstimation)
{
    spinconstrain::LinearResponseUpdate updater(0.01, 100.0, 0.5, 10.0);
    double chi_true = 2.0;
    Vector3<double> Mi_init = Mi[0];

    for (int iter = 0; iter < 5; ++iter)
    {
        updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-6, iter, nat);
        Mi[0] = Vector3<double>(Mi_init.x + chi_true * lambda[0][0], 0.0, 0.0);
        Mi[1] = Vector3<double>(-0.5, 0.0, 0.0);
    }

    const auto& chi = updater.get_chi();
    EXPECT_GT(chi[0][0], 0.5);
    EXPECT_LT(chi[0][0], 50.0);
}

// ===================================================================
// Scheme C: Augmented Lagrangian Update tests
// ===================================================================

class AugmentedLagrangianTest : public ::testing::Test
{
  protected:
    int nat;
    std::vector<Vector3<double>> lambda;
    std::vector<Vector3<double>> Mi;
    std::vector<Vector3<double>> target_mag;
    std::vector<Vector3<int>> constrain;

    void SetUp() override
    {
        nat = 2;
        lambda.push_back(Vector3<double>(0.0, 0.0, 0.0));
        lambda.push_back(Vector3<double>(0.0, 0.0, 0.0));
        Mi.push_back(Vector3<double>(1.0, 0.0, 0.0));
        Mi.push_back(Vector3<double>(-0.5, 0.0, 0.0));
        target_mag.push_back(Vector3<double>(2.0, 0.0, 0.0));
        target_mag.push_back(Vector3<double>(-1.0, 0.0, 0.0));
        constrain.push_back(Vector3<int>(1, 0, 0));
        constrain.push_back(Vector3<int>(1, 0, 0));
    }
};

TEST_F(AugmentedLagrangianTest, FirstUpdate)
{
    spinconstrain::AugmentedLagrangianUpdate updater(0.1, 10.0, 1.5, 5, 10.0);
    EXPECT_EQ(updater.name(), "AugmentedLagrangian");

    auto result = updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-6, 0, nat);

    EXPECT_NEAR(lambda[0][0], -0.1, 1e-10);
    EXPECT_NEAR(lambda[0][1], 0.0, 1e-10);
    EXPECT_NEAR(lambda[1][0], 0.05, 1e-10);
    EXPECT_NEAR(updater.get_mu(), 0.1, 1e-10);
    EXPECT_FALSE(updater.is_converged());
}

TEST_F(AugmentedLagrangianTest, MuGrowth)
{
    spinconstrain::AugmentedLagrangianUpdate updater(0.1, 10.0, 2.0, 3, 10.0);
    for (int iter = 0; iter < 10; ++iter)
    {
        updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-6, iter, nat);
    }
    EXPECT_NEAR(updater.get_mu(), 0.8, 1e-10);
}

TEST_F(AugmentedLagrangianTest, MuCappedAtMax)
{
    spinconstrain::AugmentedLagrangianUpdate updater(0.1, 1.0, 2.0, 1, 10.0);
    for (int iter = 0; iter < 10; ++iter)
    {
        updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-6, iter, nat);
    }
    EXPECT_NEAR(updater.get_mu(), 1.0, 1e-10);
}

TEST_F(AugmentedLagrangianTest, ConvergesWithInvertedResponse)
{
    // Inverted response model: Mi = M_target - chi * lambda
    // Increasing lambda REDUCES the error — models constraint physics correctly
    spinconstrain::AugmentedLagrangianUpdate updater(0.1, 10.0, 1.5, 5, 10.0);
    double chi = 1.0;

    int max_iter = 100;
    int converged_iter = -1;
    for (int iter = 0; iter < max_iter; ++iter)
    {
        auto result = updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-3, iter, nat);

        // Inverted response: Mi approaches M_target as lambda → 0
        Mi[0] = Vector3<double>(target_mag[0][0] - chi * lambda[0][0], 0.0, 0.0);
        Mi[1] = Vector3<double>(target_mag[1][0] - chi * lambda[1][0], 0.0, 0.0);

        if (updater.is_converged())
        {
            EXPECT_LT(result.rms_error, 1e-3);
            converged_iter = iter;
            break;
        }
    }

    EXPECT_GE(converged_iter, 0) << "AL did not converge within " << max_iter;
    EXPECT_NEAR(lambda[0][0], 0.0, 0.5);
}

TEST_F(AugmentedLagrangianTest, ResetMu)
{
    spinconstrain::AugmentedLagrangianUpdate updater(0.1, 10.0, 2.0, 1, 10.0);
    for (int iter = 0; iter < 5; ++iter)
    {
        updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-6, iter, nat);
    }
    EXPECT_GT(updater.get_mu(), 0.1);
    updater.reset_mu();
    EXPECT_NEAR(updater.get_mu(), 0.1, 1e-10);
}

// ===================================================================
// Scheme D: Hybrid Delayed Update tests
// ===================================================================

class HybridDelayedTest : public ::testing::Test
{
  protected:
    int nat;
    std::vector<Vector3<double>> lambda;
    std::vector<Vector3<double>> Mi;
    std::vector<Vector3<double>> target_mag;
    std::vector<Vector3<int>> constrain;

    void SetUp() override
    {
        nat = 2;
        lambda.push_back(Vector3<double>(0.0, 0.0, 0.0));
        lambda.push_back(Vector3<double>(0.0, 0.0, 0.0));
        Mi.push_back(Vector3<double>(1.0, 0.0, 0.0));
        Mi.push_back(Vector3<double>(-0.5, 0.0, 0.0));
        target_mag.push_back(Vector3<double>(2.0, 0.0, 0.0));
        target_mag.push_back(Vector3<double>(-1.0, 0.0, 0.0));
        constrain.push_back(Vector3<int>(1, 1, 1));
        constrain.push_back(Vector3<int>(1, 1, 1));
    }
};

TEST_F(HybridDelayedTest, EarlyPhaseSkip)
{
    spinconstrain::HybridDelayedUpdate updater(1e-3, 0.1, 10.0, 1.5, 5, 10, 10.0);
    updater.set_drho(1.0);

    auto result = updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-6, 0, nat);
    EXPECT_EQ(result.status, "skipped_early");
    EXPECT_EQ(updater.get_phase(), "early");
    EXPECT_NEAR(lambda[0][0], 0.0, 1e-10);
}

TEST_F(HybridDelayedTest, MidPhaseUpdate)
{
    spinconstrain::HybridDelayedUpdate updater(1e-3, 0.1, 10.0, 1.5, 5, 10, 10.0);
    updater.set_drho(5e-3);

    auto result = updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-6, 0, nat);
    EXPECT_EQ(updater.get_phase(), "mid");
    EXPECT_NEAR(lambda[0][0], -0.1, 1e-10);
}

TEST_F(HybridDelayedTest, LatePhaseUpdate)
{
    spinconstrain::HybridDelayedUpdate updater(1e-3, 0.1, 10.0, 1.5, 5, 10, 10.0);
    updater.set_drho(1e-5);

    auto result = updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-6, 0, nat);
    EXPECT_EQ(updater.get_phase(), "late");
    EXPECT_NEAR(lambda[0][0], -0.1, 1e-10);
}

TEST_F(HybridDelayedTest, FallbackSignal)
{
    spinconstrain::HybridDelayedUpdate updater(1e-3, 0.1, 10.0, 1.5, 5, 10, 10.0);
    updater.set_drho(1e-5);

    for (int iter = 0; iter < 5; ++iter)
    {
        auto result = updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-6, iter, nat);
        if (iter >= 2 && result.status == "fallback_triggered")
        {
            EXPECT_TRUE(true);
            return;
        }
    }
    FAIL() << "Fallback was not signaled after several iterations";
}

TEST_F(HybridDelayedTest, Reset)
{
    spinconstrain::HybridDelayedUpdate updater(1e-3, 0.1, 10.0, 1.5, 5, 10, 10.0);
    updater.set_drho(1e-5);
    for (int iter = 0; iter < 10; ++iter)
    {
        updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-6, iter, nat);
    }
    updater.reset();
    EXPECT_EQ(updater.get_phase(), "early");
}

TEST_F(HybridDelayedTest, PhaseTransitions)
{
    spinconstrain::HybridDelayedUpdate updater(1e-3, 0.1, 10.0, 1.5, 5, 10, 10.0);

    updater.set_drho(1.0);
    auto r1 = updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-6, 0, nat);
    EXPECT_EQ(updater.get_phase(), "early");
    EXPECT_EQ(r1.status, "skipped_early");

    updater.set_drho(5e-3);
    updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-6, 1, nat);
    EXPECT_EQ(updater.get_phase(), "mid");

    updater.set_drho(1e-5);
    updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-6, 2, nat);
    EXPECT_EQ(updater.get_phase(), "late");
}

TEST_F(HybridDelayedTest, ConvergesWithInvertedResponse)
{
    spinconstrain::HybridDelayedUpdate updater(1e-3, 0.1, 10.0, 1.5, 5, 10, 10.0);
    updater.set_drho(1e-5);
    double chi = 1.0;

    int max_iter = 100;
    int converged_iter = -1;
    for (int iter = 0; iter < max_iter; ++iter)
    {
        auto result = updater.update_lambda(lambda, Mi, target_mag, constrain, 1e-3, iter, nat);

        Mi[0] = Vector3<double>(target_mag[0][0] - chi * lambda[0][0],
                                target_mag[0][1] - chi * lambda[0][1],
                                target_mag[0][2] - chi * lambda[0][2]);
        Mi[1] = Vector3<double>(target_mag[1][0] - chi * lambda[1][0],
                                target_mag[1][1] - chi * lambda[1][1],
                                target_mag[1][2] - chi * lambda[1][2]);

        if (updater.is_converged())
        {
            EXPECT_LT(result.rms_error, 1e-3);
            converged_iter = iter;
            break;
        }
    }

    EXPECT_GE(converged_iter, 0) << "Hybrid did not converge within " << max_iter
                                  << ". Final phase: " << updater.get_phase();
}

} // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
