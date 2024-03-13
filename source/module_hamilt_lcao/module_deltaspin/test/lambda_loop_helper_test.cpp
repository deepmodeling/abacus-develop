#include "../spin_constrain.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

K_Vectors::K_Vectors()
{
}
K_Vectors::~K_Vectors()
{
}

/************************************************
 *  unit test of the functions in lambda_loop_helper.cpp
 ***********************************************/

/**
 * Tested function:
 * - SpinConstrain::check_rms_stop
 *  - check if the rms error is small enough to stop the lambda loop
 * - SpinConstrain::print_termination
 * - print termination message
 */

class SpinConstrainTest : public testing::Test
{
  protected:
    SpinConstrain<std::complex<double>, psi::DEVICE_CPU>& sc
        = SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::getScInstance();
};

TEST_F(SpinConstrainTest, PrintTermination)
{
    std::map<int, int> atomCounts = {
        {0, 1}
    };
    sc.set_nspin(4);
    sc.set_atomCounts(atomCounts);
    sc.zero_Mi();
    std::vector<ModuleBase::Vector3<double>> sc_lambda = std::vector<ModuleBase::Vector3<double>>(1, {1.0, 2.0, 3.0});
    sc.set_sc_lambda(sc_lambda.data(), 1);
    testing::internal::CaptureStdout();
    sc.print_termination();
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Inner optimization for lambda ends."));
    EXPECT_THAT(output, testing::HasSubstr("ATOM 1         0.0000000000e+00    0.0000000000e+00    0.0000000000e+00"));
    EXPECT_THAT(output, testing::HasSubstr("ATOM 1         1.0000000000e+00    2.0000000000e+00    3.0000000000e+00"));
}

TEST_F(SpinConstrainTest, CheckRmsStop)
{
    double sc_thr = 1e-6;
    int nsc = 100;
    int nsc_min = 2;
    double alpha_trial = 0.01;
    double sccut = 3.0;
    bool decay_grad_switch = 1;
    this->sc.set_input_parameters(sc_thr, nsc, nsc_min, alpha_trial, sccut, decay_grad_switch);
    testing::internal::CaptureStdout();
    EXPECT_FALSE(sc.check_rms_stop(0, 0, 1e-5));
    EXPECT_FALSE(sc.check_rms_stop(0, 11, 1e-5));
    EXPECT_TRUE(sc.check_rms_stop(0, 12, 1e-7));
    EXPECT_TRUE(sc.check_rms_stop(0, 99, 1e-5));
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Step (Outer -- Inner) =  0 -- 1           RMS = 1e-05"));
    EXPECT_THAT(output, testing::HasSubstr("Step (Outer -- Inner) =  0 -- 12          RMS = 1e-05"));
    EXPECT_THAT(output, testing::HasSubstr("Step (Outer -- Inner) =  0 -- 13          RMS = 1e-07"));
    EXPECT_THAT(output, testing::HasSubstr("Meet convergence criterion ( < 1e-06 ), exit."));
    EXPECT_THAT(output, testing::HasSubstr("Reach maximum number of steps ( 100 ), exit."));
}

TEST_F(SpinConstrainTest, PrintHeader)
{
    testing::internal::CaptureStdout();
    sc.print_header();
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Inner optimization for lambda begins ..."));
    EXPECT_THAT(output, testing::HasSubstr("Covergence criterion for the iteration: 1e-06"));
}

TEST_F(SpinConstrainTest, CheckRestriction)
{
    std::vector<ModuleBase::Vector3<double>> search = {
        {0.0, 0.0, 40}
    };
    double alpha_trial = 0.1 / 13.605698;
    testing::internal::CaptureStdout();
    sc.check_restriction(search, alpha_trial);
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("alpha after restrict = 0.075"));
    EXPECT_THAT(output, testing::HasSubstr("boundary after = 3"));
}

TEST_F(SpinConstrainTest, CalAlphaOpt)
{
    std::vector<ModuleBase::Vector3<int>> constrain = {
        {1, 1, 1}
    };
    std::vector<ModuleBase::Vector3<double>> target_mag = {
        {0.0, 0.0, 2.0}
    };
    // Set up test input data
    std::vector<ModuleBase::Vector3<double>> spin = {
        {0.0, 0.0, 0.1}
    };

    std::vector<ModuleBase::Vector3<double>> spin_plus = {
        {0.0, 0.0, 0.2}
    };

    sc.set_constrain(constrain.data(), 1);
    sc.set_target_mag(target_mag.data(), 1);

    double alpha_trial = 0.5;

    // Set up expected output data
    double expected_alpha_opt = 9.5;

    // Call the function being tested
    double actual_alpha_opt = sc.cal_alpha_opt(spin, spin_plus, alpha_trial);

    // Compare the expected and actual output
    EXPECT_NEAR(expected_alpha_opt, actual_alpha_opt, 1e-14);
}

TEST_F(SpinConstrainTest, CheckGradientDecay)
{
    // Set up some data for testing
    std::vector<ModuleBase::Vector3<double>> new_spin = {
        {0.0, 0.0, 0.1}
    };

    std::vector<ModuleBase::Vector3<double>> new_spin1 = {
        {0.0, 0.0, 10.0}
    };

    std::vector<ModuleBase::Vector3<double>> spin = {
        {0.0, 0.0, 0.2}
    };

    std::vector<ModuleBase::Vector3<double>> delta_lambda = {
        {0.0, 0.0, 1.0}
    };

    std::vector<ModuleBase::Vector3<double>> dnu_last_step = {
        {0.0, 0.0, 2.0},
    };

    std::vector<ModuleBase::Vector3<int>> constrain = {
        {0, 0, 1}
    };

    std::vector<double> decay_grad = {0.9};
    sc.set_constrain(constrain.data(), 1);
    sc.set_decay_grad(decay_grad.data(), 1);

    // Call the function to test
    EXPECT_TRUE(sc.check_gradient_decay(new_spin, spin, delta_lambda, dnu_last_step, false));
    testing::internal::CaptureStdout();
    EXPECT_FALSE(sc.check_gradient_decay(new_spin1, spin, delta_lambda, dnu_last_step, true));
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("maximum gradient appears at:"));
}