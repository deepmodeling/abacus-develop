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
    sc.set_atomCounts(atomCounts);
    sc.zero_Mi();
    testing::internal::CaptureStdout();
    sc.print_termination();
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Inner optimization for lambda ends."));
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
    EXPECT_THAT(output, testing::HasSubstr("Step (Outer -- Inner) =  0 --     1       RMS =1e-05"));
    EXPECT_THAT(output, testing::HasSubstr("Step (Outer -- Inner) =  0 --    12       RMS =1e-05"));
    EXPECT_THAT(output, testing::HasSubstr("Step (Outer -- Inner) =  0 --    13       RMS =1e-07"));
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