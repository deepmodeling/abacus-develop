#include "source_relax/ions_move_basic.h"
#include "source_relax/relax_data.h"
#include "gmock/gmock.h"
#define private public
#include "source_io/module_parameter/parameter.h"
#undef private
#include "gtest/gtest.h"
#define private public
#define protected public
#include "source_relax/bfgs_basic.h"
#undef private
#undef protected
/************************************************
 *  unit tests of class BFGS_Basic
 ***********************************************/

class BFGSBasicTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        // Initialize variables before each test
    }

    void TearDown() override
    {
        // Clean up after each test
    }

    BFGS_Basic bfgs;
};

// Test whether the allocate_basic() function can correctly allocate memory space
TEST_F(BFGSBasicTest, TestAllocate)
{
    Ions_Move_Basic::dim = 4;
    bfgs.allocate_basic();

    // Check if allocated vectors are not empty
    EXPECT_EQ(bfgs.pos.size(), 4U);
    EXPECT_EQ(bfgs.pos_p.size(), 4U);
    EXPECT_EQ(bfgs.grad.size(), 4U);
    EXPECT_EQ(bfgs.grad_p.size(), 4U);
    EXPECT_EQ(bfgs.move.size(), 4U);
    EXPECT_EQ(bfgs.move_p.size(), 4U);
}

// Test if a dimension less than or equal to 0 results in an assertion error
TEST_F(BFGSBasicTest, TestAllocateWithZeroDimension)
{
    Ions_Move_Basic::dim = 0;
    ASSERT_DEATH(bfgs.allocate_basic(), "");
}

// Test function update_inverse_hessian() assert death
TEST_F(BFGSBasicTest, UpdateInverseHessianDeath)
{
    Ions_Move_Basic::dim = 0;
    double lat0 = 1.0;
    std::ofstream ofs("test_log_update_inverse_hessian_death.log");
    ASSERT_DEATH(bfgs.update_inverse_hessian(lat0, ofs), "");
    ofs.close();
    std::remove("test_log_update_inverse_hessian_death.log");
}

// Test function update_inverse_hessian() when sdoty = 0
TEST_F(BFGSBasicTest, UpdateInverseHessianCase1)
{
    Ions_Move_Basic::dim = 3;
    double lat0 = 1.0;
    bfgs.allocate_basic();

    std::ofstream ofs("test_log_update_inverse_hessian_case1.log");
    bfgs.update_inverse_hessian(lat0, ofs);
    ofs.close();

    std::string expected_output
        = " WARINIG: unexpected behaviour in update_inverse_hessian\n Resetting bfgs history \n";
    std::ifstream ifs("test_log_update_inverse_hessian_case1.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("test_log_update_inverse_hessian_case1.log");

    EXPECT_EQ(expected_output, output);
}

// Test function update_inverse_hessian()
TEST_F(BFGSBasicTest, UpdateInverseHessianCase2)
{
    Ions_Move_Basic::dim = 3;
    double lat0 = 1.0;
    bfgs.allocate_basic();
    bfgs.pos[0] = 2.0;
    bfgs.grad[0] = 2.0;

    std::ofstream ofs("test_log_update_inverse_hessian_case2.log");
    bfgs.update_inverse_hessian(lat0, ofs);
    ofs.close();

    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 0), 0.5);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 2), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(2, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(2, 2), 0.0);
    std::remove("test_log_update_inverse_hessian_case2.log");
}

// Test function check_wolfe_conditions()
TEST_F(BFGSBasicTest, CheckWolfeConditions)
{
    Ions_Move_Basic::dim = 3;
    PARAM.input.test_relax_method = 1;
    bfgs.allocate_basic();
    bfgs.pos[0] = 2.0;
    bfgs.grad[0] = 2.0;
    bfgs.move[0] = 1.0;
    std::vector<double> etot_info = {10.0, 0.0};

    std::ofstream ofs("test_log_check_wolfe_conditions.log");
    bfgs.check_wolfe_conditions(ofs, etot_info);
    ofs.close();

    std::string expected_output
        = "                            etot - etot_p = 10\n                    relax_bfgs_w1 * dot_p = -0\n            "
          "                          dot = 0\n                    relax_bfgs_w2 * dot_p = -0\n                         "
          "   relax_bfgs_w1 = -1\n                            relax_bfgs_w2 = -1\n                                   "
          "wolfe1 = 0\n                                   wolfe2 = 0\n                            etot - etot_p = 10\n "
          "                   relax_bfgs_w1 * dot_p = -0\n                                   wolfe1 = 0\n              "
          "                     wolfe2 = 0\n                wolfe condition satisfied = 0\n";

    std::ifstream ifs("test_log_check_wolfe_conditions.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("test_log_check_wolfe_conditions.log");

    EXPECT_EQ(bfgs.wolfe_flag, false);
    EXPECT_EQ(expected_output, output);
}

// Test function reset_hessian()
TEST_F(BFGSBasicTest, ResetHessian)
{
    Ions_Move_Basic::dim = 3;
    bfgs.allocate_basic();

    bfgs.reset_hessian();

    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 2), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(2, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(2, 2), 1.0);
}

// Test function save_bfgs()
TEST_F(BFGSBasicTest, SaveBfgs)
{
    Ions_Move_Basic::dim = 2;
    bfgs.save_flag = false;
    bfgs.allocate_basic();
    bfgs.pos[0] = 1.0;
    bfgs.pos[1] = 2.0;
    bfgs.grad[0] = 3.0;
    bfgs.grad[1] = 4.0;
    bfgs.move[0] = 5.0;
    bfgs.move[1] = 6.0;

    bfgs.save_bfgs();

    EXPECT_EQ(bfgs.save_flag, true);
    EXPECT_DOUBLE_EQ(bfgs.pos[0], 1.0);
    EXPECT_DOUBLE_EQ(bfgs.pos[1], 2.0);
    EXPECT_DOUBLE_EQ(bfgs.grad[0], 3.0);
    EXPECT_DOUBLE_EQ(bfgs.grad[1], 4.0);
    EXPECT_DOUBLE_EQ(bfgs.move[0], 5.0);
    EXPECT_DOUBLE_EQ(bfgs.move[1], 6.0);
}

// Test function new_step() when update_iter == 1
TEST_F(BFGSBasicTest, NewStepCase1)
{
    Ions_Move_Basic::dim = 2;
    int update_iter = 0;
    Ions_Move_Basic::largest_grad = 0.0;
    Ions_Move_Basic::relax_bfgs_init = 0.3;
    Ions_Move_Basic::best_xxx = -0.4;
    bfgs.bfgs_ndim = 1;
    bfgs.allocate_basic();
    bfgs.grad[0] = 1.0;
    bfgs.grad[1] = 2.0;
    bfgs.inv_hess(0, 0) = -3.0;
    bfgs.inv_hess(0, 1) = -4.0;
    bfgs.inv_hess(1, 0) = -5.0;
    bfgs.inv_hess(1, 1) = -6.0;

    double lat0 = 1.0;
    std::ofstream ofs("test_log.log");
    std::vector<double> etot_info(2, 0.0);
    bfgs.new_step(lat0, update_iter, ofs, etot_info);

    EXPECT_EQ(update_iter, 1);
    EXPECT_EQ(bfgs.tr_min_hit, false);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::relax_bfgs_init, 0.2);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::best_xxx, 0.4);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, 0.2);
    EXPECT_DOUBLE_EQ(bfgs.move[0], -1.0);
    EXPECT_DOUBLE_EQ(bfgs.move[1], -2.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 1), 1.0);
}

// Test function new_step() when update_iter > 1
TEST_F(BFGSBasicTest, NewStepCase2)
{
    Ions_Move_Basic::dim = 2;
    int update_iter = 2;
    Ions_Move_Basic::largest_grad = 0.0;
    Ions_Move_Basic::relax_bfgs_init = 0.3;
    Ions_Move_Basic::best_xxx = -0.4;
    bfgs.bfgs_ndim = 1;
    bfgs.allocate_basic();
    bfgs.grad[0] = 1.0;
    bfgs.grad[1] = 2.0;
    bfgs.inv_hess(0, 0) = -3.0;
    bfgs.inv_hess(0, 1) = -4.0;
    bfgs.inv_hess(1, 0) = -5.0;
    bfgs.inv_hess(1, 1) = -6.0;

    double lat0 = 1.0;
    std::ofstream ofs("test_log.log");
    std::vector<double> etot_info(2, 0.0);
    bfgs.new_step(lat0, update_iter, ofs, etot_info);

    EXPECT_EQ(update_iter, 3);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, -1.0);
    EXPECT_DOUBLE_EQ(bfgs.move[0], -1.0);
    EXPECT_DOUBLE_EQ(bfgs.move[1], -2.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 1), 1.0);
}

// Test function new_step() when bfgs_ndim > 1
TEST_F(BFGSBasicTest, NewStepWarningQuit)
{
    Ions_Move_Basic::dim = 2;
    int update_iter = 0;
    bfgs.bfgs_ndim = 2;
    bfgs.allocate_basic();
    double lat0 = 1.0;
    std::ofstream ofs("test_log.log");
    std::vector<double> etot_info(2, 0.0);

    testing::internal::CaptureStdout();
    EXPECT_EXIT(bfgs.new_step(lat0, update_iter, ofs, etot_info), ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("bfgs_ndim > 1 not implemented yet"));
}

// Test function compute_trust_radius() case 1
TEST_F(BFGSBasicTest, ComputeTrustRadiusCase1)
{
    Ions_Move_Basic::dim = 2;
    bfgs.allocate_basic();
    bfgs.grad_p[0] = 1.0;
    bfgs.move_p[1] = 2.0;
    bfgs.inv_hess(0, 0) = -3.0;
    bfgs.inv_hess(0, 1) = -4.0;
    bfgs.inv_hess(1, 0) = -5.0;
    bfgs.inv_hess(1, 1) = -6.0;
    bfgs.wolfe_flag = true;
    bfgs.relax_bfgs_w1 = 1.0;
    std::vector<double> etot_info = {0.0, 0.0, 0.0};

    std::ofstream ofs("test_log.log");
    bfgs.compute_trust_radius(ofs, etot_info);

    EXPECT_EQ(bfgs.tr_min_hit, false);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, -1.0);
    EXPECT_DOUBLE_EQ(bfgs.move[0], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.move[1], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 0), -3.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 1), -4.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 0), -5.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 1), -6.0);
}

// Test function compute_trust_radius() case 2
TEST_F(BFGSBasicTest, ComputeTrustRadiusCase2)
{
    Ions_Move_Basic::dim = 2;
    Ions_Move_Basic::trust_radius_old = 0.0;
    Ions_Move_Basic::relax_bfgs_rmin = 100.0;
    PARAM.input.test_relax_method = 1;
    bfgs.allocate_basic();
    bfgs.grad_p[0] = 1.0;
    bfgs.move[1] = 2.0;
    bfgs.move_p[0] = 2.0;
    bfgs.inv_hess(0, 0) = -3.0;
    bfgs.inv_hess(0, 1) = -4.0;
    bfgs.inv_hess(1, 0) = -5.0;
    bfgs.inv_hess(1, 1) = -6.0;
    bfgs.wolfe_flag = false;
    bfgs.relax_bfgs_w1 = 1.0;
    bfgs.tr_min_hit = false;
    std::vector<double> etot_info = {0.0, 0.0, 0.0};

    std::ofstream ofs("test_log.log");
    bfgs.compute_trust_radius(ofs, etot_info);

    EXPECT_EQ(bfgs.tr_min_hit, true);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, 100.0);
    EXPECT_DOUBLE_EQ(bfgs.move[0], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.move[1], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 1), 1.0);
}

// Test function compute_trust_radius() warning_quit
TEST_F(BFGSBasicTest, ComputeTrustRadiusWarningQuit)
{
    Ions_Move_Basic::dim = 2;
    Ions_Move_Basic::trust_radius_old = 0.0;
    Ions_Move_Basic::relax_bfgs_rmin = 100.0;
    PARAM.input.test_relax_method = 1;
    bfgs.allocate_basic();
    bfgs.grad_p[0] = 1.0;
    bfgs.move[1] = 2.0;
    bfgs.move_p[0] = 2.0;
    bfgs.inv_hess(0, 0) = -3.0;
    bfgs.inv_hess(0, 1) = -4.0;
    bfgs.inv_hess(1, 0) = -5.0;
    bfgs.inv_hess(1, 1) = -6.0;
    bfgs.wolfe_flag = false;
    bfgs.relax_bfgs_w1 = 1.0;
    bfgs.tr_min_hit = true;
    std::vector<double> etot_info = {0.0, 0.0, 0.0};

    std::ofstream ofs("test_log.log");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(bfgs.compute_trust_radius(ofs, etot_info), ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("bfgs history already reset at previous step, we got trapped!"));
}