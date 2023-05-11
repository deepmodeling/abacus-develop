#include "module_relax/relax_old/ions_move_basic.h"
#define private public
#define protected public
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_relax/relax_old/bfgs_basic.h"
/************************************************
 *  unit tests of class BFGS_Basic
 ***********************************************/

/**
 * - Tested Functions:
 *   - BFGS_Basic::allocate_basic()
 *   - BFGS_Basic::new_step()
 *   - BFGS_Basic::reset_hessian()
 *   - BFGS_Basic::save_bfgs()
 *   - BFGS_Basic::check_move()
 *   - BFGS_Basic::update_inverse_hessian()
 *   - BFGS_Basic::check_wolfe_conditions()
 *   - BFGS_Basic::compute_trust_radius()
 */

int Ions_Move_Basic::dim = 0;
bool Ions_Move_Basic::converged = false;
double Ions_Move_Basic::largest_grad = 0.0;
int Ions_Move_Basic::update_iter = 0;
int Ions_Move_Basic::istep = 0;
double Ions_Move_Basic::ediff = 0.0;
double Ions_Move_Basic::etot = 0.0;
double Ions_Move_Basic::etot_p = 0.0;
double Ions_Move_Basic::trust_radius = 0.0;
double Ions_Move_Basic::trust_radius_old = 0.0;
double Ions_Move_Basic::relax_bfgs_rmax = -1.0;
double Ions_Move_Basic::relax_bfgs_rmin = -1.0;
double Ions_Move_Basic::relax_bfgs_init = -1.0;
double Ions_Move_Basic::best_xxx = 1.0;
int Ions_Move_Basic::out_stru = 0;

double Ions_Move_Basic::dot_func(const double *a, const double *b, const int &dim_in)
{
    double result = 0.0;
    for (int i = 0; i < dim_in; i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

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

    // Check if allocated arrays are not empty
    EXPECT_NE(nullptr, bfgs.pos);
    EXPECT_NE(nullptr, bfgs.pos_p);
    EXPECT_NE(nullptr, bfgs.grad);
    EXPECT_NE(nullptr, bfgs.grad_p);
    EXPECT_NE(nullptr, bfgs.move);
    EXPECT_NE(nullptr, bfgs.move_p);
    EXPECT_NE(nullptr, bfgs.inv_hess.c);
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
    ASSERT_DEATH(bfgs.update_inverse_hessian(lat0), "");
}

// Test function update_inverse_hessian() when sdoty = 0
TEST_F(BFGSBasicTest, UpdateInverseHessianCase1)
{
    Ions_Move_Basic::dim = 3;
    double lat0 = 1.0;
    bfgs.allocate_basic();

    GlobalV::ofs_running.open("log");
    bfgs.update_inverse_hessian(lat0);
    GlobalV::ofs_running.close();

    std::string expected_output
        = " WARINIG: unexpected behaviour in update_inverse_hessian\n Resetting bfgs history \n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    EXPECT_EQ(expected_output, output);
    std::remove("log");
}

// Test function update_inverse_hessian()
TEST_F(BFGSBasicTest, UpdateInverseHessianCase2)
{
    Ions_Move_Basic::dim = 3;
    double lat0 = 1.0;
    bfgs.allocate_basic();
    bfgs.pos[0] = 2.0;
    bfgs.grad[0] = 2.0;

    bfgs.update_inverse_hessian(lat0);

    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 0), 0.5);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(1, 2), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(2, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.inv_hess(2, 2), 0.0);
}

// Test function check_wolfe_conditions()
TEST_F(BFGSBasicTest, CheckWolfeConditions)
{
    Ions_Move_Basic::dim = 3;
    Ions_Move_Basic::etot = 10.0;
    GlobalV::test_relax_method = 1;
    bfgs.allocate_basic();
    bfgs.pos[0] = 2.0;
    bfgs.grad[0] = 2.0;
    bfgs.move[0] = 1.0;

    GlobalV::ofs_running.open("log");
    bfgs.check_wolfe_conditions();
    GlobalV::ofs_running.close();

    std::string expected_output
        = "                            etot - etot_p = 10\n                    relax_bfgs_w1 * dot_p = -0\n            "
          "                          dot = 0\n                    relax_bfgs_w2 * dot_p = -0\n                         "
          "   relax_bfgs_w1 = -1\n                            relax_bfgs_w2 = -1\n                                   "
          "wolfe1 = 0\n                                   wolfe2 = 0\n                            etot - etot_p = 10\n "
          "                   relax_bfgs_w1 * dot_p = -0\n                                   wolfe1 = 0\n              "
          "                     wolfe2 = 0\n                wolfe condition satisfied = 0\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    EXPECT_EQ(bfgs.wolfe_flag, false);
    EXPECT_EQ(expected_output, output);
    std::remove("log");
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
    Ions_Move_Basic::update_iter = 0;
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
    bfgs.new_step(lat0);

    EXPECT_EQ(Ions_Move_Basic::update_iter, 1);
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
    Ions_Move_Basic::update_iter = 2;
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
    bfgs.new_step(lat0);

    EXPECT_EQ(Ions_Move_Basic::update_iter, 3);
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
    bfgs.bfgs_ndim = 2;
    bfgs.allocate_basic();
    double lat0 = 1.0;

    testing::internal::CaptureStdout();
    EXPECT_EXIT(bfgs.new_step(lat0), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("bfgs_ndim > 1 not implemented yet"));
}

// Test function compute_trust_radius() case 1
TEST_F(BFGSBasicTest, ComputeTrustRadiusCase1)
{
    Ions_Move_Basic::dim = 2;
    Ions_Move_Basic::etot = 0.0;
    Ions_Move_Basic::etot_p = 0.0;
    bfgs.allocate_basic();
    bfgs.grad_p[0] = 1.0;
    bfgs.move_p[1] = 2.0;
    bfgs.inv_hess(0, 0) = -3.0;
    bfgs.inv_hess(0, 1) = -4.0;
    bfgs.inv_hess(1, 0) = -5.0;
    bfgs.inv_hess(1, 1) = -6.0;
    bfgs.wolfe_flag = true;
    bfgs.relax_bfgs_w1 = 1.0;

    bfgs.compute_trust_radius();

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
    Ions_Move_Basic::etot = 0.0;
    Ions_Move_Basic::etot_p = 0.0;
    Ions_Move_Basic::trust_radius_old = 0.0;
    Ions_Move_Basic::relax_bfgs_rmin = 100.0;
    GlobalV::test_relax_method = 1;
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

    bfgs.compute_trust_radius();

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
    Ions_Move_Basic::etot = 0.0;
    Ions_Move_Basic::etot_p = 0.0;
    Ions_Move_Basic::trust_radius_old = 0.0;
    Ions_Move_Basic::relax_bfgs_rmin = 100.0;
    GlobalV::test_relax_method = 1;
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

    testing::internal::CaptureStdout();
    EXPECT_EXIT(bfgs.compute_trust_radius(), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("bfgs history already reset at previous step, we got trapped!"));
}