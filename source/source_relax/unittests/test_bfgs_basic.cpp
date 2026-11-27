#include "source_io/module_parameter/parameter.h"
#include "source_relax/bfgs_basic.h"
#include "source_relax/ions_move_basic.h"
#include "source_relax/relax_data.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
/************************************************
 *  unit tests of class BFGS_Basic
 ***********************************************/

class BFGSBasicTest : public ::testing::Test
{
  public:
    int test_relax_method = 0;

  protected:
    void SetUp() override
    {
        Ions_Move_Basic::trust_radius_old = 0.1;
        Ions_Move_Basic::relax_bfgs_rmax = 0.2;
        Ions_Move_Basic::relax_bfgs_rmin = 1.e-5;
        BFGS_Basic::relax_bfgs_w1 = 0.01;
        BFGS_Basic::relax_bfgs_w2 = 0.5;
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
    bfgs.allocate_basic_for_testing();

    // Check if allocated vectors are not empty
    EXPECT_EQ(bfgs.get_pos().size(), 4U);
    EXPECT_EQ(bfgs.get_pos_p().size(), 4U);
    EXPECT_EQ(bfgs.get_grad().size(), 4U);
    EXPECT_EQ(bfgs.get_grad_p().size(), 4U);
    EXPECT_EQ(bfgs.get_move().size(), 4U);
    EXPECT_EQ(bfgs.get_move_p().size(), 4U);
}

// Test if a dimension less than or equal to 0 results in an assertion error
TEST_F(BFGSBasicTest, TestAllocateWithZeroDimension)
{
    Ions_Move_Basic::dim = 0;
    ASSERT_DEATH(bfgs.allocate_basic_for_testing(), "");
}

// Test function update_inverse_hessian() assert death
TEST_F(BFGSBasicTest, UpdateInverseHessianDeath)
{
    Ions_Move_Basic::dim = 0;
    std::ofstream ofs("test_log_update_inverse_hessian_death.log");
    ASSERT_DEATH(bfgs.update_inverse_hessian_for_testing(ofs), "");
    ofs.close();
    std::remove("test_log_update_inverse_hessian_death.log");
}

// Test function update_inverse_hessian() when sdoty = 0
TEST_F(BFGSBasicTest, UpdateInverseHessianCase1)
{
    Ions_Move_Basic::dim = 3;
    bfgs.allocate_basic_for_testing();

    std::ofstream ofs("test_log_update_inverse_hessian_case1.log");
    bfgs.update_inverse_hessian_for_testing(ofs);
    ofs.close();

    std::string expected_output = " Skipping inverse BFGS update with insufficient positive curvature.\n";
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
    bfgs.allocate_basic_for_testing();
    bfgs.get_move_p()[0] = 1.0;
    bfgs.get_grad()[0] = 2.0;

    std::ofstream ofs("test_log_update_inverse_hessian_case2.log");
    bfgs.update_inverse_hessian_for_testing(ofs);
    ofs.close();

    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(0, 0), 0.5);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(1, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(1, 2), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(2, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(2, 2), 0.0);
    std::remove("test_log_update_inverse_hessian_case2.log");
}

// Test function check_wolfe_conditions()
TEST_F(BFGSBasicTest, CheckWolfeConditions)
{
    Ions_Move_Basic::dim = 1;
    bfgs.allocate_basic_for_testing();
    bfgs.relax_bfgs_w1 = 0.25;
    bfgs.get_grad_p()[0] = -2.0;
    bfgs.get_move_p()[0] = 0.25;
    std::vector<double> energies = {0.875, 1.0};
    std::ofstream ofs;

    // Both inequalities are satisfied exactly at the boundary.
    bfgs.get_grad()[0] = 1.0;
    bfgs.check_wolfe_conditions_for_testing(ofs, energies);
    EXPECT_TRUE(bfgs.get_wolfe_flag());

    bfgs.get_grad()[0] = -1.5;
    bfgs.check_wolfe_conditions_for_testing(ofs, energies);
    EXPECT_FALSE(bfgs.get_wolfe_flag());

    bfgs.get_grad()[0] = 0.0;
    energies[0] = 0.9;
    bfgs.check_wolfe_conditions_for_testing(ofs, energies);
    EXPECT_FALSE(bfgs.get_wolfe_flag());
}

// Test function reset_hessian()
TEST_F(BFGSBasicTest, ResetHessian)
{
    Ions_Move_Basic::dim = 3;
    bfgs.allocate_basic_for_testing();

    bfgs.reset_hessian_for_testing();

    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(1, 2), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(2, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(2, 2), 1.0);
}

// Test function save_bfgs()
TEST_F(BFGSBasicTest, SaveBfgs)
{
    bfgs.get_save_flag() = false;
    Ions_Move_Basic::dim = 2;
    bfgs.allocate_basic_for_testing();
    bfgs.get_pos()[0] = 1.0;
    bfgs.get_pos()[1] = 2.0;
    bfgs.get_grad()[0] = 3.0;
    bfgs.get_grad()[1] = 4.0;
    bfgs.get_move()[0] = 5.0;
    bfgs.get_move()[1] = 6.0;

    bfgs.save_bfgs_for_testing();

    EXPECT_EQ(bfgs.get_save_flag(), true);
    EXPECT_DOUBLE_EQ(bfgs.get_pos_p()[0], 1.0);
    EXPECT_DOUBLE_EQ(bfgs.get_pos_p()[1], 2.0);
    EXPECT_DOUBLE_EQ(bfgs.get_grad_p()[0], 3.0);
    EXPECT_DOUBLE_EQ(bfgs.get_grad_p()[1], 4.0);
    EXPECT_DOUBLE_EQ(bfgs.get_move()[0], 5.0);
    EXPECT_DOUBLE_EQ(bfgs.get_move()[1], 6.0);
}

// Test function new_step() when update_iter == 1
TEST_F(BFGSBasicTest, NewStepCase1)
{
    int update_iter = 0;
    Ions_Move_Basic::largest_grad = 0.0;
    Ions_Move_Basic::relax_bfgs_init = 0.3;
    Ions_Move_Basic::best_xxx = -0.4;
    bfgs.get_bfgs_ndim() = 1;
    Ions_Move_Basic::dim = 2;
    bfgs.allocate_basic_for_testing();
    bfgs.get_grad()[0] = 1.0;
    bfgs.get_grad()[1] = 2.0;
    bfgs.get_inv_hess()(0, 0) = -3.0;
    bfgs.get_inv_hess()(0, 1) = -4.0;
    bfgs.get_inv_hess()(1, 0) = -5.0;
    bfgs.get_inv_hess()(1, 1) = -6.0;

    std::ofstream ofs("test_log.log");
    std::vector<double> etot_info(2, 0.0);
    bfgs.new_step_for_testing(update_iter, ofs, etot_info, test_relax_method);

    EXPECT_EQ(update_iter, 1);
    EXPECT_EQ(bfgs.get_tr_min_hit(), false);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::relax_bfgs_init, 0.3);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::best_xxx, -0.4);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, 0.2);
    EXPECT_DOUBLE_EQ(bfgs.get_move()[0], -1.0);
    EXPECT_DOUBLE_EQ(bfgs.get_move()[1], -2.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(1, 1), 1.0);
}

// Test function new_step() when update_iter > 1
TEST_F(BFGSBasicTest, NewStepCase2)
{
    int update_iter = 2;
    Ions_Move_Basic::largest_grad = 0.0;
    Ions_Move_Basic::relax_bfgs_init = 0.3;
    Ions_Move_Basic::best_xxx = -0.4;
    bfgs.get_bfgs_ndim() = 1;
    Ions_Move_Basic::dim = 2;
    bfgs.allocate_basic_for_testing();
    bfgs.get_grad()[0] = 1.0;
    bfgs.get_grad()[1] = 2.0;
    bfgs.get_inv_hess()(0, 0) = -3.0;
    bfgs.get_inv_hess()(0, 1) = -4.0;
    bfgs.get_inv_hess()(1, 0) = -5.0;
    bfgs.get_inv_hess()(1, 1) = -6.0;

    const std::string log_path = testing::TempDir() + "bfgs_basic_zero_curvature.log";
    std::ofstream ofs(log_path);
    std::vector<double> etot_info(2, 0.0);
    bfgs.new_step_for_testing(update_iter, ofs, etot_info, test_relax_method);

    ofs.close();
    std::ifstream log(log_path);
    const std::string output((std::istreambuf_iterator<char>(log)), std::istreambuf_iterator<char>());
    log.close();
    std::remove(log_path.c_str());
    EXPECT_THAT(output, testing::HasSubstr("Skipping inverse BFGS update with insufficient positive curvature"));
    EXPECT_EQ(update_iter, 3);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, 0.1);
    EXPECT_DOUBLE_EQ(bfgs.get_move()[0], -1.0);
    EXPECT_DOUBLE_EQ(bfgs.get_move()[1], -2.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(1, 1), 1.0);
}

// Test function new_step() when bfgs_ndim > 1
TEST_F(BFGSBasicTest, NewStepWarningQuit)
{
    int update_iter = 0;
    bfgs.get_bfgs_ndim() = 2;
    Ions_Move_Basic::dim = 2;
    bfgs.allocate_basic_for_testing();
    std::ofstream ofs("test_log.log");
    std::vector<double> etot_info(2, 0.0);

    testing::internal::CaptureStdout();
    EXPECT_EXIT(
        bfgs.new_step_for_testing(update_iter, ofs, etot_info, test_relax_method), ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("bfgs_ndim > 1 not implemented yet"));
}

TEST_F(BFGSBasicTest, AcceptedInteriorStepKeepsRadius)
{
    Ions_Move_Basic::dim = 1;
    bfgs.allocate_basic_for_testing();
    bfgs.get_grad_p()[0] = -1.0;
    bfgs.get_move_p()[0] = 0.01;
    bfgs.get_wolfe_flag() = true;
    std::vector<double> energies = {-0.01, 0.0};
    std::ofstream ofs;
    bfgs.compute_trust_radius_for_testing(ofs, energies, test_relax_method);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, 0.1);
}

TEST_F(BFGSBasicTest, AcceptedBoundaryStepCanGrowRadius)
{
    Ions_Move_Basic::dim = 1;
    bfgs.allocate_basic_for_testing();
    bfgs.get_grad_p()[0] = -1.0;
    bfgs.get_move_p()[0] = 0.1;
    std::vector<double> energies = {-0.05, 0.0};
    std::ofstream ofs;
    bfgs.compute_trust_radius_for_testing(ofs, energies, test_relax_method);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, 0.15);

    bfgs.get_wolfe_flag() = true;
    Ions_Move_Basic::relax_bfgs_rmax = 0.18;
    bfgs.compute_trust_radius_for_testing(ofs, energies, test_relax_method);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, 0.18);
}
