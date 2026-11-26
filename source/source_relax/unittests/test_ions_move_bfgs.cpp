#include "for_test.h"
#include "source_io/module_parameter/parameter.h"
#include "source_relax/ions_move_basic.h"
#include "source_relax/ions_move_bfgs.h"
#include "source_relax/relax_criteria.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>

class IonicBfgsStepTest : public ::testing::TestWithParam<double>
{
  protected:
    UnitCell cell;
    Ions_Move_BFGS optimizer;
    Relax_Criteria criteria;
    ModuleBase::matrix force{2, 3};
    std::vector<double> energies = {0.0, 0.0};
    std::ofstream log;
    int step = 0;
    int iteration = 0;

    void SetUp() override
    {
        // Two equivalent descriptions of an 8-Bohr cubic cell.
        cell.lat0 = GetParam();
        cell.latvec.Identity();
        cell.latvec *= 8.0 / cell.lat0;
        cell.GT = cell.latvec.Inverse();
        cell.atoms[0].tau[0] = {};
        cell.atoms[0].tau[1] = {};
        Ions_Move_Basic::dim = 6;
        Ions_Move_Basic::relax_bfgs_init = 0.5;
        Ions_Move_Basic::relax_bfgs_rmax = 0.2;
        Ions_Move_Basic::relax_bfgs_rmin = 1.e-5;
        Ions_Move_Basic::best_xxx = 1.0;
        BFGS_Basic::relax_bfgs_w1 = 0.01;
        BFGS_Basic::relax_bfgs_w2 = 0.5;
        criteria.force_thr = 1.e-8;
        criteria.out_level = "m";
        optimizer.allocate();
    }

    void TearDown() override { delete[] cell.atoms; }

    double x() const { return cell.atoms[0].tau[0].x * cell.lat0; }

    bool evaluate(double energy, double fx)
    {
        force(0, 0) = fx;
        ++step;
        return optimizer.start(cell, force, energy, step, iteration, log, energies, criteria);
    }
};

TEST_P(IonicBfgsStepTest, SoftHarmonicModeCrossesBoundary)
{
    // E = 0.05 (x - 8.076)^2 Ry, x in Bohr. The two Newton steps fit inside rmax.
    cell.atoms[0].tau[0].x = 7.996 / cell.lat0;
    Ions_Move_Basic::relax_bfgs_rmin = 0.01;
    const auto harmonic = [&]()
    {
        const double dx = std::remainder(x() - 0.076, 8.0);
        return evaluate(0.05 * dx * dx, -0.1 * dx);
    };
    ASSERT_FALSE(harmonic());
    EXPECT_NEAR(x(), 0.004, 1.e-12);
    ASSERT_FALSE(harmonic());
    EXPECT_NEAR(x(), 0.076, 1.e-12);
    EXPECT_NEAR(8.0 * cell.atoms[0].dis[0].x, 0.072, 1.e-12);
    EXPECT_EQ(cell.atoms[0].boundary_shift[0].x, 0);
    EXPECT_TRUE(harmonic());
}

TEST_P(IonicBfgsStepTest, RetryUsesWorkAlongAppliedStep)
{
    // E = 4 (x - 1.05)^2 Ry. The first step hits rmax, interpolation gives the minimum.
    cell.atoms[0].tau[0].x = 1.0 / cell.lat0;
    const auto harmonic = [&]()
    {
        const double dx = x() - 1.05;
        return evaluate(4.0 * dx * dx, -8.0 * dx);
    };
    ASSERT_FALSE(harmonic());
    EXPECT_NEAR(x(), 1.2, 1.e-12);
    ASSERT_FALSE(harmonic());
    EXPECT_NEAR(x(), 1.05, 1.e-12);
    EXPECT_TRUE(harmonic());
}

TEST_P(IonicBfgsStepTest, NegativeCurvatureKeepsADescentStep)
{
    // E = -x^2/4 has negative curvature along the step.
    cell.atoms[0].tau[0].x = 0.1 / cell.lat0;
    ASSERT_FALSE(evaluate(-0.0025, 0.05));
    EXPECT_NEAR(x(), 0.15, 1.e-12);
    ASSERT_FALSE(evaluate(-0.005625, 0.075));
    EXPECT_NEAR(x(), 0.3, 1.e-12);
}

TEST_P(IonicBfgsStepTest, RejectedTrialsKeepAcceptedPoint)
{
    cell.atoms[0].tau[0].x = 7.96 / cell.lat0;
    ASSERT_FALSE(evaluate(0.0, 0.12));
    EXPECT_NEAR(x(), 0.08, 1.e-12);
    // A flat, higher-energy trial still needs to be rejected.
    ASSERT_FALSE(evaluate(0.0144, 0.0));
    EXPECT_NEAR(x(), 7.99, 1.e-12);
    ASSERT_FALSE(evaluate(0.0036, -0.5));
    EXPECT_NEAR(x(), 7.9675, 1.e-12);
    EXPECT_EQ(iteration, 1);
}

INSTANTIATE_TEST_SUITE_P(LatticeConstant, IonicBfgsStepTest, ::testing::Values(1.0, 8.0));

/************************************************
 *  unit tests of class Ions_Move_BFGS
 ***********************************************/

// Define a fixture for the tests
class IonsMoveBFGSTest : public ::testing::Test
{
  public:
    Relax_Criteria criteria;

  protected:
    Ions_Move_BFGS bfgs;
    int update_iter;

    virtual void SetUp()
    {
        // Initialize variables before each test
        Ions_Move_Basic::dim = 6;
        Ions_Move_Basic::trust_radius = 0.2;
        Ions_Move_Basic::relax_bfgs_init = 0.5;
        Ions_Move_Basic::relax_bfgs_rmax = 0.2;
        Ions_Move_Basic::relax_bfgs_rmin = 1.e-5;
        Ions_Move_Basic::best_xxx = 1.0;
        update_iter = 0;
    }

    virtual void TearDown()
    {
        // Clean up after each test
    }
};

// Test the allocate() function case 1
TEST_F(IonsMoveBFGSTest, AllocateCase1)
{
    // Initilize data
    bfgs.get_init_done() = true;
    bfgs.get_save_flag() = true;

    // Call the function being tested
    bfgs.allocate();

    // Check that the expected results
    EXPECT_EQ(bfgs.get_init_done(), true);
    EXPECT_EQ(bfgs.get_save_flag(), true);
}

// Test the allocate() function case 2
TEST_F(IonsMoveBFGSTest, AllocateCase2)
{
    // Initilize data
    bfgs.get_init_done() = false;
    bfgs.get_save_flag() = true;

    // Call the function being tested
    bfgs.allocate();

    // Check that the expected results
    EXPECT_EQ(bfgs.get_init_done(), true);
    EXPECT_EQ(bfgs.get_save_flag(), false);
}

// Test the start() function case 1
TEST_F(IonsMoveBFGSTest, StartCase1)
{
    // Initilize data
    UnitCell ucell;
    ModuleBase::matrix force(2, 3);
    double energy_in = 0.0;
    const int istep = 1;
    bfgs.get_init_done() = false;
    bfgs.get_save_flag() = true;
    std::vector<double> etot_info(2, 0.0);

    // Call the function being tested
    bfgs.allocate();
    std::ofstream ofs("test_start_case1.log");
    bfgs.start(ucell, force, energy_in, istep, update_iter, ofs, etot_info, criteria);
    ofs.close();

    // Check the results
    std::ifstream ifs("test_start_case1.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("test_start_case1.log");

    EXPECT_THAT(output, testing::HasSubstr("update iteration"));
}

// Test the start() function case 2
TEST_F(IonsMoveBFGSTest, StartCase2)
{
    // Initilize data
    UnitCell ucell;
    // Initialize UnitCell with 2 atoms
    ucell.ntype = 1;
    ucell.nat = 2;
    ucell.atoms = new Atom[ucell.ntype];
    ucell.atoms[0].na = 2;
    ucell.atoms[0].tau = std::vector<ModuleBase::Vector3<double>>(2);
    ucell.atoms[0].taud = std::vector<ModuleBase::Vector3<double>>(2);
    ucell.atoms[0].mbl = std::vector<ModuleBase::Vector3<int>>(2, {1, 1, 1});
    ucell.atoms[0].tau[0].x = 0.0;
    ucell.atoms[0].tau[0].y = 0.0;
    ucell.atoms[0].tau[0].z = 0.0;
    ucell.atoms[0].tau[1].x = 1.0;
    ucell.atoms[0].tau[1].y = 0.0;
    ucell.atoms[0].tau[1].z = 0.0;
    ucell.lat0 = 1.0;
    ucell.set_atom_flag = true;

    // Initialize PARAM
    criteria.force_thr = 1.0e-3;
    criteria.force_thr_ev = criteria.force_thr * 13.6058 / 0.529177;
    criteria.test_relax_method = 1;
    criteria.out_level = "ie";

    // Initialize istep
    const int istep = 1;

    ModuleBase::matrix force(2, 3);
    force(0, 0) = 10.0;
    double energy_in = 0.0;
    bfgs.get_init_done() = false;
    bfgs.get_save_flag() = true;
    std::vector<double> etot_info(2, 0.0);

    // Call the function being tested
    bfgs.allocate();
    std::ofstream ofs("test_start_case2.log");
    bfgs.start(ucell, force, energy_in, istep, update_iter, ofs, etot_info, criteria);
    ofs.close();

    // Check the results
    std::ifstream ifs("test_start_case2.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("test_start_case2.log");

    // Reporting moved to IonCellOptimizer::relax_step; check_converged no longer
    // prints convergence info to the log. Assert state instead.
    EXPECT_THAT(output, testing::Not(testing::HasSubstr("Ion relaxation is not converged yet")));
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 10.0);

    // Clean up
    delete[] ucell.atoms;
}

// Test the restart_bfgs() function case 1
TEST_F(IonsMoveBFGSTest, RestartBfgsCase1)
{
    // Initilize data
    bfgs.get_init_done() = false;
    criteria.test_relax_method = 1;
    bfgs.allocate();
    bfgs.get_save_flag() = true;
    bfgs.get_move_p()[0] = 0.05;

    std::ofstream ofs;
    bfgs.restart_bfgs_for_testing(update_iter, ofs, criteria.test_relax_method);

    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius_old, 0.2);
    EXPECT_NEAR(bfgs.get_move_p()[0], 0.05, 1.e-12);
}

// Test the restart_bfgs() function case 2
TEST_F(IonsMoveBFGSTest, RestartBfgsCase2)
{
    // Initilize data
    bfgs.get_init_done() = false;
    bfgs.allocate();
    criteria.test_relax_method = 1;
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        bfgs.get_move_p()[i] = 1.0;
        bfgs.get_pos()[i] = i;
        bfgs.get_pos_p()[i] = i;
    }

    // Call the function being tested
    std::ofstream ofs("test_restart_bfgs_case2.log");
    bfgs.restart_bfgs_for_testing(update_iter, ofs, criteria.test_relax_method);
    ofs.close();
    std::remove("test_restart_bfgs_case2.log");

    // Check the results
    EXPECT_DOUBLE_EQ(update_iter, 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_tr_min_hit(), false);
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        EXPECT_DOUBLE_EQ(bfgs.get_pos_p()[i], 0.0);
        EXPECT_DOUBLE_EQ(bfgs.get_grad_p()[i], 0.0);
        EXPECT_DOUBLE_EQ(bfgs.get_move_p()[i], 0.0);
        for (int j = 0; j < Ions_Move_Basic::dim; ++j)
        {
            if (i == j)
            {
                EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(i, j), 1.0);
            }
            else
            {
                EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(i, j), 0.0);
            }
        }
    }
}

// Test the bfgs_routine() function case 1
TEST_F(IonsMoveBFGSTest, BfgsRoutineCase1)
{
    // Initilize data
    bfgs.get_init_done() = false;
    bfgs.allocate();
    bfgs.get_tr_min_hit() = false;
    criteria.test_relax_method = 1;
    criteria.out_level = "ie";
    const int istep = 1;
    std::vector<double> etot_info = {1.0, 0.9, 0.1};
    Ions_Move_Basic::relax_bfgs_rmin = 0.1;
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        bfgs.get_move_p()[i] = 0.0;
        bfgs.get_grad_p()[i] = i;
        bfgs.get_pos_p()[i] = i;
        bfgs.get_pos()[i] = i + bfgs.get_move_p()[i];
    }

    // Call the function being tested
    std::ofstream ofs("test_bfgs_routine_case1.log");
    testing::internal::CaptureStdout();
    bfgs.bfgs_routine_for_testing(istep, update_iter, ofs, etot_info, criteria.out_level, criteria.test_relax_method);
    std::string std_outout = testing::internal::GetCapturedStdout();
    ofs.close();

    // Check the results
    std::ifstream ifs("test_bfgs_routine_case1.log");
    std::string ofs_output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("test_bfgs_routine_case1.log");

    EXPECT_THAT(ofs_output, ::testing::HasSubstr("reset bfgs history"));
    EXPECT_EQ(" BFGS TRUST (Bohr)    : 0.1\n", std_outout);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, 0.1);
    EXPECT_DOUBLE_EQ(etot_info[0], 0.9);
    EXPECT_DOUBLE_EQ(bfgs.get_tr_min_hit(), true);
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        EXPECT_NEAR(bfgs.get_move()[i], -0.1 * i / std::sqrt(55.0), 1.e-12);
        EXPECT_DOUBLE_EQ(bfgs.get_pos()[i], i);
        EXPECT_DOUBLE_EQ(bfgs.get_grad()[i], i);
        for (int j = 0; j < Ions_Move_Basic::dim; ++j)
        {
            if (i == j)
            {
                EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(i, j), 1.0);
            }
            else
            {
                EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(i, j), 0.0);
            }
        }
    }
}

// Test the bfgs_routine() function case 2
TEST_F(IonsMoveBFGSTest, BfgsRoutineCase2)
{
    // Initilize data
    bfgs.get_init_done() = false;
    bfgs.allocate();
    bfgs.get_tr_min_hit() = false;
    criteria.test_relax_method = 0;
    criteria.out_level = "none";
    const int istep = 1;
    std::vector<double> etot_info = {1.0, 0.9, 0.1};
    Ions_Move_Basic::relax_bfgs_rmax = 10.0;
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        bfgs.get_move_p()[i] = i;
        bfgs.get_grad_p()[i] = i;
        bfgs.get_pos_p()[i] = i;
        bfgs.get_pos()[i] = i + bfgs.get_move_p()[i];
    }

    // Call the function being tested
    std::ofstream ofs("test_bfgs_routine_case2.log");
    testing::internal::CaptureStdout();
    bfgs.bfgs_routine_for_testing(istep, update_iter, ofs, etot_info, criteria.out_level, criteria.test_relax_method);
    std::string std_outout = testing::internal::GetCapturedStdout();
    ofs.close();

    // Check the results
    std::ifstream ifs("test_bfgs_routine_case2.log");
    std::string ofs_output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("test_bfgs_routine_case2.log");

    std::string expected_ofs = " quadratic interpolation is impossible.\n";
    std::string expected_std = "";

    EXPECT_THAT(ofs_output, ::testing::HasSubstr(expected_ofs));
    EXPECT_EQ(expected_std, std_outout);

    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, 0.5 * std::sqrt(55.0));
    EXPECT_DOUBLE_EQ(etot_info[0], 0.9);
    EXPECT_FALSE(bfgs.get_tr_min_hit());
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        EXPECT_DOUBLE_EQ(bfgs.get_move()[i], 0.5 * i);
        EXPECT_DOUBLE_EQ(bfgs.get_pos()[i], i);
        EXPECT_DOUBLE_EQ(bfgs.get_grad()[i], i);
        for (int j = 0; j < Ions_Move_Basic::dim; ++j)
        {
            EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(i, j), 0.0);
        }
    }
}

// Test the bfgs_routine() function case 3
TEST_F(IonsMoveBFGSTest, BfgsRoutineCase3)
{
    // Initilize data
    const int istep = 1;
    std::vector<double> etot_info = {0.9, 1.0, -0.1};
    update_iter = 0;
    Ions_Move_Basic::largest_grad = 0.0;
    Ions_Move_Basic::relax_bfgs_init = 0.3;
    Ions_Move_Basic::best_xxx = -0.4;
    bfgs.get_init_done() = false;
    bfgs.allocate();
    bfgs.get_bfgs_ndim() = 1;
    bfgs.get_grad()[0] = 1.0;
    bfgs.get_grad()[1] = 2.0;
    bfgs.get_inv_hess()(0, 0) = -3.0;
    bfgs.get_inv_hess()(0, 1) = -4.0;
    bfgs.get_inv_hess()(1, 0) = -5.0;
    bfgs.get_inv_hess()(1, 1) = -6.0;

    // Call the function being tested
    std::ofstream ofs("test_bfgs_routine_case3.log");
    bfgs.bfgs_routine_for_testing(istep, update_iter, ofs, etot_info, criteria.out_level, criteria.test_relax_method);
    ofs.close();

    // Check the results
    std::ifstream ifs("test_bfgs_routine_case3.log");
    std::string ofs_output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("test_bfgs_routine_case3.log");

    std::string expected_ofs = " check the norm of new move 20.2485 (Bohr)\n Uphill move : resetting bfgs history\n";

    EXPECT_THAT(ofs_output, ::testing::HasSubstr(expected_ofs));
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, 0.2);
    EXPECT_DOUBLE_EQ(etot_info[0], 0.9);
    EXPECT_DOUBLE_EQ(bfgs.get_tr_min_hit(), false);
    EXPECT_NEAR(bfgs.get_move()[0], -0.089442719099991588, 1e-12);
    EXPECT_NEAR(bfgs.get_move()[1], -0.17888543819998318, 1e-12);
    EXPECT_DOUBLE_EQ(bfgs.get_move()[2], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_move()[3], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_move()[4], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_move()[5], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_pos()[0], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_pos()[1], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_pos()[2], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_pos()[3], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_pos()[4], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_pos()[5], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_grad()[0], 1.0);
    EXPECT_DOUBLE_EQ(bfgs.get_grad()[1], 2.0);
    EXPECT_DOUBLE_EQ(bfgs.get_grad()[2], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_grad()[3], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_grad()[4], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.get_grad()[5], 0.0);
}

// Test the bfgs_routine() function warning quit 1
TEST_F(IonsMoveBFGSTest, BfgsRoutineWarningQuit1)
{
    // Initilize data
    bfgs.get_init_done() = false;
    bfgs.allocate();
    bfgs.get_tr_min_hit() = true;
    criteria.test_relax_method = 1;
    criteria.out_level = "ie";
    const int istep = 1;
    std::vector<double> etot_info = {1.0, 0.9, 0.1};
    Ions_Move_Basic::relax_bfgs_rmin = 0.1;
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        bfgs.get_move_p()[i] = 0.0;
        bfgs.get_grad_p()[i] = i;
        bfgs.get_pos_p()[i] = i;
        bfgs.get_pos()[i] = i + bfgs.get_move_p()[i];
    }

    // Check the results
    std::ofstream ofs("test_bfgs_routine_warning_quit1.log");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(
        bfgs.bfgs_routine_for_testing(
            istep, update_iter, ofs, etot_info, criteria.out_level, criteria.test_relax_method),
        ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    ofs.close();
    std::remove("test_bfgs_routine_warning_quit1.log");
    EXPECT_THAT(output, testing::HasSubstr("trust radius is too small! Break down."));
}

// Test the bfgs_routine() function warning quit 2
TEST_F(IonsMoveBFGSTest, BfgsRoutineWarningQuit2)
{
    // Initilize data
    bfgs.get_init_done() = false;
    bfgs.allocate();
    bfgs.get_tr_min_hit() = false;
    criteria.test_relax_method = 1;
    criteria.out_level = "ie";
    const int istep = 1;
    std::vector<double> etot_info = {1.0, 0.9, 0.1};
    Ions_Move_Basic::relax_bfgs_rmin = 0.1;

    // Check the results
    std::ofstream ofs("test_bfgs_routine_warning_quit2.log");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(
        bfgs.bfgs_routine_for_testing(
            istep, update_iter, ofs, etot_info, criteria.out_level, criteria.test_relax_method),
        ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    ofs.close();
    std::remove("test_bfgs_routine_warning_quit2.log");
    EXPECT_THAT(output, testing::HasSubstr("BFGS: move-length unreasonably short"));
}
