#include "for_test.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#define private public
#define protected public
#include "source_io/module_parameter/parameter.h"
#include "source_relax/ions_move_basic.h"
#include "source_relax/ions_move_bfgs.h"
#undef private
#undef protected

/************************************************
 *  unit tests of class Ions_Move_BFGS
 ***********************************************/

// Define a fixture for the tests
class IonsMoveBFGSTest : public ::testing::Test
{
  protected:
    Ions_Move_BFGS bfgs;
    int update_iter;

    virtual void SetUp()
    {
        // Initialize variables before each test
        Ions_Move_Basic::dim = 6;
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
    bfgs.init_done = true;
    bfgs.save_flag = true;

    // Call the function being tested
    bfgs.allocate();

    // Check that the expected results
    EXPECT_EQ(bfgs.init_done, true);
    EXPECT_EQ(bfgs.save_flag, true);
}

// Test the allocate() function case 2
TEST_F(IonsMoveBFGSTest, AllocateCase2)
{
    // Initilize data
    bfgs.init_done = false;
    bfgs.save_flag = true;

    // Call the function being tested
    bfgs.allocate();

    // Check that the expected results
    EXPECT_EQ(bfgs.init_done, true);
    EXPECT_EQ(bfgs.save_flag, false);
}

// Test the start() function case 1
TEST_F(IonsMoveBFGSTest, StartCase1)
{
    // Initilize data
    UnitCell ucell;
    ModuleBase::matrix force(2, 3);
    double energy_in = 0.0;
    const int istep = 1;
    bfgs.init_done = false;
    bfgs.save_flag = true;
    std::vector<double> etot_info(2, 0.0);

    // Call the function being tested
    bfgs.allocate();
    std::ofstream ofs("test_start_case1.log");
    bfgs.start(ucell, force, energy_in, istep, update_iter, ofs, etot_info);
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
    ucell.atoms[0].tau[0].x = 0.0; ucell.atoms[0].tau[0].y = 0.0; ucell.atoms[0].tau[0].z = 0.0;
    ucell.atoms[0].tau[1].x = 1.0; ucell.atoms[0].tau[1].y = 0.0; ucell.atoms[0].tau[1].z = 0.0;
    ucell.lat0 = 1.0;
    ucell.set_atom_flag = true;

    // Initialize PARAM
    PARAM.input.force_thr = 1.0e-3;
    PARAM.input.force_thr_ev = PARAM.input.force_thr * 13.6058 / 0.529177;
    PARAM.input.test_relax_method = 1;
    PARAM.input.out_level = "ie";

    // Initialize istep
    const int istep = 1;

    ModuleBase::matrix force(2, 3);
    force(0, 0) = 10.0;
    double energy_in = 0.0;
    bfgs.init_done = false;
    bfgs.save_flag = true;
    std::vector<double> etot_info(2, 0.0);

    // Call the function being tested
    bfgs.allocate();
    std::ofstream ofs("test_start_case2.log");
    bfgs.start(ucell, force, energy_in, istep, update_iter, ofs, etot_info);
    ofs.close();

    // Check the results
    std::ifstream ifs("test_start_case2.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("test_start_case2.log");

    EXPECT_THAT(output, testing::HasSubstr("Ion relaxation is not converged yet"));

    // Clean up
    delete[] ucell.atoms;
}

// Test the restart_bfgs() function case 1
TEST_F(IonsMoveBFGSTest, RestartBfgsCase1)
{
    // Initilize data
    bfgs.init_done = false;
    PARAM.input.test_relax_method = 1;
    double lat0 = 1.0;
    bfgs.allocate();
    bfgs.save_flag = true;
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        bfgs.move_p[i] = 1.0;
        bfgs.pos[i] = 1.0;
        bfgs.pos_p[i] = i;
    }

    // Call the function being tested
    std::ofstream ofs("test_restart_bfgs_case1.log");
    bfgs.restart_bfgs(lat0, update_iter, ofs);
    ofs.close();

    // Check the results
    std::ifstream ifs("test_restart_bfgs_case1.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("test_restart_bfgs_case1.log");

    std::string expected_output = "                  trust_radius_old (bohr) = 2.44949\n";

    EXPECT_THAT(output, testing::HasSubstr(expected_output));
    EXPECT_NEAR(Ions_Move_Basic::trust_radius_old, 2.4494897427831779, 1e-12);
    EXPECT_DOUBLE_EQ(bfgs.move_p[0], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.move_p[1], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.move_p[2], 0.0);
    EXPECT_NEAR(bfgs.move_p[3], -0.40824829046386307, 1e-12);
    EXPECT_NEAR(bfgs.move_p[4], -0.81649658092772615, 1e-12);
    EXPECT_NEAR(bfgs.move_p[5], -1.2247448713915892, 1e-12);
}

// Test the restart_bfgs() function case 2
TEST_F(IonsMoveBFGSTest, RestartBfgsCase2)
{
    // Initilize data
    bfgs.init_done = false;
    bfgs.allocate();
    PARAM.input.test_relax_method = 1;
    double lat0 = 1.0;
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        bfgs.move_p[i] = 1.0;
        bfgs.pos[i] = i;
        bfgs.pos_p[i] = i;
    }

    // Call the function being tested
    std::ofstream ofs("test_restart_bfgs_case2.log");
    bfgs.restart_bfgs(lat0, update_iter, ofs);
    ofs.close();
    std::remove("test_restart_bfgs_case2.log");

    // Check the results
    EXPECT_DOUBLE_EQ(update_iter, 0.0);
    EXPECT_DOUBLE_EQ(bfgs.tr_min_hit, false);
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        EXPECT_DOUBLE_EQ(bfgs.pos_p[i], 0.0);
        EXPECT_DOUBLE_EQ(bfgs.grad_p[i], 0.0);
        EXPECT_DOUBLE_EQ(bfgs.move_p[i], 0.0);
        for (int j = 0; j < Ions_Move_Basic::dim; ++j)
        {
            if (i == j)
            {
                EXPECT_DOUBLE_EQ(bfgs.inv_hess(i, j), 1.0);
            }
            else
            {
                EXPECT_DOUBLE_EQ(bfgs.inv_hess(i, j), 0.0);
            }
        }
    }
}

// Test the bfgs_routine() function case 1
TEST_F(IonsMoveBFGSTest, BfgsRoutineCase1)
{
    // Initilize data
    bfgs.init_done = false;
    bfgs.allocate();
    bfgs.tr_min_hit = false;
    PARAM.input.test_relax_method = 1;
    PARAM.input.out_level = "ie";
    double lat0 = 1.0;
    const int istep = 1;
    std::vector<double> etot_info = {1.0, 0.9, 0.1};
    Ions_Move_Basic::relax_bfgs_rmin = 1.0;
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        bfgs.move_p[i] = 0.0;
        bfgs.grad_p[i] = i;
        bfgs.pos_p[i] = i;
    }

    // Call the function being tested
    std::ofstream ofs("test_bfgs_routine_case1.log");
    testing::internal::CaptureStdout();
    bfgs.bfgs_routine(lat0, istep, update_iter, ofs, etot_info);
    std::string std_outout = testing::internal::GetCapturedStdout();
    ofs.close();

    // Check the results
    std::ifstream ifs("test_bfgs_routine_case1.log");
    std::string ofs_output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("test_bfgs_routine_case1.log");

    std::string expected_ofs
        = "                                     dE0s = 0\n                                      den = 0.1\n            "
          "    interpolated trust radius = 0\ntrust_radius = 0\nrelax_bfgs_rmin = 1\nrelax_bfgs_rmax = -1\n "
          "trust_radius < relax_bfgs_rmin, reset bfgs history.\n";
    std::string expected_std = " BFGS TRUST (Bohr)    : 1\n";

    EXPECT_THAT(ofs_output, ::testing::HasSubstr(expected_ofs));
    EXPECT_EQ(expected_std, std_outout);

    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, 1.0);
    EXPECT_DOUBLE_EQ(etot_info[0], 0.9);
    EXPECT_DOUBLE_EQ(bfgs.tr_min_hit, true);
    EXPECT_NEAR(bfgs.move[0], 0.0, 1e-12);
    EXPECT_NEAR(bfgs.move[1], -0.13483997249264842, 1e-12);
    EXPECT_NEAR(bfgs.move[2], -0.26967994498529685, 1e-12);
    EXPECT_NEAR(bfgs.move[3], -0.40451991747794525, 1e-12);
    EXPECT_NEAR(bfgs.move[4], -0.5393598899705937, 1e-12);
    EXPECT_NEAR(bfgs.move[5], -0.67419986246324215, 1e-12);
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        EXPECT_DOUBLE_EQ(bfgs.pos[i], i);
        EXPECT_DOUBLE_EQ(bfgs.grad[i], i);
        for (int j = 0; j < Ions_Move_Basic::dim; ++j)
        {
            if (i == j)
                EXPECT_DOUBLE_EQ(bfgs.inv_hess(i, j), 1.0);
            else
                EXPECT_DOUBLE_EQ(bfgs.inv_hess(i, j), 0.0);
        }
    }
}

// Test the bfgs_routine() function case 2
TEST_F(IonsMoveBFGSTest, BfgsRoutineCase2)
{
    // Initilize data
    bfgs.init_done = false;
    bfgs.allocate();
    bfgs.tr_min_hit = false;
    PARAM.input.test_relax_method = 0;
    PARAM.input.out_level = "none";
    double lat0 = 1.0;
    const int istep = 1;
    std::vector<double> etot_info = {1.0, 0.9, 0.1};
    Ions_Move_Basic::relax_bfgs_rmin = -1.0;
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        bfgs.move_p[i] = i;
        bfgs.grad_p[i] = i;
        bfgs.pos_p[i] = i;
    }

    // Call the function being tested
    std::ofstream ofs("test_bfgs_routine_case2.log");
    testing::internal::CaptureStdout();
    bfgs.bfgs_routine(lat0, istep, update_iter, ofs, etot_info);
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

    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, -0.5);
    EXPECT_DOUBLE_EQ(etot_info[0], 0.9);
    EXPECT_DOUBLE_EQ(bfgs.tr_min_hit, false);
    EXPECT_NEAR(bfgs.move[0], 0.0, 1e-12);
    EXPECT_NEAR(bfgs.move[1], 0.067419986246324212, 1e-12);
    EXPECT_NEAR(bfgs.move[2], 0.13483997249264842, 1e-12);
    EXPECT_NEAR(bfgs.move[3], 0.20225995873897262, 1e-12);
    EXPECT_NEAR(bfgs.move[4], 0.26967994498529685, 1e-12);
    EXPECT_NEAR(bfgs.move[5], 0.33709993123162107, 1e-12);
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        EXPECT_DOUBLE_EQ(bfgs.pos[i], i);
        EXPECT_DOUBLE_EQ(bfgs.grad[i], i);
        for (int j = 0; j < Ions_Move_Basic::dim; ++j)
        {
            EXPECT_DOUBLE_EQ(bfgs.inv_hess(i, j), 0.0);
        }
    }
}

// Test the bfgs_routine() function case 3
TEST_F(IonsMoveBFGSTest, BfgsRoutineCase3)
{
    // Initilize data
    double lat0 = 1.0;
    const int istep = 1;
    std::vector<double> etot_info = {0.9, 1.0, -0.1};
    update_iter = 0;
    Ions_Move_Basic::largest_grad = 0.0;
    Ions_Move_Basic::relax_bfgs_init = 0.3;
    Ions_Move_Basic::best_xxx = -0.4;
    bfgs.init_done = false;
    bfgs.allocate();
    bfgs.bfgs_ndim = 1;
    bfgs.grad[0] = 1.0;
    bfgs.grad[1] = 2.0;
    bfgs.inv_hess(0, 0) = -3.0;
    bfgs.inv_hess(0, 1) = -4.0;
    bfgs.inv_hess(1, 0) = -5.0;
    bfgs.inv_hess(1, 1) = -6.0;

    // Call the function being tested
    std::ofstream ofs("test_bfgs_routine_case3.log");
    bfgs.bfgs_routine(lat0, istep, update_iter, ofs, etot_info);
    ofs.close();

    // Check the results
    std::ifstream ifs("test_bfgs_routine_case3.log");
    std::string ofs_output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("test_bfgs_routine_case3.log");

    std::string expected_ofs = " check the norm of new move 410 (Bohr)\n Uphill move : resetting bfgs history\n";

     EXPECT_THAT(ofs_output, ::testing::HasSubstr(expected_ofs));
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, 0.2);
    EXPECT_DOUBLE_EQ(etot_info[0], 0.9);
    EXPECT_DOUBLE_EQ(bfgs.tr_min_hit, false);
    EXPECT_NEAR(bfgs.move[0], -0.089442719099991588, 1e-12);
    EXPECT_NEAR(bfgs.move[1], -0.17888543819998318, 1e-12);
    EXPECT_DOUBLE_EQ(bfgs.move[2], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.move[3], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.move[4], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.move[5], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.pos[0], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.pos[1], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.pos[2], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.pos[3], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.pos[4], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.pos[5], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.grad[0], 1.0);
    EXPECT_DOUBLE_EQ(bfgs.grad[1], 2.0);
    EXPECT_DOUBLE_EQ(bfgs.grad[2], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.grad[3], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.grad[4], 0.0);
    EXPECT_DOUBLE_EQ(bfgs.grad[5], 0.0);
}

// Test the bfgs_routine() function warning quit 1
TEST_F(IonsMoveBFGSTest, BfgsRoutineWarningQuit1)
{
    // Initilize data
    bfgs.init_done = false;
    bfgs.allocate();
    bfgs.tr_min_hit = true;
    PARAM.input.test_relax_method = 1;
    PARAM.input.out_level = "ie";
    double lat0 = 1.0;
    const int istep = 1;
    std::vector<double> etot_info = {1.0, 0.9, 0.1};
    Ions_Move_Basic::relax_bfgs_rmin = 1.0;
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        bfgs.move_p[i] = 0.0;
        bfgs.grad_p[i] = i;
        bfgs.pos_p[i] = i;
    }

    // Check the results
    std::ofstream ofs("test_bfgs_routine_warning_quit1.log");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(bfgs.bfgs_routine(lat0, istep, update_iter, ofs, etot_info), ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    ofs.close();
    std::remove("test_bfgs_routine_warning_quit1.log");
    EXPECT_THAT(output, testing::HasSubstr("trust radius is too small! Break down."));
}

// Test the bfgs_routine() function warning quit 2
TEST_F(IonsMoveBFGSTest, BfgsRoutineWarningQuit2)
{
    // Initilize data
    bfgs.init_done = false;
    bfgs.allocate();
    bfgs.tr_min_hit = false;
    PARAM.input.test_relax_method = 1;
    PARAM.input.out_level = "ie";
    double lat0 = 1.0;
    const int istep = 1;
    std::vector<double> etot_info = {1.0, 0.9, 0.1};
    Ions_Move_Basic::relax_bfgs_rmin = 1.0;

    // Check the results
    std::ofstream ofs("test_bfgs_routine_warning_quit2.log");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(bfgs.bfgs_routine(lat0, istep, update_iter, ofs, etot_info), ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    ofs.close();
    std::remove("test_bfgs_routine_warning_quit2.log");
    EXPECT_THAT(output, testing::HasSubstr("BFGS: move-length unreasonably short"));
}
