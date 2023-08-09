#include "for_test.h"
#define private public
#define protected public
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_relax/relax_old/ions_move_basic.h"
#include "module_relax/relax_old/ions_move_bfgs.h"
/************************************************
 *  unit tests of class Ions_Move_BFGS
 ***********************************************/

/**
 * - Tested Functions:
 *   - Ions_Move_BFGS::allocate()
 *   - Ions_Move_BFGS::start()
 *   - Ions_Move_BFGS::bfgs_routine()
 *   - Ions_Move_BFGS::restart_bfgs()
 */

// Define a fixture for the tests
class IonsMoveBFGSTest : public ::testing::Test
{
  protected:
    Ions_Move_BFGS bfgs;

    virtual void SetUp()
    {
        // Initialize variables before each test
        Ions_Move_Basic::dim = 6;
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
    bfgs.init_done = false;
    bfgs.save_flag = true;

    // Call the function being tested
    bfgs.allocate();
    GlobalV::ofs_running.open("log");
    bfgs.start(ucell, force, energy_in);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_THAT(output, testing::HasSubstr("update iteration"));
}

// Test the start() function case 2
TEST_F(IonsMoveBFGSTest, StartCase2)
{
    // Initilize data
    UnitCell ucell;
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 10.0;
    double energy_in = 0.0;
    bfgs.init_done = false;
    bfgs.save_flag = true;

    // Call the function being tested
    bfgs.allocate();
    GlobalV::ofs_running.open("log");
    bfgs.start(ucell, force, energy_in);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_THAT(output, testing::HasSubstr("Ion relaxation is not converged yet"));
}

// Test the restart_bfgs() function case 1
TEST_F(IonsMoveBFGSTest, RestartBfgsCase1)
{
    // Initilize data
    bfgs.init_done = false;
    GlobalV::test_relax_method = 1;
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
    GlobalV::ofs_running.open("log");
    bfgs.restart_bfgs(lat0);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    std::string expected_output = "                  trust_radius_old (bohr) = 2.44949\n";

    EXPECT_EQ(output, expected_output);
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
    GlobalV::test_relax_method = 1;
    double lat0 = 1.0;
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        bfgs.move_p[i] = 1.0;
        bfgs.pos[i] = i;
        bfgs.pos_p[i] = i;
    }

    // Call the function being tested
    bfgs.restart_bfgs(lat0);

    // Check the results
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::update_iter, 0.0);
    EXPECT_DOUBLE_EQ(bfgs.tr_min_hit, false);
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        EXPECT_DOUBLE_EQ(bfgs.pos_p[i], 0.0);
        EXPECT_DOUBLE_EQ(bfgs.grad_p[i], 0.0);
        EXPECT_DOUBLE_EQ(bfgs.move_p[i], 0.0);
        for (int j = 0; j < Ions_Move_Basic::dim; ++j)
        {
            if (i == j)
                EXPECT_DOUBLE_EQ(bfgs.inv_hess(i, j), 1.0);
            else
                EXPECT_DOUBLE_EQ(bfgs.inv_hess(i, j), 0.0);
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
    GlobalV::test_relax_method = 1;
    GlobalV::OUT_LEVEL = "ie";
    double lat0 = 1.0;
    Ions_Move_Basic::etot = 1.0;
    Ions_Move_Basic::etot_p = 0.9;
    Ions_Move_Basic::relax_bfgs_rmin = 1.0;
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        bfgs.move_p[i] = 0.0;
        bfgs.grad_p[i] = i;
        bfgs.pos_p[i] = i;
    }

    // Call the function being tested
    GlobalV::ofs_running.open("log");
    testing::internal::CaptureStdout();
    bfgs.bfgs_routine(lat0);
    std::string std_outout = testing::internal::GetCapturedStdout();
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string ofs_output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    std::string expected_ofs
        = "                                     dE0s = 0\n                                      den = 0.1\n            "
          "    interpolated trust radius = 0\ntrust_radius = 0\nrelax_bfgs_rmin = 1\nrelax_bfgs_rmax = -1\n "
          "trust_radius < relax_bfgs_rmin, reset bfgs history.\n                                    istep = 0\n        "
          "                 update iteration = 0\n";
    std::string expected_std = " BFGS TRUST (Bohr)    : 1\n";

    EXPECT_EQ(expected_ofs, ofs_output);
    EXPECT_EQ(expected_std, std_outout);

    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, 1.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::etot, 0.9);
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
    GlobalV::test_relax_method = 0;
    GlobalV::OUT_LEVEL = "none";
    double lat0 = 1.0;
    Ions_Move_Basic::etot = 1.0;
    Ions_Move_Basic::etot_p = 0.9;
    Ions_Move_Basic::relax_bfgs_rmin = -1.0;
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        bfgs.move_p[i] = i;
        bfgs.grad_p[i] = i;
        bfgs.pos_p[i] = i;
    }

    // Call the function being tested
    GlobalV::ofs_running.open("log");
    testing::internal::CaptureStdout();
    bfgs.bfgs_routine(lat0);
    std::string std_outout = testing::internal::GetCapturedStdout();
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string ofs_output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    std::string expected_ofs = " quadratic interpolation is impossible.\n                                    istep = "
                               "0\n                         update iteration = 0\n";
    std::string expected_std = "";

    EXPECT_EQ(expected_ofs, ofs_output);
    EXPECT_EQ(expected_std, std_outout);

    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, -0.5);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::etot, 0.9);
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
    Ions_Move_Basic::etot = 0.9;
    Ions_Move_Basic::etot_p = 1.0;
    Ions_Move_Basic::update_iter = 0;
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
    GlobalV::ofs_running.open("log");
    bfgs.bfgs_routine(lat0);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string ofs_output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    std::string expected_ofs = " check the norm of new move 410 (Bohr)\n Uphill move : resetting bfgs history\n        "
                               "                            istep = 0\n                         update iteration = 1\n";

    EXPECT_EQ(expected_ofs, ofs_output);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, 0.2);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::etot, 0.9);
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
    GlobalV::test_relax_method = 1;
    GlobalV::OUT_LEVEL = "ie";
    double lat0 = 1.0;
    Ions_Move_Basic::etot = 1.0;
    Ions_Move_Basic::etot_p = 0.9;
    Ions_Move_Basic::relax_bfgs_rmin = 1.0;
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        bfgs.move_p[i] = 0.0;
        bfgs.grad_p[i] = i;
        bfgs.pos_p[i] = i;
    }

    // Check the results
    testing::internal::CaptureStdout();
    EXPECT_EXIT(bfgs.bfgs_routine(lat0), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("trust radius is too small! Break down."));
}

// Test the bfgs_routine() function warning quit 2
TEST_F(IonsMoveBFGSTest, BfgsRoutineWarningQuit2)
{
    // Initilize data
    bfgs.init_done = false;
    bfgs.allocate();
    bfgs.tr_min_hit = false;
    GlobalV::test_relax_method = 1;
    GlobalV::OUT_LEVEL = "ie";
    double lat0 = 1.0;
    Ions_Move_Basic::etot = 1.0;
    Ions_Move_Basic::etot_p = 0.9;
    Ions_Move_Basic::relax_bfgs_rmin = 1.0;

    // Check the results
    testing::internal::CaptureStdout();
    EXPECT_EXIT(bfgs.bfgs_routine(lat0), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("BFGS: move-length unreasonably short"));
}

#include "mpi.h"
int main(int argc, char **argv)
{
        MPI_Init(&argc, &argv);
        testing::InitGoogleTest(&argc, argv);

        int nproc;
        int my_rank;

        MPI_Comm_size(MPI_COMM_WORLD,&nproc);
        MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

        int result = RUN_ALL_TESTS();

        MPI_Finalize();

        return result;
}