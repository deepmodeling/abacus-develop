#include "for_test.h"
#define private public
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_relax/relax_old/ions_move_basic.h"
#include "module_relax/relax_old/ions_move_sd.h"
/************************************************
 *  unit tests of class Ions_Move_SD
 ***********************************************/

/**
 * - Tested Functions:
 *   - Ions_Move_SD::allocate()
 *   - Ions_Move_SD::start()
 *   - Ions_Move_SD::cal_tradius_sd()
 */

class IonsMoveSDTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        // Initialize variables before each test
        Ions_Move_Basic::dim = 6;
        Ions_Move_Basic::update_iter = 5;
        im_sd.allocate();
    }

    void TearDown() override
    {
        // Clean up after each test
    }

    Ions_Move_SD im_sd;
};

// Test whether the allocate() function can correctly allocate memory space
TEST_F(IonsMoveSDTest, TestAllocate)
{
    Ions_Move_Basic::dim = 4;
    im_sd.allocate();

    // Check if allocated arrays are not empty
    EXPECT_NE(nullptr, im_sd.grad_saved);
    EXPECT_NE(nullptr, im_sd.pos_saved);
}

// Test if a dimension less than or equal to 0 results in an assertion error
TEST_F(IonsMoveSDTest, TestAllocateWithZeroDimension)
{
    Ions_Move_Basic::dim = 0;
    ASSERT_DEATH(im_sd.allocate(), "");
}

// Check that the arrays are correctly initialized to 0
TEST_F(IonsMoveSDTest, TestAllocateAndInitialize)
{
    Ions_Move_Basic::dim = 4;
    im_sd.allocate();

    // Check that the arrays are correctly initialized to 0
    EXPECT_DOUBLE_EQ(im_sd.grad_saved[0], 0.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[1], 0.0);
    EXPECT_DOUBLE_EQ(im_sd.grad_saved[2], 0.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[3], 0.0);
}

// Test function start() when converged
TEST_F(IonsMoveSDTest, TestStartConverged)
{
    // setup data
    Ions_Move_Basic::istep = 1;
    Ions_Move_Basic::converged = true;
    UnitCell ucell;
    ModuleBase::matrix force(2, 3);
    double etot = 0.0;

    // call function
    GlobalV::ofs_running.open("log");
    im_sd.start(ucell, force, etot);
    GlobalV::ofs_running.close();

    // Check output
    std::string expected_output = " largest force is 0, no movement is possible.\n it may converged, otherwise no "
                                  "movement of atom is allowed.\n end of geometry optimization\n                       "
                                  "             istep = 1\n                         update iteration = 5\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_EQ(expected_output, output);
    EXPECT_EQ(Ions_Move_Basic::converged, true);
    EXPECT_EQ(Ions_Move_Basic::update_iter, 5);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.0);
    EXPECT_DOUBLE_EQ(im_sd.energy_saved, 0.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[0], 0.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[1], 10.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[2], 20.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[3], 30.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[4], 40.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[5], 50.0);
}

// Test function start() when nor converged
TEST_F(IonsMoveSDTest, TestStartNotConverged)
{
    // setup data
    Ions_Move_Basic::istep = 1;
    Ions_Move_Basic::converged = true;
    UnitCell ucell;
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 1.0;
    double etot = 0.0;

    // call function
    GlobalV::ofs_running.open("log");
    im_sd.start(ucell, force, etot);
    GlobalV::ofs_running.close();

    // Check output
    std::string expected_output = "\n Ion relaxation is not converged yet (threshold is 0.0257111)\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_EQ(expected_output, output);
    EXPECT_EQ(Ions_Move_Basic::converged, false);
    EXPECT_EQ(Ions_Move_Basic::update_iter, 6);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 1.0);
    EXPECT_DOUBLE_EQ(im_sd.energy_saved, 0.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[0], -1.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[1], 10.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[2], 20.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[3], 30.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[4], 40.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[5], 50.0);
    EXPECT_DOUBLE_EQ(im_sd.grad_saved[0], -1.0);
    EXPECT_DOUBLE_EQ(im_sd.grad_saved[1], 0.0);
    EXPECT_DOUBLE_EQ(im_sd.grad_saved[2], 0.0);
    EXPECT_DOUBLE_EQ(im_sd.grad_saved[3], 0.0);
    EXPECT_DOUBLE_EQ(im_sd.grad_saved[4], 0.0);
    EXPECT_DOUBLE_EQ(im_sd.grad_saved[5], 0.0);
}

// Test function cal_tradius_sd() case 1
TEST_F(IonsMoveSDTest, CalTradiusSdCase1)
{
    // setup data
    Ions_Move_Basic::istep = 1;
    GlobalV::OUT_LEVEL = "ie";

    // call function
    testing::internal::CaptureStdout();
    im_sd.cal_tradius_sd();
    std::string std_outout = testing::internal::GetCapturedStdout();

    // Check the results
    std::string expected_std = " SD RADIUS (Bohr)     : -1\n";
    EXPECT_EQ(expected_std, std_outout);
    EXPECT_EQ(Ions_Move_Basic::trust_radius, -1.0);
}

// Test function cal_tradius_sd() case 2
TEST_F(IonsMoveSDTest, CalTradiusSdCase2)
{
    // setup data
    Ions_Move_Basic::istep = 2;
    Ions_Move_Basic::ediff = -1.0;
    GlobalV::OUT_LEVEL = "m";

    // call function
    im_sd.cal_tradius_sd();

    // Check the results
    EXPECT_EQ(Ions_Move_Basic::trust_radius, -1.0);
}

// Test function cal_tradius_sd() case 3
TEST_F(IonsMoveSDTest, CalTradiusSdCase3)
{
    // setup data
    Ions_Move_Basic::istep = 2;
    Ions_Move_Basic::ediff = 1.0;
    GlobalV::OUT_LEVEL = "m";

    // call function
    im_sd.cal_tradius_sd();

    // Check the results
    EXPECT_EQ(Ions_Move_Basic::trust_radius, -0.5);
}

// Test function cal_tradius_sd() warning quit
TEST_F(IonsMoveSDTest, CalTradiusWraningQuit)
{
    // setup data
    Ions_Move_Basic::istep = 0;

    // Check the results
    testing::internal::CaptureStdout();
    EXPECT_EXIT(im_sd.cal_tradius_sd(), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("istep < 1!"));
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