#include <regex>
#include "for_test.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "source_io/module_parameter/parameter.h"
#include "source_relax/ions_move_basic.h"
#include "source_relax/ions_move_sd.h"
#undef private

/************************************************
 *  unit tests of class Ions_Move_SD
 ***********************************************/

class IonsMoveSDTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        // Initialize variables before each test
        Ions_Move_Basic::dim = 6;
        update_iter = 5;
        im_sd.allocate();
        PARAM.input.force_thr = 0.001;
    }

    void TearDown() override
    {
        // Clean up after each test
    }

    Ions_Move_SD im_sd;
    int update_iter;
};

// Test whether the allocate() function can correctly allocate memory space
TEST_F(IonsMoveSDTest, TestAllocate)
{
    Ions_Move_Basic::dim = 4;
    im_sd.allocate();

    // Check if allocated vectors are not empty
    EXPECT_EQ(im_sd.grad_saved.size(), 4U);
    EXPECT_EQ(im_sd.pos_saved.size(), 4U);
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
    const int istep = 1;
    UnitCell ucell;
    ModuleBase::matrix force(2, 3);
    double etot = 0.0;
    std::vector<double> etot_info(2, 0.0);

    // call function
    std::ofstream ofs("test_sd_start_converged.log");
    im_sd.start(ucell, force, etot, istep, update_iter, ofs, etot_info);
    ofs.close();

    // Check output
    std::string expected_output = "\n Largest force is 0 eV/Angstrom while threshold is -1 eV/Angstrom\n"
                                  " largest force is 0, no movement is possible.\n it may converged, otherwise no "
                                  "movement of atom is allowed.\n end of geometry optimization\n                       "
                                  "             istep = 1\n                         update iteration = 5\n";
    std::ifstream ifs("test_sd_start_converged.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("test_sd_start_converged.log");

    std::regex pattern(R"(==> .*::.*\t[\d\.]+ GB\t\d+ s\n )");
    output = std::regex_replace(output, pattern, "");
    EXPECT_THAT(output, testing::HasSubstr(expected_output));
    EXPECT_EQ(update_iter, 5);
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
    const int istep = 1;
    UnitCell ucell;
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 1.0;
    double etot = 0.0;
    std::vector<double> etot_info(2, 0.0);
    for (int it = 0; it < ucell.ntype; it++)
    {
        Atom* atom = &ucell.atoms[it];
        for (int ia = 0; ia < atom->na; ia++)
        {
            for (int ik = 0; ik < 3; ++ik)
            {
                atom->tau[ia][ik] = (ik + 1)/3;
                atom->mbl[ia][ik] = 1;
            }
        }
    }

    // call function
    std::ofstream ofs("test_sd_start_not_converged.log");
    im_sd.start(ucell, force, etot, istep, update_iter, ofs, etot_info);
    ofs.close();

    // Check output
    std::string expected_output = "\n Largest force is 25.7111 eV/Angstrom while threshold is -1 eV/Angstrom\n\n"
                                  " Ion relaxation is not converged yet (threshold is 0.0257111)\n";
    std::ifstream ifs("test_sd_start_not_converged.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("test_sd_start_not_converged.log");

    EXPECT_THAT(output, testing::HasSubstr(expected_output));
    EXPECT_EQ(update_iter, 6);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 1.0);
    EXPECT_DOUBLE_EQ(im_sd.energy_saved, 0.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[0], -1.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[1], 0.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[2], 10.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[3], 0.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[4], 0.0);
    EXPECT_DOUBLE_EQ(im_sd.pos_saved[5], 10.0);
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
    const int istep = 1;
    PARAM.input.out_level = "ie";
    std::vector<double> etot_info(2, 0.0);

    // call function
    testing::internal::CaptureStdout();
    im_sd.cal_tradius_sd(istep, etot_info);
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
    const int istep = 2;
    std::vector<double> etot_info = {0.0, 1.0};
    PARAM.input.out_level = "m";

    // call function
    im_sd.cal_tradius_sd(istep, etot_info);

    // Check the results
    EXPECT_EQ(Ions_Move_Basic::trust_radius, -1.0);
}

// Test function cal_tradius_sd() case 3
TEST_F(IonsMoveSDTest, CalTradiusSdCase3)
{
    // setup data
    const int istep = 2;
    std::vector<double> etot_info = {1.0, 0.0};
    PARAM.input.out_level = "m";

    // call function
    im_sd.cal_tradius_sd(istep, etot_info);

    // Check the results
    EXPECT_EQ(Ions_Move_Basic::trust_radius, -0.5);
}

// Test function cal_tradius_sd() warning quit
TEST_F(IonsMoveSDTest, CalTradiusWraningQuit)
{
    // setup data
    const int istep = 0;
    std::vector<double> etot_info(2, 0.0);

    // Check the results
    testing::internal::CaptureStdout();
    EXPECT_EXIT(im_sd.cal_tradius_sd(istep, etot_info), ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("istep < 1!"));
}
