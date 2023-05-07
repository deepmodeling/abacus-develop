#include "module_relax/relax_old/ions_move_basic.h"

#include "for_test.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
/************************************************
 *  unit tests of namespace Ions_Move_Basic
 ***********************************************/

/**
 * - Tested Functions:
 *   - Ions_Move_Basic::setup_gradient()
 *   - Ions_Move_Basic::move_atoms()
 *   - Ions_Move_Basic::check_converged()
 *   - Ions_Move_Basic::terminate()
 *   - Ions_Move_Basic::setup_etot()
 *   - Ions_Move_Basic::dot_func()
 *   - Ions_Move_Basic::third_order()
 */

// Define a fixture for the tests
class IonsMoveBasicTest : public ::testing::Test
{
  protected:
    UnitCell ucell;
    double move[6], pos[6], grad[6];
    ModuleBase::matrix force;

    virtual void SetUp()
    {
        // Initialize variables before each test
        force.create(ucell.nat, 3);
        force(0, 0) = 1.0;
        force(0, 1) = 2.0;
        force(0, 2) = 3.0;
        force(1, 0) = 4.0;
        force(1, 1) = 5.0;
        force(1, 2) = 6.0;
    }

    virtual void TearDown()
    {
        // Clean up after each test
    }
};

// Test the setup_gradient() function
TEST_F(IonsMoveBasicTest, SetupGradient)
{
    // Call the function being tested
    Ions_Move_Basic::dim = 6;
    Ions_Move_Basic::setup_gradient(ucell, force, pos, grad);

    // Check that the expected positions and gradients were generated
    EXPECT_DOUBLE_EQ(pos[0], 0.0);
    EXPECT_DOUBLE_EQ(pos[1], 10.0);
    EXPECT_DOUBLE_EQ(pos[2], 20.0);
    EXPECT_DOUBLE_EQ(pos[3], 30.0);
    EXPECT_DOUBLE_EQ(pos[4], 40.0);
    EXPECT_DOUBLE_EQ(pos[5], 50.0);
    EXPECT_DOUBLE_EQ(grad[0], -10.0);
    EXPECT_DOUBLE_EQ(grad[1], -20.0);
    EXPECT_DOUBLE_EQ(grad[2], -30.0);
    EXPECT_DOUBLE_EQ(grad[3], -40.0);
    EXPECT_DOUBLE_EQ(grad[4], -50.0);
    EXPECT_DOUBLE_EQ(grad[5], -60.0);
}

// Test the move_atoms() function
TEST_F(IonsMoveBasicTest, MoveAtoms)
{
    // Initialize data
    Ions_Move_Basic::dim = 6;
    GlobalV::test_relax_method = 1;
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        pos[i] = 0.0;
        move[i] = i;
    }

    // Call the function being tested
    GlobalV::ofs_running.open("log");
    Ions_Move_Basic::move_atoms(ucell, move, pos);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string expected_output = "\n movement of ions (unit is Bohr) : \n         Atom              x              y  "
                                  "            z\n       move_1              0              1              2\n       "
                                  "move_2              3              4              5\n";
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_EQ(output, expected_output);
    EXPECT_DOUBLE_EQ(pos[0], 0.0);
    EXPECT_DOUBLE_EQ(pos[1], 1.0);
    EXPECT_DOUBLE_EQ(pos[2], 2.0);
    EXPECT_DOUBLE_EQ(pos[3], 3.0);
    EXPECT_DOUBLE_EQ(pos[4], 4.0);
    EXPECT_DOUBLE_EQ(pos[5], 5.0);
}

// Test the check_converged() function case 1
TEST_F(IonsMoveBasicTest, CheckConvergedCase1)
{
    // Initialize data
    Ions_Move_Basic::dim = 6;
    Ions_Move_Basic::update_iter = 1;
    GlobalV::test_relax_method = 1;
    GlobalV::OUT_LEVEL = "ie";
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        grad[i] = 0.0;
    }

    // Call the function being tested
    GlobalV::ofs_running.open("log");
    testing::internal::CaptureStdout();
    Ions_Move_Basic::check_converged(ucell, grad);
    std::string std_outout = testing::internal::GetCapturedStdout();
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string ofs_output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    std::string expected_ofs
        = "                    old total energy (ry) = 0\n                    new total energy (ry) = 0\n              "
          "     energy difference (ry) = 0\n               largest gradient (ry/bohr) = 0\n largest force is 0, no "
          "movement is possible.\n it may converged, otherwise no movement of atom is allowed.\n";
    std::string expected_std = " ETOT DIFF (eV)       : 0\n LARGEST GRAD (eV/A)  : 0\n";

    EXPECT_EQ(expected_ofs, ofs_output);
    EXPECT_EQ(expected_std, std_outout);
    EXPECT_EQ(Ions_Move_Basic::update_iter, 1);
    EXPECT_EQ(Ions_Move_Basic::converged, true);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.0);
}

// Test the check_converged() function case 2
TEST_F(IonsMoveBasicTest, CheckConvergedCase2)
{
    // Initialize data
    Ions_Move_Basic::dim = 6;
    Ions_Move_Basic::update_iter = 1;
    Ions_Move_Basic::ediff = 0.0;
    GlobalV::test_relax_method = 1;
    GlobalV::OUT_LEVEL = "ie";
    GlobalV::FORCE_THR = 1.0;
    grad[0] = 1.0;

    // Call the function being tested
    GlobalV::ofs_running.open("log");
    testing::internal::CaptureStdout();
    Ions_Move_Basic::check_converged(ucell, grad);
    std::string std_outout = testing::internal::GetCapturedStdout();
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string ofs_output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    std::string expected_ofs
        = "                    old total energy (ry) = 0\n                    new total energy (ry) = 0\n              "
          "     energy difference (ry) = 0\n               largest gradient (ry/bohr) = 0.1\n\n Ion relaxation is "
          "converged!\n\n Energy difference (Ry) = 0\n\n Largest gradient is (eV/A) = 2.57111\n";
    std::string expected_std = " ETOT DIFF (eV)       : 0\n LARGEST GRAD (eV/A)  : 2.57111\n";

    EXPECT_EQ(expected_ofs, ofs_output);
    EXPECT_EQ(expected_std, std_outout);
    EXPECT_EQ(Ions_Move_Basic::update_iter, 2);
    EXPECT_EQ(Ions_Move_Basic::converged, true);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.1);
}

// Test the check_converged() function case 3
TEST_F(IonsMoveBasicTest, CheckConvergedCase3)
{
    // Initialize data
    Ions_Move_Basic::dim = 6;
    Ions_Move_Basic::update_iter = 1;
    Ions_Move_Basic::ediff = 1.0;
    GlobalV::test_relax_method = 1;
    GlobalV::OUT_LEVEL = "ie";
    GlobalV::FORCE_THR = 1.0;
    grad[0] = 1.0;

    // Call the function being tested
    GlobalV::ofs_running.open("log");
    testing::internal::CaptureStdout();
    Ions_Move_Basic::check_converged(ucell, grad);
    std::string std_outout = testing::internal::GetCapturedStdout();
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string ofs_output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    std::string expected_ofs
        = "                    old total energy (ry) = 0\n                    new total energy (ry) = 0\n              "
          "     energy difference (ry) = 1\n               largest gradient (ry/bohr) = 0.1\n\n Ion relaxation is not "
          "converged yet (threshold is 25.7111)\n";
    std::string expected_std = " ETOT DIFF (eV)       : 13.6057\n LARGEST GRAD (eV/A)  : 2.57111\n";

    EXPECT_EQ(expected_ofs, ofs_output);
    EXPECT_EQ(expected_std, std_outout);
    EXPECT_EQ(Ions_Move_Basic::update_iter, 1);
    EXPECT_EQ(Ions_Move_Basic::converged, false);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.1);
}

// Test the terminate() function when converged
TEST_F(IonsMoveBasicTest, TerminateConverged)
{
    // Initialize data
    Ions_Move_Basic::converged = true;
    Ions_Move_Basic::istep = 2;
    Ions_Move_Basic::update_iter = 5;

    // Call the function being tested
    GlobalV::ofs_running.open("log");
    Ions_Move_Basic::terminate(ucell);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string ofs_output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    std::string expected_ofs = " end of geometry optimization\n                                    istep = 2\n         "
                               "                update iteration = 5\n";

    EXPECT_EQ(expected_ofs, ofs_output);
}

// Test the terminate() function when not converged
TEST_F(IonsMoveBasicTest, TerminateNotConverged)
{
    // Initialize data
    Ions_Move_Basic::converged = false;

    // Call the function being tested
    GlobalV::ofs_running.open("log");
    Ions_Move_Basic::terminate(ucell);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string ofs_output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    std::string expected_ofs = " the maximum number of steps has been reached.\n end of geometry optimization.\n";

    EXPECT_EQ(expected_ofs, ofs_output);
}

// Test the setup_etot() function case 1
TEST_F(IonsMoveBasicTest, SetupEtotCase1)
{
    // Initialize data
    Ions_Move_Basic::istep = 1;
    Ions_Move_Basic::etot_p = 1.0;
    Ions_Move_Basic::etot = 2.0;
    Ions_Move_Basic::ediff = 0.0;
    double energy_in = 3.0;
    bool judgement = true;

    // Call the function being tested
    Ions_Move_Basic::setup_etot(energy_in, judgement);

    // Check the results
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::etot_p, 3.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::etot, 3.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::ediff, 0.0);
}

// Test the setup_etot() function case 2
TEST_F(IonsMoveBasicTest, SetupEtotCase2)
{
    // Initialize data
    Ions_Move_Basic::istep = 2;
    Ions_Move_Basic::etot_p = 4.0;
    Ions_Move_Basic::etot = 2.0;
    Ions_Move_Basic::ediff = 0.0;
    double energy_in = 3.0;
    bool judgement = true;

    // Call the function being tested
    Ions_Move_Basic::setup_etot(energy_in, judgement);

    // Check the results
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::etot_p, 3.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::etot, 3.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::ediff, -1.0);
}

// Test the setup_etot() function case 3
TEST_F(IonsMoveBasicTest, SetupEtotCase3)
{
    // Initialize data
    Ions_Move_Basic::istep = 2;
    Ions_Move_Basic::etot_p = 1.0;
    Ions_Move_Basic::etot = 2.0;
    Ions_Move_Basic::ediff = 0.0;
    double energy_in = 3.0;
    bool judgement = true;

    // Call the function being tested
    Ions_Move_Basic::setup_etot(energy_in, judgement);

    // Check the results
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::etot_p, 1.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::etot, 3.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::ediff, 0.0);
}

// Test the setup_etot() function case 4
TEST_F(IonsMoveBasicTest, SetupEtotCase4)
{
    // Initialize data
    Ions_Move_Basic::istep = 2;
    Ions_Move_Basic::etot_p = 1.0;
    Ions_Move_Basic::etot = 2.0;
    Ions_Move_Basic::ediff = 0.0;
    double energy_in = 3.0;
    bool judgement = false;

    // Call the function being tested
    Ions_Move_Basic::setup_etot(energy_in, judgement);

    // Check the results
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::etot_p, 2.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::etot, 3.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::ediff, 1.0);
}

// Test the dot_func() function
TEST_F(IonsMoveBasicTest, DotFunc)
{
    // Initialize data
    double dim_in = 3;
    double a[3] = {1.0, 2.0, 3.0};
    double b[3] = {1.0, 2.0, 3.0};

    // Call the function being tested
    double result = Ions_Move_Basic::dot_func(a, b, dim_in);

    // Check the results
    EXPECT_DOUBLE_EQ(result, 14.0);
}