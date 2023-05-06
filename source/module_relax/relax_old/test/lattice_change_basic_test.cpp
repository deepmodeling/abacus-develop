#include "module_relax/relax_old/lattice_change_basic.h"

#include "for_test.h"
#include "gtest/gtest.h"
/************************************************
 *  unit tests of namespace Lattice_Change_Basic
 ***********************************************/

/**
 * - Tested Functions:
 *   - Lattice_Change_Basic::setup_gradient()
 *   - Lattice_Change_Basic::change_lattice()
 *   - Lattice_Change_Basic::check_converged()
 *   - Lattice_Change_Basic::terminate()
 *   - Lattice_Change_Basic::setup_etot()
 */

// Define a fixture for the tests
class LatticeChangeBasicTest : public ::testing::Test
{
  protected:
    ModuleBase::matrix stress;
    UnitCell ucell;
    double lat[9], grad[9], move[9];

    virtual void SetUp()
    {
        // Initialize variables before each test
        stress.create(3, 3);
    }

    virtual void TearDown()
    {
        // Clean up after each test
    }
};

// Test the setup_gradient function with fixed_axes is volume
TEST_F(LatticeChangeBasicTest, SetupGradientVolume)
{
    // Initialize variables
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    stress(0, 0) = 1.0;
    stress(0, 1) = 2.0;
    stress(0, 2) = 3.0;
    stress(1, 0) = 4.0;
    stress(1, 1) = 5.0;
    stress(1, 2) = 6.0;
    stress(2, 0) = 7.0;
    stress(2, 1) = 8.0;
    stress(2, 2) = 9.0;

    Lattice_Change_Basic::fixed_axes = "volume";

    // Call setup_gradient method
    Lattice_Change_Basic::setup_gradient(ucell, lat, grad, stress);

    // Check expected values for stress
    EXPECT_DOUBLE_EQ(stress(0, 0), -4.0);
    EXPECT_DOUBLE_EQ(stress(1, 1), 0.0);
    EXPECT_DOUBLE_EQ(stress(2, 2), 4.0);

    // Check expected values for lat
    EXPECT_DOUBLE_EQ(lat[0], 10.0);
    EXPECT_DOUBLE_EQ(lat[1], 0.0);
    EXPECT_DOUBLE_EQ(lat[2], 0.0);
    EXPECT_DOUBLE_EQ(lat[3], 0.0);
    EXPECT_DOUBLE_EQ(lat[4], 10.0);
    EXPECT_DOUBLE_EQ(lat[5], 0.0);
    EXPECT_DOUBLE_EQ(lat[6], 0.0);
    EXPECT_DOUBLE_EQ(lat[7], 0.0);
    EXPECT_DOUBLE_EQ(lat[8], 10.0);

    // Check expected values for grad
    EXPECT_DOUBLE_EQ(grad[0], 40.0);
    EXPECT_DOUBLE_EQ(grad[1], -20.0);
    EXPECT_DOUBLE_EQ(grad[2], -30.0);
    EXPECT_DOUBLE_EQ(grad[3], -40.0);
    EXPECT_DOUBLE_EQ(grad[4], 0.0);
    EXPECT_DOUBLE_EQ(grad[5], -60.0);
    EXPECT_DOUBLE_EQ(grad[6], -70.0);
    EXPECT_DOUBLE_EQ(grad[7], -80.0);
    EXPECT_DOUBLE_EQ(grad[8], -40.0);
}

// Test the setup_gradient function with fixed_axes is not volume
TEST_F(LatticeChangeBasicTest, SetupGradientNone)
{
    // Initialize variables
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    stress(0, 0) = 1.0;
    stress(0, 1) = 2.0;
    stress(0, 2) = 3.0;
    stress(1, 0) = 4.0;
    stress(1, 1) = 5.0;
    stress(1, 2) = 6.0;
    stress(2, 0) = 7.0;
    stress(2, 1) = 8.0;
    stress(2, 2) = 9.0;

    Lattice_Change_Basic::fixed_axes = "None";

    // Call setup_gradient method
    Lattice_Change_Basic::setup_gradient(ucell, lat, grad, stress);

    // Check expected values for grad
    EXPECT_DOUBLE_EQ(grad[0], -10.0);
    EXPECT_DOUBLE_EQ(grad[1], -20.0);
    EXPECT_DOUBLE_EQ(grad[2], -30.0);
    EXPECT_DOUBLE_EQ(grad[3], -40.0);
    EXPECT_DOUBLE_EQ(grad[4], -50.0);
    EXPECT_DOUBLE_EQ(grad[5], -60.0);
    EXPECT_DOUBLE_EQ(grad[6], -70.0);
    EXPECT_DOUBLE_EQ(grad[7], -80.0);
    EXPECT_DOUBLE_EQ(grad[8], -90.0);
}

TEST_F(LatticeChangeBasicTest, ChangeLattice)
{
    // Initialize variables
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    lat[0] = 1.0;
    lat[1] = 0.0;
    lat[2] = 0.0;
    lat[3] = 0.0;
    lat[4] = 2.0;
    lat[5] = 0.0;
    lat[6] = 0.0;
    lat[7] = 0.0;
    lat[8] = 3.0;

    move[0] = 1.0;
    move[1] = 0.0;
    move[2] = 0.0;
    move[3] = 0.0;
    move[4] = 2.0;
    move[5] = 0.0;
    move[6] = 0.0;
    move[7] = 0.0;
    move[8] = 3.0;

    // Call change_lattice method
    Lattice_Change_Basic::change_lattice(ucell, move, lat);

    // Check expected values for ucell after lattice change
    EXPECT_DOUBLE_EQ(ucell.latvec.e11, 0.2);
    EXPECT_DOUBLE_EQ(ucell.latvec.e12, 0.0);
    EXPECT_DOUBLE_EQ(ucell.latvec.e13, 0.0);
    EXPECT_DOUBLE_EQ(ucell.latvec.e21, 0.0);
    EXPECT_DOUBLE_EQ(ucell.latvec.e22, 0.4);
    EXPECT_DOUBLE_EQ(ucell.latvec.e23, 0.0);
    EXPECT_DOUBLE_EQ(ucell.latvec.e31, 0.0);
    EXPECT_DOUBLE_EQ(ucell.latvec.e32, 0.0);
    EXPECT_DOUBLE_EQ(ucell.latvec.e33, 0.6);

    EXPECT_DOUBLE_EQ(ucell.a1.x, 0.2);
    EXPECT_DOUBLE_EQ(ucell.a1.y, 0.0);
    EXPECT_DOUBLE_EQ(ucell.a1.z, 0.0);
    EXPECT_DOUBLE_EQ(ucell.a2.x, 0.0);
    EXPECT_DOUBLE_EQ(ucell.a2.y, 0.4);
    EXPECT_DOUBLE_EQ(ucell.a2.z, 0.0);
    EXPECT_DOUBLE_EQ(ucell.a3.x, 0.0);
    EXPECT_DOUBLE_EQ(ucell.a3.y, 0.0);
    EXPECT_DOUBLE_EQ(ucell.a3.z, 0.6);

    EXPECT_DOUBLE_EQ(ucell.omega, 48.0);

    EXPECT_DOUBLE_EQ(ucell.GT.e11, 5.0);
    EXPECT_DOUBLE_EQ(ucell.GT.e12, 0.0);
    EXPECT_DOUBLE_EQ(ucell.GT.e13, 0.0);
    EXPECT_DOUBLE_EQ(ucell.GT.e21, 0.0);
    EXPECT_DOUBLE_EQ(ucell.GT.e22, 2.5);
    EXPECT_DOUBLE_EQ(ucell.GT.e23, 0.0);
    EXPECT_DOUBLE_EQ(ucell.GT.e31, 0.0);
    EXPECT_DOUBLE_EQ(ucell.GT.e32, 0.0);
    EXPECT_NEAR(ucell.GT.e33, 1.666666666666667, 1e-12);

    EXPECT_DOUBLE_EQ(ucell.G.e11, 5.0);
    EXPECT_DOUBLE_EQ(ucell.G.e12, 0.0);
    EXPECT_DOUBLE_EQ(ucell.G.e13, 0.0);
    EXPECT_DOUBLE_EQ(ucell.G.e21, 0.0);
    EXPECT_DOUBLE_EQ(ucell.G.e22, 2.5);
    EXPECT_DOUBLE_EQ(ucell.G.e23, 0.0);
    EXPECT_DOUBLE_EQ(ucell.G.e31, 0.0);
    EXPECT_DOUBLE_EQ(ucell.G.e32, 0.0);
    EXPECT_NEAR(ucell.G.e33, 1.666666666666667, 1e-12);

    EXPECT_DOUBLE_EQ(ucell.GGT.e11, 25.0);
    EXPECT_DOUBLE_EQ(ucell.GGT.e12, 0.0);
    EXPECT_DOUBLE_EQ(ucell.GGT.e13, 0.0);
    EXPECT_DOUBLE_EQ(ucell.GGT.e21, 0.0);
    EXPECT_DOUBLE_EQ(ucell.GGT.e22, 6.25);
    EXPECT_DOUBLE_EQ(ucell.GGT.e23, 0.0);
    EXPECT_DOUBLE_EQ(ucell.GGT.e31, 0.0);
    EXPECT_DOUBLE_EQ(ucell.GGT.e32, 0.0);
    EXPECT_NEAR(ucell.GGT.e33, 2.7777777777777786, 1e-12);

    EXPECT_DOUBLE_EQ(ucell.invGGT.e11, 0.04);
    EXPECT_DOUBLE_EQ(ucell.invGGT.e12, 0.0);
    EXPECT_DOUBLE_EQ(ucell.invGGT.e13, 0.0);
    EXPECT_DOUBLE_EQ(ucell.invGGT.e21, 0.0);
    EXPECT_DOUBLE_EQ(ucell.invGGT.e22, 0.16);
    EXPECT_DOUBLE_EQ(ucell.invGGT.e23, 0.0);
    EXPECT_DOUBLE_EQ(ucell.invGGT.e31, 0.0);
    EXPECT_DOUBLE_EQ(ucell.invGGT.e32, 0.0);
    EXPECT_DOUBLE_EQ(ucell.invGGT.e33, 0.36);
}

// Test for check_converged when ucell.lc[0] == 1 && ucell.lc[1] == 1 && ucell.lc[2] == 1, but not converged
TEST_F(LatticeChangeBasicTest, CheckConvergedCase1)
{
    // Set up test data
    Lattice_Change_Basic::update_iter = 0;
    GlobalV::STRESS_THR = 10.0;
    GlobalV::ofs_running.open("log");
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    stress(0, 0) = 1.0;
    stress(0, 1) = 2.0;
    stress(0, 2) = 3.0;
    stress(1, 0) = 4.0;
    stress(1, 1) = 5.0;
    stress(1, 2) = 6.0;
    stress(2, 0) = 7.0;
    stress(2, 1) = 8.0;
    stress(2, 2) = 9.0;

    // Call the function under test
    Lattice_Change_Basic::check_converged(ucell, stress, grad);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string expected_output = "\n Lattice relaxation is not converged yet (threshold is 10)\n";
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_EQ(output, expected_output);
    EXPECT_EQ(Lattice_Change_Basic::update_iter, 0);
    EXPECT_NEAR(Lattice_Change_Basic::largest_grad, 1323947.0517790401, 1e-12);
    EXPECT_FALSE(Lattice_Change_Basic::converged);

    ifs.close();
    std::remove("log");
}

// Test for check_converged when ucell.lc[0] == 1 && ucell.lc[1] == 1 && ucell.lc[2] == 1 && largest_grad == 0
TEST_F(LatticeChangeBasicTest, CheckConvergedCase2)
{
    // Set up test data
    Lattice_Change_Basic::update_iter = 0;
    GlobalV::STRESS_THR = 10.0;
    GlobalV::ofs_running.open("log");
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    stress(0, 0) = 0.0;
    stress(0, 1) = 0.0;
    stress(0, 2) = 0.0;
    stress(1, 0) = 0.0;
    stress(1, 1) = 0.0;
    stress(1, 2) = 0.0;
    stress(2, 0) = 0.0;
    stress(2, 1) = 0.0;
    stress(2, 2) = 0.0;

    // Call the function under test
    Lattice_Change_Basic::check_converged(ucell, stress, grad);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string expected_output = " largest stress is 0, no movement is possible.\n it may converged, otherwise no "
                                  "movement of lattice parameters is allowed.\n";
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_EQ(output, expected_output);
    EXPECT_EQ(Lattice_Change_Basic::update_iter, 0);
    EXPECT_DOUBLE_EQ(Lattice_Change_Basic::largest_grad, 0.0);
    EXPECT_TRUE(Lattice_Change_Basic::converged);

    ifs.close();
    std::remove("log");
}

// Test for check_converged when ucell.lc[0] == 1 && ucell.lc[1] == 1 && ucell.lc[2] == 1, and converged
TEST_F(LatticeChangeBasicTest, CheckConvergedCase3)
{
    // Set up test data
    Lattice_Change_Basic::update_iter = 0;
    GlobalV::STRESS_THR = 10.0;
    GlobalV::ofs_running.open("log");
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    stress(0, 0) = 0.000001;
    stress(0, 1) = 0.0;
    stress(0, 2) = 0.0;
    stress(1, 0) = 0.0;
    stress(1, 1) = 0.0;
    stress(1, 2) = 0.0;
    stress(2, 0) = 0.0;
    stress(2, 1) = 0.0;
    stress(2, 2) = 0.0;

    // Call the function under test
    Lattice_Change_Basic::check_converged(ucell, stress, grad);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string expected_output = "\n Lattice relaxation is converged!\n\n Largest gradient is = 0.147105\n";
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_EQ(output, expected_output);
    EXPECT_EQ(Lattice_Change_Basic::update_iter, 1);
    EXPECT_NEAR(Lattice_Change_Basic::largest_grad, 0.14710522797544887, 1e-12);
    EXPECT_TRUE(Lattice_Change_Basic::converged);

    ifs.close();
    std::remove("log");
}

// Test for check_converged when ucell.lc != 1, but not converged
TEST_F(LatticeChangeBasicTest, CheckConvergedCase4)
{
    // Set up test data
    Lattice_Change_Basic::update_iter = 0;
    GlobalV::STRESS_THR = 10.0;
    GlobalV::ofs_running.open("log");
    ucell.lc[0] = 0;
    ucell.lc[1] = 0;
    ucell.lc[2] = 0;
    grad[0] = 1.0;
    grad[1] = 1.0;
    grad[2] = 1.0;
    grad[3] = 1.0;
    grad[4] = 1.0;
    grad[5] = 1.0;
    grad[6] = 1.0;
    grad[7] = 1.0;
    grad[8] = 1.0;

    // Call the function under test
    Lattice_Change_Basic::check_converged(ucell, stress, grad);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string expected_output = "\n Lattice relaxation is not converged yet (threshold is 10)\n";
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_EQ(output, expected_output);
    EXPECT_EQ(Lattice_Change_Basic::update_iter, 0);
    EXPECT_NEAR(Lattice_Change_Basic::largest_grad, 147105.22797544891, 1e-12);
    EXPECT_FALSE(Lattice_Change_Basic::converged);

    ifs.close();
    std::remove("log");
}

// Test for check_converged when ucell.lc != 1, and largest_grad == 0
TEST_F(LatticeChangeBasicTest, CheckConvergedCase5)
{
    // Set up test data
    Lattice_Change_Basic::update_iter = 0;
    GlobalV::STRESS_THR = 10.0;
    GlobalV::ofs_running.open("log");
    ucell.lc[0] = 0;
    ucell.lc[1] = 0;
    ucell.lc[2] = 0;
    grad[0] = 0.0;
    grad[1] = 0.0;
    grad[2] = 0.0;
    grad[3] = 0.0;
    grad[4] = 0.0;
    grad[5] = 0.0;
    grad[6] = 0.0;
    grad[7] = 0.0;
    grad[8] = 0.0;

    // Call the function under test
    Lattice_Change_Basic::check_converged(ucell, stress, grad);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string expected_output = " largest stress is 0, no movement is possible.\n it may converged, otherwise no "
                                  "movement of lattice parameters is allowed.\n";
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_EQ(output, expected_output);
    EXPECT_EQ(Lattice_Change_Basic::update_iter, 0);
    EXPECT_DOUBLE_EQ(Lattice_Change_Basic::largest_grad, 0.0);
    EXPECT_TRUE(Lattice_Change_Basic::converged);

    ifs.close();
    std::remove("log");
}

// Test for check_converged when ucell.lc != 1, and converged
TEST_F(LatticeChangeBasicTest, CheckConvergedCase6)
{
    // Set up test data
    Lattice_Change_Basic::update_iter = 0;
    GlobalV::STRESS_THR = 10.0;
    GlobalV::ofs_running.open("log");
    ucell.lc[0] = 0;
    ucell.lc[1] = 0;
    ucell.lc[2] = 0;
    grad[0] = 0.000001;
    grad[1] = 0.0;
    grad[2] = 0.0;
    grad[3] = 0.0;
    grad[4] = 0.0;
    grad[5] = 0.0;
    grad[6] = 0.0;
    grad[7] = 0.0;
    grad[8] = 0.0;

    // Call the function under test
    Lattice_Change_Basic::check_converged(ucell, stress, grad);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string expected_output = "\n Lattice relaxation is converged!\n\n Largest gradient is = 0.147105\n";
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_EQ(output, expected_output);
    EXPECT_EQ(Lattice_Change_Basic::update_iter, 1);
    EXPECT_NEAR(Lattice_Change_Basic::largest_grad, 0.14710522797544887, 1e-12);
    EXPECT_TRUE(Lattice_Change_Basic::converged);

    ifs.close();
    std::remove("log");
}

TEST_F(LatticeChangeBasicTest, TerminateConverged)
{
    Lattice_Change_Basic::converged = true;
    Lattice_Change_Basic::stress_step = 5;
    Lattice_Change_Basic::update_iter = 10;

    std::string expected_output = " end of lattice optimization\n                              stress_step = 5\n       "
                                  "                  update iteration = 10\n";

    GlobalV::ofs_running.open("log");
    Lattice_Change_Basic::terminate();
    GlobalV::ofs_running.close();

    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    ifs.close();
    std::remove("log");
}

TEST_F(LatticeChangeBasicTest, TerminateNotConverged)
{
    Lattice_Change_Basic::converged = false;

    std::string expected_output = " the maximum number of steps has been reached.\n end of lattice optimization.\n";

    GlobalV::ofs_running.open("log");
    Lattice_Change_Basic::terminate();
    GlobalV::ofs_running.close();

    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    ifs.close();
    std::remove("log");
}

TEST_F(LatticeChangeBasicTest, SetupEtotStressStep1)
{
    Lattice_Change_Basic::stress_step = 1;
    double energy_in = 100.0;

    Lattice_Change_Basic::setup_etot(energy_in, true);

    EXPECT_DOUBLE_EQ(energy_in, Lattice_Change_Basic::etot_p);
    EXPECT_DOUBLE_EQ(energy_in, Lattice_Change_Basic::etot);
    EXPECT_DOUBLE_EQ(0.0, Lattice_Change_Basic::ediff);
}

TEST_F(LatticeChangeBasicTest, SetupEtotJudgementTrueHigherEnergy)
{
    Lattice_Change_Basic::stress_step = 2;
    double energy_in = 90.0;
    Lattice_Change_Basic::etot_p = 100.0;

    Lattice_Change_Basic::setup_etot(energy_in, true);

    EXPECT_DOUBLE_EQ(90.0, Lattice_Change_Basic::etot);
    EXPECT_DOUBLE_EQ(-10.0, Lattice_Change_Basic::ediff);
}

TEST_F(LatticeChangeBasicTest, SetupEtotJudgementTrueLowerEnergy)
{
    Lattice_Change_Basic::stress_step = 2;
    double energy_in = 100.0;
    Lattice_Change_Basic::etot_p = 90.0;

    Lattice_Change_Basic::setup_etot(energy_in, true);

    EXPECT_DOUBLE_EQ(100.0, Lattice_Change_Basic::etot);
    EXPECT_DOUBLE_EQ(0.0, Lattice_Change_Basic::ediff);
}

TEST_F(LatticeChangeBasicTest, SetupEtotJudgementFalse)
{
    Lattice_Change_Basic::stress_step = 2;
    double energy_in = 80.0;
    Lattice_Change_Basic::etot_p = 90.0;
    Lattice_Change_Basic::etot = 100.0;

    Lattice_Change_Basic::setup_etot(energy_in, false);

    EXPECT_DOUBLE_EQ(100.0, Lattice_Change_Basic::etot_p);
    EXPECT_DOUBLE_EQ(80.0, Lattice_Change_Basic::etot);
    EXPECT_DOUBLE_EQ(-20.0, Lattice_Change_Basic::ediff);
}
