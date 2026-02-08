#include "source_relax/lattice_change_basic.h"
#include "mock_remake_cell.h"

#include "for_test.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#define private public
#include "source_io/module_parameter/parameter.h"
#undef private

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
        // Reset mock state before each test
        unitcell::reset_remake_cell_mock();
        // Reset fixed_ibrav to default
        PARAM.input.fixed_ibrav = false;
    }

    virtual void TearDown()
    {
        // Clean up after each test
        unitcell::reset_remake_cell_mock();
        PARAM.input.fixed_ibrav = false;
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
    PARAM.input.stress_thr = 10.0;
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
    std::string expected_output = "\n Geometry relaxation is not converged because threshold is 10 kbar\n";
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
    PARAM.input.stress_thr = 10.0;
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
    std::string expected_output = " Largest stress is 0, movement is impossible.\n";
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
    PARAM.input.stress_thr = 10.0;
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
    std::string expected_output = "\n Geometry relaxation is converged!\n\n Largest stress is 0.147105 kbar while threshold is 10 kbar\n";
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
    PARAM.input.stress_thr = 10.0;
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
    std::string expected_output = "\n Geometry relaxation is not converged because threshold is 10 kbar\n";
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
    PARAM.input.stress_thr = 10.0;
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
    std::string expected_output = " Largest stress is 0, movement is impossible.\n";
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
    PARAM.input.stress_thr = 10.0;
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
    std::string expected_output = "\n Geometry relaxation is converged!\n\n Largest stress is 0.147105 kbar while threshold is 10 kbar\n";
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

// ============================================================================
// NEW TESTS FOR SHAPE CONSTRAINT, VOLUME RESCALING, AND FIXED_IBRAV
// ============================================================================

// Test the setup_gradient function with fixed_axes = "shape"
TEST_F(LatticeChangeBasicTest, SetupGradientShape)
{
    // Initialize variables
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;

    // Non-isotropic stress tensor
    stress(0, 0) = 1.0;
    stress(0, 1) = 2.0;
    stress(0, 2) = 3.0;
    stress(1, 0) = 4.0;
    stress(1, 1) = 5.0;
    stress(1, 2) = 6.0;
    stress(2, 0) = 7.0;
    stress(2, 1) = 8.0;
    stress(2, 2) = 9.0;

    Lattice_Change_Basic::fixed_axes = "shape";

    // Call setup_gradient method
    Lattice_Change_Basic::setup_gradient(ucell, lat, grad, stress);

    // Check that stress becomes isotropic (only diagonal, all equal)
    // Average pressure = (1 + 5 + 9) / 3 = 5
    EXPECT_DOUBLE_EQ(stress(0, 0), 5.0);
    EXPECT_DOUBLE_EQ(stress(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(stress(2, 2), 5.0);

    // Off-diagonal elements should be zero
    EXPECT_DOUBLE_EQ(stress(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(stress(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(stress(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(stress(1, 2), 0.0);
    EXPECT_DOUBLE_EQ(stress(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(stress(2, 1), 0.0);
}

// Test volume constraint rescaling in change_lattice
TEST_F(LatticeChangeBasicTest, ChangeLatticeVolumeRescaling)
{
    // Initialize variables
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    ucell.lat0 = 10.0;

    // Set initial lattice (cubic, volume = 1000)
    ucell.latvec.e11 = 1.0;
    ucell.latvec.e12 = 0.0;
    ucell.latvec.e13 = 0.0;
    ucell.latvec.e21 = 0.0;
    ucell.latvec.e22 = 1.0;
    ucell.latvec.e23 = 0.0;
    ucell.latvec.e31 = 0.0;
    ucell.latvec.e32 = 0.0;
    ucell.latvec.e33 = 1.0;

    ucell.omega = 1000.0; // Initial volume

    lat[0] = 10.0;
    lat[1] = 0.0;
    lat[2] = 0.0;
    lat[3] = 0.0;
    lat[4] = 10.0;
    lat[5] = 0.0;
    lat[6] = 0.0;
    lat[7] = 0.0;
    lat[8] = 10.0;

    // Apply a move that would change volume (expand by 10%)
    move[0] = 1.0;
    move[1] = 0.0;
    move[2] = 0.0;
    move[3] = 0.0;
    move[4] = 1.0;
    move[5] = 0.0;
    move[6] = 0.0;
    move[7] = 0.0;
    move[8] = 1.0;

    Lattice_Change_Basic::fixed_axes = "volume";

    // Call change_lattice method
    Lattice_Change_Basic::change_lattice(ucell, move, lat);

    // Check that volume is preserved (should still be 1000)
    EXPECT_NEAR(ucell.omega, 1000.0, 1e-8);

    // Check that lattice vectors were rescaled uniformly
    double expected_scale = std::pow(1000.0 / 1331.0, 1.0/3.0); // (old_vol / new_vol)^(1/3)
    EXPECT_NEAR(ucell.latvec.e11, 1.1 * expected_scale, 1e-10);
    EXPECT_NEAR(ucell.latvec.e22, 1.1 * expected_scale, 1e-10);
    EXPECT_NEAR(ucell.latvec.e33, 1.1 * expected_scale, 1e-10);
}

// Test volume constraint with non-cubic cell
TEST_F(LatticeChangeBasicTest, ChangeLatticeVolumeRescalingNonCubic)
{
    // Initialize variables
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    ucell.lat0 = 10.0;

    // Set initial lattice (non-cubic, volume = 1200)
    ucell.latvec.e11 = 1.0;
    ucell.latvec.e12 = 0.0;
    ucell.latvec.e13 = 0.0;
    ucell.latvec.e21 = 0.0;
    ucell.latvec.e22 = 1.2;
    ucell.latvec.e23 = 0.0;
    ucell.latvec.e31 = 0.0;
    ucell.latvec.e32 = 0.0;
    ucell.latvec.e33 = 1.0;

    ucell.omega = 1200.0; // Initial volume

    lat[0] = 10.0;
    lat[1] = 0.0;
    lat[2] = 0.0;
    lat[3] = 0.0;
    lat[4] = 12.0;
    lat[5] = 0.0;
    lat[6] = 0.0;
    lat[7] = 0.0;
    lat[8] = 10.0;

    // Apply a move that would change volume
    move[0] = 0.5;
    move[1] = 0.0;
    move[2] = 0.0;
    move[3] = 0.0;
    move[4] = -0.5;
    move[5] = 0.0;
    move[6] = 0.0;
    move[7] = 0.0;
    move[8] = 0.3;

    Lattice_Change_Basic::fixed_axes = "volume";

    // Call change_lattice method
    Lattice_Change_Basic::change_lattice(ucell, move, lat);

    // Check that volume is preserved
    EXPECT_NEAR(ucell.omega, 1200.0, 1e-8);
}

// Test change_lattice without volume constraint (should change volume)
TEST_F(LatticeChangeBasicTest, ChangeLatticeNoVolumeConstraint)
{
    // Initialize variables
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    ucell.lat0 = 10.0;

    // Set initial lattice
    ucell.latvec.e11 = 1.0;
    ucell.latvec.e12 = 0.0;
    ucell.latvec.e13 = 0.0;
    ucell.latvec.e21 = 0.0;
    ucell.latvec.e22 = 1.0;
    ucell.latvec.e23 = 0.0;
    ucell.latvec.e31 = 0.0;
    ucell.latvec.e32 = 0.0;
    ucell.latvec.e33 = 1.0;

    ucell.omega = 1000.0; // Initial volume

    lat[0] = 10.0;
    lat[1] = 0.0;
    lat[2] = 0.0;
    lat[3] = 0.0;
    lat[4] = 10.0;
    lat[5] = 0.0;
    lat[6] = 0.0;
    lat[7] = 0.0;
    lat[8] = 10.0;

    // Apply a move that changes volume
    move[0] = 1.0;
    move[1] = 0.0;
    move[2] = 0.0;
    move[3] = 0.0;
    move[4] = 1.0;
    move[5] = 0.0;
    move[6] = 0.0;
    move[7] = 0.0;
    move[8] = 1.0;

    Lattice_Change_Basic::fixed_axes = "None";

    // Call change_lattice method
    Lattice_Change_Basic::change_lattice(ucell, move, lat);

    // Check that volume DID change (should be 1331)
    EXPECT_NEAR(ucell.omega, 1331.0, 1e-8);

    // Check lattice vectors
    EXPECT_DOUBLE_EQ(ucell.latvec.e11, 1.1);
    EXPECT_DOUBLE_EQ(ucell.latvec.e22, 1.1);
    EXPECT_DOUBLE_EQ(ucell.latvec.e33, 1.1);
}

// Test fixed_ibrav with simple cubic lattice
TEST_F(LatticeChangeBasicTest, ChangeLatticeFixedIbravSimpleCubic)
{
    // Initialize variables
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    ucell.lat0 = 10.0;
    ucell.latName = "sc";

    // Set initial lattice (slightly distorted cubic)
    ucell.latvec.e11 = 1.0;
    ucell.latvec.e12 = 0.01; // Small distortion
    ucell.latvec.e13 = 0.0;
    ucell.latvec.e21 = 0.0;
    ucell.latvec.e22 = 1.0;
    ucell.latvec.e23 = 0.0;
    ucell.latvec.e31 = 0.0;
    ucell.latvec.e32 = 0.0;
    ucell.latvec.e33 = 1.0;

    lat[0] = 10.0;
    lat[1] = 0.1;
    lat[2] = 0.0;
    lat[3] = 0.0;
    lat[4] = 10.0;
    lat[5] = 0.0;
    lat[6] = 0.0;
    lat[7] = 0.0;
    lat[8] = 10.0;

    // Apply a move
    move[0] = 0.1;
    move[1] = 0.0;
    move[2] = 0.0;
    move[3] = 0.0;
    move[4] = 0.1;
    move[5] = 0.0;
    move[6] = 0.0;
    move[7] = 0.0;
    move[8] = 0.1;

    PARAM.input.fixed_ibrav = true;
    Lattice_Change_Basic::fixed_axes = "None";

    // Verify remake_cell was not called yet
    EXPECT_FALSE(unitcell::was_remake_cell_called());

    // Call change_lattice method
    Lattice_Change_Basic::change_lattice(ucell, move, lat);

    // Verify remake_cell was called
    EXPECT_TRUE(unitcell::was_remake_cell_called());

    // Check that lattice is now perfect cubic (all diagonal, equal)
    // This is enforced by the mock remake_cell function
    EXPECT_NEAR(ucell.latvec.e11, ucell.latvec.e22, 1e-10);
    EXPECT_NEAR(ucell.latvec.e22, ucell.latvec.e33, 1e-10);
    EXPECT_NEAR(ucell.latvec.e12, 0.0, 1e-10);
    EXPECT_NEAR(ucell.latvec.e13, 0.0, 1e-10);
    EXPECT_NEAR(ucell.latvec.e21, 0.0, 1e-10);
    EXPECT_NEAR(ucell.latvec.e23, 0.0, 1e-10);
    EXPECT_NEAR(ucell.latvec.e31, 0.0, 1e-10);
    EXPECT_NEAR(ucell.latvec.e32, 0.0, 1e-10);

    // Reset for other tests
    PARAM.input.fixed_ibrav = false;
}

// Test fixed_ibrav with FCC lattice
TEST_F(LatticeChangeBasicTest, ChangeLatticeFixedIbravFCC)
{
    // Initialize variables
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    ucell.lat0 = 10.0;
    ucell.latName = "fcc";

    // Set initial lattice (slightly distorted FCC)
    double celldm = 1.0;
    ucell.latvec.e11 = -celldm + 0.01; // Small distortion
    ucell.latvec.e12 = 0.0;
    ucell.latvec.e13 = celldm;
    ucell.latvec.e21 = 0.0;
    ucell.latvec.e22 = celldm;
    ucell.latvec.e23 = celldm;
    ucell.latvec.e31 = -celldm;
    ucell.latvec.e32 = celldm;
    ucell.latvec.e33 = 0.0;

    lat[0] = (-celldm + 0.01) * ucell.lat0;
    lat[1] = 0.0;
    lat[2] = celldm * ucell.lat0;
    lat[3] = 0.0;
    lat[4] = celldm * ucell.lat0;
    lat[5] = celldm * ucell.lat0;
    lat[6] = -celldm * ucell.lat0;
    lat[7] = celldm * ucell.lat0;
    lat[8] = 0.0;

    // Apply a small move
    for (int i = 0; i < 9; i++) move[i] = 0.01 * ucell.lat0;

    PARAM.input.fixed_ibrav = true;
    Lattice_Change_Basic::fixed_axes = "None";

    // Verify remake_cell was not called yet
    EXPECT_FALSE(unitcell::was_remake_cell_called());

    // Call change_lattice method
    Lattice_Change_Basic::change_lattice(ucell, move, lat);

    // Verify remake_cell was called
    EXPECT_TRUE(unitcell::was_remake_cell_called());

    // Check that lattice maintains FCC structure
    // For FCC: a1 = (-a, 0, a), a2 = (0, a, a), a3 = (-a, a, 0)
    // All should have same magnitude
    double mag1 = std::sqrt(ucell.latvec.e11*ucell.latvec.e11 +
                           ucell.latvec.e12*ucell.latvec.e12 +
                           ucell.latvec.e13*ucell.latvec.e13);
    double mag2 = std::sqrt(ucell.latvec.e21*ucell.latvec.e21 +
                           ucell.latvec.e22*ucell.latvec.e22 +
                           ucell.latvec.e23*ucell.latvec.e23);
    double mag3 = std::sqrt(ucell.latvec.e31*ucell.latvec.e31 +
                           ucell.latvec.e32*ucell.latvec.e32 +
                           ucell.latvec.e33*ucell.latvec.e33);

    EXPECT_NEAR(mag1, mag2, 1e-10);
    EXPECT_NEAR(mag2, mag3, 1e-10);

    // Check FCC structure: e11 should be negative, e13 positive, e12 = 0
    EXPECT_LT(ucell.latvec.e11, 0.0);
    EXPECT_GT(ucell.latvec.e13, 0.0);
    EXPECT_NEAR(ucell.latvec.e12, 0.0, 1e-10);

    // Reset for other tests
    PARAM.input.fixed_ibrav = false;
}

// Test combination of fixed_axes = "volume" and fixed_ibrav
TEST_F(LatticeChangeBasicTest, ChangeLatticeVolumeAndIbrav)
{
    // Initialize variables
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    ucell.lat0 = 10.0;
    ucell.latName = "sc";

    // Set initial lattice (cubic)
    ucell.latvec.e11 = 1.0;
    ucell.latvec.e12 = 0.0;
    ucell.latvec.e13 = 0.0;
    ucell.latvec.e21 = 0.0;
    ucell.latvec.e22 = 1.0;
    ucell.latvec.e23 = 0.0;
    ucell.latvec.e31 = 0.0;
    ucell.latvec.e32 = 0.0;
    ucell.latvec.e33 = 1.0;

    ucell.omega = 1000.0; // Initial volume

    lat[0] = 10.0;
    lat[1] = 0.0;
    lat[2] = 0.0;
    lat[3] = 0.0;
    lat[4] = 10.0;
    lat[5] = 0.0;
    lat[6] = 0.0;
    lat[7] = 0.0;
    lat[8] = 10.0;

    // Apply a move that would change volume and break symmetry
    move[0] = 1.0;
    move[1] = 0.1; // Try to break symmetry
    move[2] = 0.0;
    move[3] = 0.0;
    move[4] = 0.8;
    move[5] = 0.0;
    move[6] = 0.0;
    move[7] = 0.0;
    move[8] = 1.2;

    PARAM.input.fixed_ibrav = true;
    Lattice_Change_Basic::fixed_axes = "volume";

    // Verify remake_cell was not called yet
    EXPECT_FALSE(unitcell::was_remake_cell_called());

    // Call change_lattice method
    Lattice_Change_Basic::change_lattice(ucell, move, lat);

    // Verify remake_cell was called (should be called before volume rescaling)
    EXPECT_TRUE(unitcell::was_remake_cell_called());

    // Check that volume is preserved
    EXPECT_NEAR(ucell.omega, 1000.0, 1e-8);

    // Check that lattice is cubic (all diagonal, equal)
    EXPECT_NEAR(ucell.latvec.e11, ucell.latvec.e22, 1e-10);
    EXPECT_NEAR(ucell.latvec.e22, ucell.latvec.e33, 1e-10);
    EXPECT_NEAR(ucell.latvec.e12, 0.0, 1e-10);
    EXPECT_NEAR(ucell.latvec.e13, 0.0, 1e-10);
    EXPECT_NEAR(ucell.latvec.e21, 0.0, 1e-10);
    EXPECT_NEAR(ucell.latvec.e23, 0.0, 1e-10);
    EXPECT_NEAR(ucell.latvec.e31, 0.0, 1e-10);
    EXPECT_NEAR(ucell.latvec.e32, 0.0, 1e-10);

    // Reset for other tests
    PARAM.input.fixed_ibrav = false;
}

// Test axis constraint with fixed_axes = "a"
TEST_F(LatticeChangeBasicTest, SetupGradientAxisA)
{
    // Initialize variables
    ucell.lc[0] = 0; // First axis fixed
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

    Lattice_Change_Basic::fixed_axes = "a";

    // Call setup_gradient method
    Lattice_Change_Basic::setup_gradient(ucell, lat, grad, stress);

    // Check that gradient for first lattice vector is zero
    EXPECT_DOUBLE_EQ(grad[0], 0.0);
    EXPECT_DOUBLE_EQ(grad[1], 0.0);
    EXPECT_DOUBLE_EQ(grad[2], 0.0);

    // Check that gradients for other vectors are non-zero
    EXPECT_NE(grad[3], 0.0);
    EXPECT_NE(grad[4], 0.0);
    EXPECT_NE(grad[5], 0.0);
    EXPECT_NE(grad[6], 0.0);
    EXPECT_NE(grad[7], 0.0);
    EXPECT_NE(grad[8], 0.0);
}

// Test that fixed axis doesn't move in change_lattice
TEST_F(LatticeChangeBasicTest, ChangeLatticeFixedAxisA)
{
    // Initialize variables
    ucell.lc[0] = 0; // First axis fixed
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    ucell.lat0 = 10.0;

    // Set initial lattice
    ucell.latvec.e11 = 1.0;
    ucell.latvec.e12 = 0.1;
    ucell.latvec.e13 = 0.2;
    ucell.latvec.e21 = 0.0;
    ucell.latvec.e22 = 1.0;
    ucell.latvec.e23 = 0.0;
    ucell.latvec.e31 = 0.0;
    ucell.latvec.e32 = 0.0;
    ucell.latvec.e33 = 1.0;

    // Save initial first lattice vector
    double initial_e11 = ucell.latvec.e11;
    double initial_e12 = ucell.latvec.e12;
    double initial_e13 = ucell.latvec.e13;

    lat[0] = 10.0;
    lat[1] = 1.0;
    lat[2] = 2.0;
    lat[3] = 0.0;
    lat[4] = 10.0;
    lat[5] = 0.0;
    lat[6] = 0.0;
    lat[7] = 0.0;
    lat[8] = 10.0;

    // Try to move all lattice vectors
    move[0] = 1.0;
    move[1] = 0.5;
    move[2] = 0.3;
    move[3] = 0.5;
    move[4] = 1.0;
    move[5] = 0.0;
    move[6] = 0.0;
    move[7] = 0.0;
    move[8] = 1.0;

    Lattice_Change_Basic::fixed_axes = "a";

    // Call change_lattice method
    Lattice_Change_Basic::change_lattice(ucell, move, lat);

    // Check that first lattice vector didn't change
    EXPECT_DOUBLE_EQ(ucell.latvec.e11, initial_e11);
    EXPECT_DOUBLE_EQ(ucell.latvec.e12, initial_e12);
    EXPECT_DOUBLE_EQ(ucell.latvec.e13, initial_e13);

    // Check that other lattice vectors did change
    EXPECT_NE(ucell.latvec.e21, 0.0);
    EXPECT_NE(ucell.latvec.e22, 1.0);
    EXPECT_NE(ucell.latvec.e33, 1.0);
}

// Test that remake_cell is NOT called when fixed_ibrav = false
TEST_F(LatticeChangeBasicTest, ChangeLatticeNoFixedIbrav)
{
    // Initialize variables
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    ucell.lat0 = 10.0;
    ucell.latName = "sc";

    // Set initial lattice
    ucell.latvec.e11 = 1.0;
    ucell.latvec.e12 = 0.01; // Small distortion
    ucell.latvec.e13 = 0.0;
    ucell.latvec.e21 = 0.0;
    ucell.latvec.e22 = 1.0;
    ucell.latvec.e23 = 0.0;
    ucell.latvec.e31 = 0.0;
    ucell.latvec.e32 = 0.0;
    ucell.latvec.e33 = 1.0;

    lat[0] = 10.0;
    lat[1] = 0.1;
    lat[2] = 0.0;
    lat[3] = 0.0;
    lat[4] = 10.0;
    lat[5] = 0.0;
    lat[6] = 0.0;
    lat[7] = 0.0;
    lat[8] = 10.0;

    // Apply a move
    move[0] = 0.1;
    move[1] = 0.0;
    move[2] = 0.0;
    move[3] = 0.0;
    move[4] = 0.1;
    move[5] = 0.0;
    move[6] = 0.0;
    move[7] = 0.0;
    move[8] = 0.1;

    PARAM.input.fixed_ibrav = false; // Explicitly set to false
    Lattice_Change_Basic::fixed_axes = "None";

    // Verify remake_cell was not called yet
    EXPECT_FALSE(unitcell::was_remake_cell_called());

    // Call change_lattice method
    Lattice_Change_Basic::change_lattice(ucell, move, lat);

    // Verify remake_cell was NOT called
    EXPECT_FALSE(unitcell::was_remake_cell_called());

    // Check that distortion remains (e12 should still be non-zero)
    EXPECT_NE(ucell.latvec.e12, 0.0);
}
