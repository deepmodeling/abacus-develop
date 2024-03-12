#include "for_test.h"
#define private public
#include "gtest/gtest.h"
#include "module_relax/relax_old/lattice_change_basic.h"
#include "module_relax/relax_old/lattice_change_cg.h"
/************************************************
 *  unit tests of class Lattice_Change_CG
 ***********************************************/

/**
 * - Tested Functions:
 *   - Lattice_Change_CG::allocate()
 *   - Lattice_Change_CG::start()
 *   - Lattice_Change_CG::setup_cg_grad()
 *   - Lattice_Change_CG::setup_move()
 *   - Lattice_Change_CG::Brent()
 *   - Lattice_Change_CG::f_cal()
 *   - Lattice_Change_CG::third_order()
 */

class LatticeChangeCGTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        // Initialize variables before each test
        Lattice_Change_Basic::dim = 9;
        Lattice_Change_Basic::stress_step = 1;
        Lattice_Change_Basic::update_iter = 5;
        Lattice_Change_Basic::converged = true;
        lc_cg.allocate();
    }

    void TearDown() override
    {
        // Clean up after each test
    }

    Lattice_Change_CG lc_cg;
};

// Test whether the allocate() function can correctly allocate memory space
TEST_F(LatticeChangeCGTest, TestAllocate)
{
    Lattice_Change_Basic::dim = 4;
    lc_cg.allocate();

    // Check if allocated arrays are not empty
    EXPECT_NE(nullptr, lc_cg.lat0);
    EXPECT_NE(nullptr, lc_cg.grad0);
    EXPECT_NE(nullptr, lc_cg.cg_grad0);
    EXPECT_NE(nullptr, lc_cg.move0);
}

// Test if a dimension less than or equal to 0 results in an assertion error
TEST_F(LatticeChangeCGTest, TestAllocateWithZeroDimension)
{
    Lattice_Change_Basic::dim = 0;
    ASSERT_DEATH(lc_cg.allocate(), "");
}

// Check that the arrays are correctly initialized to 0
TEST_F(LatticeChangeCGTest, TestAllocateAndInitialize)
{
    Lattice_Change_Basic::dim = 3;
    lc_cg.allocate();

    // Check that the arrays are correctly initialized to 0
    EXPECT_DOUBLE_EQ(0.0, lc_cg.lat0[0]);
    EXPECT_DOUBLE_EQ(0.0, lc_cg.grad0[1]);
    EXPECT_DOUBLE_EQ(0.0, lc_cg.cg_grad0[2]);
    EXPECT_DOUBLE_EQ(0.0, lc_cg.move0[0]);
}

// Test function start() when converged
TEST_F(LatticeChangeCGTest, TestStartConverged)
{
    // setup data
    UnitCell ucell;
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    ModuleBase::matrix stress(3, 3);
    double etot = 0.0;

    // call function
    GlobalV::ofs_running.open("log");
    lc_cg.start(ucell, stress, etot);
    GlobalV::ofs_running.close();

    // Check output
    std::string expected_output
        = " largest stress is 0, no movement is possible.\n it may converged, otherwise no movement of lattice "
          "parameters is allowed.\n end of lattice optimization\n                              stress_step = 1\n       "
          "                  update iteration = 5\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    ifs.close();
    std::remove("log");
}

// Test function start() sd branch
TEST_F(LatticeChangeCGTest, TestStartSd)
{
    // setup data
    UnitCell ucell;
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    ModuleBase::matrix stress(3, 3);
    stress(0, 0) = 0.01;
    double etot = 0.0;

    // call function
    GlobalV::ofs_running.open("log");
    lc_cg.start(ucell, stress, etot);
    GlobalV::ofs_running.close();

    // Check output
    std::string expected_output = "\n Lattice relaxation is not converged yet (threshold is 0.5)\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    EXPECT_DOUBLE_EQ(Lattice_Change_Basic::lattice_change_ini, 0.01);
    ifs.close();
    std::remove("log");
}

// Test function start() trial branch with goto
TEST_F(LatticeChangeCGTest, TestStartTrialGoto)
{
    // setup data
    UnitCell ucell;
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    ModuleBase::matrix stress(3, 3);
    stress(0, 1) = 0.01;
    double etot = 0.0;

    // call function
    lc_cg.move0[0] = 1.0;
    lc_cg.start(ucell, stress, etot);
    Lattice_Change_Basic::stress_step = 2;
    lc_cg.move0[0] = 10.0;
    GlobalV::ofs_running.open("log");
    lc_cg.start(ucell, stress, etot);
    GlobalV::ofs_running.close();

    // Check output
    std::string expected_output = "\n Lattice relaxation is not converged yet (threshold is 0.5)\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    EXPECT_NEAR(Lattice_Change_Basic::lattice_change_ini, 10.000004999998749, 1e-12);
    ifs.close();
    std::remove("log");
}

// Test function start() trial branch without goto
TEST_F(LatticeChangeCGTest, TestStartTrial)
{
    // setup data
    UnitCell ucell;
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    ModuleBase::matrix stress(3, 3);
    stress(0, 1) = 0.01;
    double etot = 0.0;

    // call function
    lc_cg.start(ucell, stress, etot);
    Lattice_Change_Basic::stress_step = 2;
    GlobalV::ofs_running.open("log");
    lc_cg.start(ucell, stress, etot);
    GlobalV::ofs_running.close();

    // Check output
    std::string expected_output = "\n Lattice relaxation is not converged yet (threshold is 0.5)\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    EXPECT_NEAR(Lattice_Change_Basic::lattice_change_ini, 70.000034999991243, 1e-12);
    ifs.close();
    std::remove("log");
}

// Test function start() no trial branch with goto case 1
TEST_F(LatticeChangeCGTest, TestStartNoTrialGotoCase1)
{
    // setup data
    UnitCell ucell;
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    ModuleBase::matrix stress(3, 3);
    stress(0, 1) = 0.01;
    double etot = 0.0;

    // call function
    lc_cg.start(ucell, stress, etot);
    Lattice_Change_Basic::stress_step = 2;
    lc_cg.start(ucell, stress, etot);
    GlobalV::ofs_running.open("log");
    lc_cg.start(ucell, stress, etot);
    GlobalV::ofs_running.close();

    // Check output
    std::string expected_output = "\n Lattice relaxation is not converged yet (threshold is 0.5)\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    EXPECT_NEAR(Lattice_Change_Basic::lattice_change_ini, 490.00024499993867, 1e-12);
    ifs.close();
    std::remove("log");
}

// Test function start() no trial branch with goto case 2
TEST_F(LatticeChangeCGTest, TestStartNoTrialGotoCase2)
{
    // setup data
    UnitCell ucell;
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    ModuleBase::matrix stress(3, 3);
    stress(0, 1) = 0.01;
    double etot = 0.0;

    // call function
    lc_cg.move0[0] = 0.1;
    lc_cg.start(ucell, stress, etot);
    Lattice_Change_Basic::stress_step = 2;
    lc_cg.start(ucell, stress, etot);
    GlobalV::ofs_running.open("log");
    lc_cg.move0[0] = 0.1;
    stress(0, 1) = 0.0001;
    lc_cg.start(ucell, stress, etot);
    GlobalV::ofs_running.close();

    // Check output
    std::string expected_output = "\n Lattice relaxation is not converged yet (threshold is 0.5)\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    EXPECT_NEAR(Lattice_Change_Basic::lattice_change_ini, 3430.0017149995706, 1e-12);
    ifs.close();
    std::remove("log");
}

// Test function start() no trial branch without goto
TEST_F(LatticeChangeCGTest, TestStartNoTrial)
{
    // setup data
    UnitCell ucell;
    ucell.lc[0] = 1;
    ucell.lc[1] = 1;
    ucell.lc[2] = 1;
    ModuleBase::matrix stress(3, 3);
    stress(0, 1) = 0.01;
    double etot = 0.0;

    // call function
    lc_cg.move0[0] = 1.0;
    lc_cg.start(ucell, stress, etot);
    Lattice_Change_Basic::stress_step = 2;
    lc_cg.move0[0] = 10.0;
    lc_cg.start(ucell, stress, etot);
    GlobalV::ofs_running.open("log");
    lc_cg.start(ucell, stress, etot);
    GlobalV::ofs_running.close();

    // Check output
    std::string expected_output = "\n Lattice relaxation is not converged yet (threshold is 0.5)\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    EXPECT_NEAR(Lattice_Change_Basic::lattice_change_ini, 96040.106328872833, 1e-12);
    ifs.close();
    std::remove("log");
}

// Test function setup_cg_grad() when ncggrad is multiple of 10000
TEST_F(LatticeChangeCGTest, SetupCgGradNcggradIsMultipleOf10000)
{
    double grad[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    double grad0[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    double cggrad[9] = {9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
    double cggrad0[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    int ncggrad = 50000; // multiple of 10000
    int flag = 0;

    lc_cg.setup_cg_grad(grad, grad0, cggrad, cggrad0, ncggrad, flag);

    EXPECT_DOUBLE_EQ(cggrad[0], grad[0]);
    EXPECT_DOUBLE_EQ(cggrad[1], grad[1]);
    EXPECT_DOUBLE_EQ(cggrad[2], grad[2]);
    EXPECT_DOUBLE_EQ(cggrad[3], grad[3]);
    EXPECT_DOUBLE_EQ(cggrad[4], grad[4]);
    EXPECT_DOUBLE_EQ(cggrad[5], grad[5]);
    EXPECT_DOUBLE_EQ(cggrad[6], grad[6]);
    EXPECT_DOUBLE_EQ(cggrad[7], grad[7]);
    EXPECT_DOUBLE_EQ(cggrad[8], grad[8]);
}

// Test function setup_cg_grad() when ncggrad is not multiple of 10000, gamma1 < 0.5
TEST_F(LatticeChangeCGTest, SetupCgGradNcggradIsNotMultipleOf10000Case1)
{
    double grad[9] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double grad0[9] = {4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double cggrad[9] = {4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double cggrad0[9] = {4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int ncggrad = 100;
    int flag = 0;

    lc_cg.setup_cg_grad(grad, grad0, cggrad, cggrad0, ncggrad, flag);

    EXPECT_DOUBLE_EQ(cggrad[0], 1.25);
    EXPECT_DOUBLE_EQ(cggrad[1], 0.0);
    EXPECT_DOUBLE_EQ(cggrad[2], 0.0);
    EXPECT_DOUBLE_EQ(cggrad[3], 0.0);
    EXPECT_DOUBLE_EQ(cggrad[4], 0.0);
    EXPECT_DOUBLE_EQ(cggrad[5], 0.0);
    EXPECT_DOUBLE_EQ(cggrad[6], 0.0);
    EXPECT_DOUBLE_EQ(cggrad[7], 0.0);
    EXPECT_DOUBLE_EQ(cggrad[8], 0.0);
}

// Test function setup_cg_grad() when ncggrad is not multiple of 10000, gamma1 >= 0.5
TEST_F(LatticeChangeCGTest, SetupCgGradNcggradIsNotMultipleOf10000Case2)
{
    double grad[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    double grad0[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    double cggrad[9] = {9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
    double cggrad0[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    int ncggrad = 100;
    int flag = 0;

    lc_cg.setup_cg_grad(grad, grad0, cggrad, cggrad0, ncggrad, flag);

    EXPECT_DOUBLE_EQ(cggrad[0], grad[0]);
    EXPECT_DOUBLE_EQ(cggrad[1], grad[1]);
    EXPECT_DOUBLE_EQ(cggrad[2], grad[2]);
    EXPECT_DOUBLE_EQ(cggrad[3], grad[3]);
    EXPECT_DOUBLE_EQ(cggrad[4], grad[4]);
    EXPECT_DOUBLE_EQ(cggrad[5], grad[5]);
    EXPECT_DOUBLE_EQ(cggrad[6], grad[6]);
    EXPECT_DOUBLE_EQ(cggrad[7], grad[7]);
    EXPECT_DOUBLE_EQ(cggrad[8], grad[8]);
}

// Test function third_order() case 1
TEST_F(LatticeChangeCGTest, ThirdOrderCase1)
{
    double e0 = 1.0;
    double e1 = 1.0;
    double fa = 10.0;
    double fb = -9.99;
    double x = 1.0;
    double bestX = -1.0; // arbitrary initial value

    lc_cg.third_order(e0, e1, fa, fb, x, bestX);

    EXPECT_DOUBLE_EQ(bestX, x * fb / (fa - fb));
}

// Test function third_order() case 2
TEST_F(LatticeChangeCGTest, ThirdOrderCase2)
{
    double e0 = 1.0;
    double e1 = 1.0;
    double fa = -10.0;
    double fb = 9.9;
    double x = 1.0;
    double bestX = -1.0; // arbitrary initial value

    lc_cg.third_order(e0, e1, fa, fb, x, bestX);

    EXPECT_DOUBLE_EQ(bestX, x * fb / (fa - fb));
}

// Test function third_order() case 3
TEST_F(LatticeChangeCGTest, ThirdOrderCase3)
{
    double e0 = 1.0;
    double e1 = 1.0;
    double fa = 10.0;
    double fb = -10.1;
    double x = 1.0;
    double bestX = -1.0; // arbitrary initial value

    lc_cg.third_order(e0, e1, fa, fb, x, bestX);

    EXPECT_DOUBLE_EQ(bestX, x * fb / (fa - fb));
}

// Test function Brent() case 1
TEST_F(LatticeChangeCGTest, BrentCase1)
{
    double fa = 2.0;
    double fb = 1.0;
    double fc = 1.0;
    double xa = -3.0;
    double xb = 2.0;
    double xc = 1.0;
    double best_x = 0.0;
    double xpt = 0.0;

    lc_cg.Brent(fa, fb, fc, xa, xb, xc, best_x, xpt);

    EXPECT_DOUBLE_EQ(fa, 2.0);
    EXPECT_DOUBLE_EQ(fb, 1.0);
    EXPECT_DOUBLE_EQ(fc, 1.0);
    EXPECT_DOUBLE_EQ(xa, -3.0);
    EXPECT_DOUBLE_EQ(xb, 1.0);
    EXPECT_DOUBLE_EQ(xc, 4.0);
    EXPECT_DOUBLE_EQ(best_x, 4.0);
    EXPECT_DOUBLE_EQ(xpt, 4.0);
}

// Test function Brent() case 2
TEST_F(LatticeChangeCGTest, BrentCase2)
{
    double fa = -2.0;
    double fb = 3.0;
    double fc = -4.0;
    double xa = 1.0;
    double xb = 2.0;
    double xc = 3.0;
    double best_x = 0.0;
    double xpt = 0.0;

    lc_cg.Brent(fa, fb, fc, xa, xb, xc, best_x, xpt);

    EXPECT_DOUBLE_EQ(fa, -4.0);
    EXPECT_DOUBLE_EQ(fb, 3.0);
    EXPECT_DOUBLE_EQ(fc, -4.0);
    EXPECT_DOUBLE_EQ(xa, 3.0);
    EXPECT_DOUBLE_EQ(xb, 2.0);
    EXPECT_NEAR(xc, 1.2046663545568725, 1e-12);
    EXPECT_NEAR(best_x, 1.2046663545568725, 1e-12);
    EXPECT_NEAR(xpt, 1.2046663545568725, 1e-12);
}

// Test function Brent() case 3
TEST_F(LatticeChangeCGTest, BrentCase3)
{
    double fa = 1.0;
    double fb = -3.0;
    double fc = -4.0;
    double xa = 3.0;
    double xb = 2.0;
    double xc = 1.0;
    double best_x = 0.0;
    double xpt = 0.0;

    lc_cg.Brent(fa, fb, fc, xa, xb, xc, best_x, xpt);

    EXPECT_DOUBLE_EQ(fa, 1.0);
    EXPECT_DOUBLE_EQ(fb, -4.0);
    EXPECT_DOUBLE_EQ(fc, -4.0);
    EXPECT_DOUBLE_EQ(xa, 3.0);
    EXPECT_DOUBLE_EQ(xb, 1.0);
    EXPECT_NEAR(xc, 2.8081429669660172, 1e-12);
    EXPECT_NEAR(best_x, 2.8081429669660172, 1e-12);
    EXPECT_NEAR(xpt, 2.8081429669660172, 1e-12);
}

// Test function Brent() case 4
TEST_F(LatticeChangeCGTest, BrentCase4)
{
    double fa = 2.0;
    double fb = -3.0;
    double fc = 4.0;
    double xa = 0.0;
    double xb = 2.0;
    double xc = 1.0;
    double best_x = 0.0;
    double xpt = 0.0;

    lc_cg.Brent(fa, fb, fc, xa, xb, xc, best_x, xpt);

    EXPECT_DOUBLE_EQ(fa, 4.0);
    EXPECT_DOUBLE_EQ(fb, -3.0);
    EXPECT_DOUBLE_EQ(fc, 4.0);
    EXPECT_DOUBLE_EQ(xa, 1.0);
    EXPECT_DOUBLE_EQ(xb, 2.0);
    EXPECT_DOUBLE_EQ(xc, 2.0);
    EXPECT_DOUBLE_EQ(best_x, 2.0);
    EXPECT_DOUBLE_EQ(xpt, 2.0);
}

// Test function f_cal()
TEST_F(LatticeChangeCGTest, Fcal)
{
    Lattice_Change_Basic::dim = 9;
    double g0[9] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double g1[9] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double f_value;

    lc_cg.f_cal(g0, g1, Lattice_Change_Basic::dim, f_value);

    EXPECT_DOUBLE_EQ(f_value, 3.0);
}

// Test function setup_move()
TEST_F(LatticeChangeCGTest, SetupMove)
{
    Lattice_Change_Basic::dim = 9;
    double trust_radius = 1.0;
    double cg_gradn[9] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double move[9];

    lc_cg.setup_move(move, cg_gradn, trust_radius);

    EXPECT_DOUBLE_EQ(move[0], -1.0);
    EXPECT_DOUBLE_EQ(move[1], -1.0);
    EXPECT_DOUBLE_EQ(move[2], -1.0);
    EXPECT_DOUBLE_EQ(move[3], -1.0);
    EXPECT_DOUBLE_EQ(move[4], -1.0);
    EXPECT_DOUBLE_EQ(move[5], -1.0);
    EXPECT_DOUBLE_EQ(move[6], -1.0);
    EXPECT_DOUBLE_EQ(move[7], -1.0);
    EXPECT_DOUBLE_EQ(move[8], -1.0);
}
