#include "for_test.h"
#define private public
#include "gtest/gtest.h"
#include "module_relax/relax_old/ions_move_basic.h"
#include "module_relax/relax_old/ions_move_cg.h"
/************************************************
 *  unit tests of class Ions_Move_CG
 ***********************************************/

/**
 * - Tested Functions:
 *   - Ions_Move_CG::allocate()
 *   - Ions_Move_CG::start()
 *   - Ions_Move_CG::setup_cg_grad()
 *   - Ions_Move_CG::setup_move()
 *   - Ions_Move_CG::Brent()
 *   - Ions_Move_CG::f_cal()
 *   - Ions_Move_CG::third_order()
 */

class IonsMoveCGTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        // Initialize variables before each test
        Ions_Move_Basic::dim = 6;
        Ions_Move_Basic::update_iter = 5;
        im_cg.allocate();
    }

    void TearDown() override
    {
        // Clean up after each test
    }

    Ions_Move_CG im_cg;
};

// Test whether the allocate() function can correctly allocate memory space
TEST_F(IonsMoveCGTest, TestAllocate)
{
    Ions_Move_Basic::dim = 4;
    im_cg.allocate();

    // Check if allocated arrays are not empty
    EXPECT_NE(nullptr, im_cg.pos0);
    EXPECT_NE(nullptr, im_cg.grad0);
    EXPECT_NE(nullptr, im_cg.cg_grad0);
    EXPECT_NE(nullptr, im_cg.move0);
}

// Test if a dimension less than or equal to 0 results in an assertion error
TEST_F(IonsMoveCGTest, TestAllocateWithZeroDimension)
{
    Ions_Move_Basic::dim = 0;
    ASSERT_DEATH(im_cg.allocate(), "");
}

// Check that the arrays are correctly initialized to 0
TEST_F(IonsMoveCGTest, TestAllocateAndInitialize)
{
    Ions_Move_Basic::dim = 3;
    im_cg.allocate();

    // Check that the arrays are correctly initialized to 0
    EXPECT_DOUBLE_EQ(0.0, im_cg.pos0[0]);
    EXPECT_DOUBLE_EQ(0.0, im_cg.grad0[1]);
    EXPECT_DOUBLE_EQ(0.0, im_cg.cg_grad0[2]);
    EXPECT_DOUBLE_EQ(0.0, im_cg.move0[0]);
}

// Test function start() when converged
TEST_F(IonsMoveCGTest, TestStartConverged)
{
    // setup data
    Ions_Move_Basic::istep = 1;
    Ions_Move_Basic::converged = true;
    UnitCell ucell;
    ModuleBase::matrix force(2, 3);
    double etot = 0.0;

    // call function
    GlobalV::ofs_running.open("log");
    im_cg.start(ucell, force, etot);
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
}

// Test function start() sd branch
TEST_F(IonsMoveCGTest, TestStartSd)
{
    // setup data
    Ions_Move_Basic::istep = 1;
    Ions_Move_Basic::converged = false;
    GlobalV::RELAX_METHOD = "cg_bfgs";
    Ions_Move_CG::RELAX_CG_THR = 100.0;
    UnitCell ucell;
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 0.01;
    double etot = 0.0;

    // call function
    GlobalV::ofs_running.open("log");
    im_cg.start(ucell, force, etot);
    GlobalV::ofs_running.close();

    // Check output
    std::string expected_output = "\n Ion relaxation is not converged yet (threshold is 0.0257111)\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_EQ(expected_output, output);
    EXPECT_EQ(Ions_Move_Basic::converged, false);
    EXPECT_EQ(Ions_Move_Basic::update_iter, 5);
    EXPECT_EQ(GlobalV::RELAX_METHOD, "bfgs");
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.01);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::best_xxx, -1.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::relax_bfgs_init, 1.0);
}

// Test function start() trial branch with goto
TEST_F(IonsMoveCGTest, TestStartTrialGoto)
{
    // setup data
    Ions_Move_Basic::istep = 1;
    Ions_Move_Basic::converged = false;
    UnitCell ucell;
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 0.1;
    double etot = 0.0;

    // call function
    im_cg.move0[0] = 1.0;
    im_cg.start(ucell, force, etot);
    Ions_Move_Basic::istep = 2;
    im_cg.move0[0] = 10.0;
    force(0, 0) = 0.001;
    GlobalV::ofs_running.open("log");
    im_cg.start(ucell, force, etot);
    GlobalV::ofs_running.close();

    // Check output
    std::string expected_output = "\n Ion relaxation is not converged yet (threshold is 0.0257111)\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_EQ(expected_output, output);
    EXPECT_EQ(Ions_Move_Basic::converged, false);
    EXPECT_EQ(Ions_Move_Basic::update_iter, 5);
    EXPECT_EQ(GlobalV::RELAX_METHOD, "bfgs");
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.001);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::best_xxx, -1.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::relax_bfgs_init, 10.0);
}

// Test function start() trial branch without goto
TEST_F(IonsMoveCGTest, TestStartTrial)
{
    // setup data
    Ions_Move_Basic::istep = 1;
    Ions_Move_Basic::converged = false;
    UnitCell ucell;
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 0.01;
    double etot = 0.0;

    // call function
    im_cg.move0[0] = 1.0;
    im_cg.start(ucell, force, etot);
    Ions_Move_Basic::istep = 2;
    im_cg.move0[0] = 10.0;
    GlobalV::ofs_running.open("log");
    im_cg.start(ucell, force, etot);
    GlobalV::ofs_running.close();

    // Check output
    std::string expected_output = "\n Ion relaxation is not converged yet (threshold is 0.0257111)\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_EQ(expected_output, output);
    EXPECT_EQ(Ions_Move_Basic::converged, false);
    EXPECT_EQ(Ions_Move_Basic::update_iter, 5);
    EXPECT_EQ(GlobalV::RELAX_METHOD, "bfgs");
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.01);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::best_xxx, -1.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::relax_bfgs_init, 70.0);
}

// Test function start() no trial branch with goto case 1
TEST_F(IonsMoveCGTest, TestStartNoTrialGotoCase1)
{
    // setup data
    Ions_Move_Basic::istep = 1;
    Ions_Move_Basic::converged = false;
    UnitCell ucell;
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 0.1;
    double etot = 0.0;

    // call function
    im_cg.move0[0] = 1.0;
    im_cg.start(ucell, force, etot);
    Ions_Move_Basic::istep = 2;
    im_cg.start(ucell, force, etot);
    im_cg.move0[0] = 1.0;
    force(0, 0) = 0.001;
    GlobalV::ofs_running.open("log");
    im_cg.start(ucell, force, etot);
    GlobalV::ofs_running.close();

    // Check output
    std::string expected_output = "\n Ion relaxation is not converged yet (threshold is 0.0257111)\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_EQ(expected_output, output);
    EXPECT_EQ(Ions_Move_Basic::converged, false);
    EXPECT_EQ(Ions_Move_Basic::update_iter, 5);
    EXPECT_EQ(GlobalV::RELAX_METHOD, "bfgs");
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.001);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::best_xxx, -1.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::relax_bfgs_init, 490.0);
}

// Test function start() no trial branch with goto case 2
TEST_F(IonsMoveCGTest, TestStartNoTrialGotoCase2)
{
    // setup data
    Ions_Move_Basic::istep = 1;
    Ions_Move_Basic::converged = false;
    UnitCell ucell;
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 0.01;
    double etot = 0.0;

    // call function
    im_cg.move0[0] = 1.0;
    im_cg.start(ucell, force, etot);
    Ions_Move_Basic::istep = 2;
    im_cg.move0[0] = 10.0;
    im_cg.start(ucell, force, etot);
    GlobalV::ofs_running.open("log");
    im_cg.start(ucell, force, etot);
    GlobalV::ofs_running.close();

    // Check output
    std::string expected_output = "\n Ion relaxation is not converged yet (threshold is 0.0257111)\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_EQ(expected_output, output);
    EXPECT_EQ(Ions_Move_Basic::converged, false);
    EXPECT_EQ(Ions_Move_Basic::update_iter, 5);
    EXPECT_EQ(GlobalV::RELAX_METHOD, "bfgs");
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.01);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::best_xxx, -1.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::relax_bfgs_init, 70.0);
}

// Test function start() no trial branch without goto
TEST_F(IonsMoveCGTest, TestStartNoTrial)
{
    // setup data
    Ions_Move_Basic::istep = 1;
    Ions_Move_Basic::converged = false;
    UnitCell ucell;
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 0.01;
    double etot = 0.0;

    // call function
    im_cg.move0[0] = 1.0;
    im_cg.start(ucell, force, etot);
    Ions_Move_Basic::istep = 2;
    im_cg.move0[0] = 1.0;
    force(0, 0) = 0.001;
    im_cg.start(ucell, force, etot);
    GlobalV::ofs_running.open("log");
    im_cg.start(ucell, force, etot);
    GlobalV::ofs_running.close();

    // Check output
    std::string expected_output = "\n Ion relaxation is not converged yet (threshold is 0.0257111)\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_EQ(expected_output, output);
    EXPECT_EQ(Ions_Move_Basic::converged, false);
    EXPECT_EQ(Ions_Move_Basic::update_iter, 5);
    EXPECT_EQ(GlobalV::RELAX_METHOD, "bfgs");
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.001);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::best_xxx, -1.0);
    EXPECT_NEAR(Ions_Move_Basic::relax_bfgs_init, 1.2345679012345678, 1e-12);
}

// Test function setup_cg_grad() when ncggrad is multiple of 10000
TEST_F(IonsMoveCGTest, SetupCgGradNcggradIsMultipleOf10000)
{
    double grad[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double grad0[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double cggrad[6] = {9.0, 8.0, 7.0, 6.0, 5.0, 4.0};
    double cggrad0[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    int ncggrad = 50000; // multiple of 10000
    int flag = 0;

    im_cg.setup_cg_grad(grad, grad0, cggrad, cggrad0, ncggrad, flag);

    EXPECT_DOUBLE_EQ(cggrad[0], grad[0]);
    EXPECT_DOUBLE_EQ(cggrad[1], grad[1]);
    EXPECT_DOUBLE_EQ(cggrad[2], grad[2]);
    EXPECT_DOUBLE_EQ(cggrad[3], grad[3]);
    EXPECT_DOUBLE_EQ(cggrad[4], grad[4]);
    EXPECT_DOUBLE_EQ(cggrad[5], grad[5]);
}

// Test function setup_cg_grad() when ncggrad is not multiple of 10000, gamma1 < 0.5
TEST_F(IonsMoveCGTest, SetupCgGradNcggradIsNotMultipleOf10000Case1)
{
    double grad[6] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double grad0[6] = {4.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double cggrad[6] = {4.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double cggrad0[6] = {4.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int ncggrad = 100;
    int flag = 0;

    im_cg.setup_cg_grad(grad, grad0, cggrad, cggrad0, ncggrad, flag);

    EXPECT_DOUBLE_EQ(cggrad[0], 1.25);
    EXPECT_DOUBLE_EQ(cggrad[1], 0.0);
    EXPECT_DOUBLE_EQ(cggrad[2], 0.0);
    EXPECT_DOUBLE_EQ(cggrad[3], 0.0);
    EXPECT_DOUBLE_EQ(cggrad[4], 0.0);
    EXPECT_DOUBLE_EQ(cggrad[5], 0.0);
}

// Test function setup_cg_grad() when ncggrad is not multiple of 10000, gamma1 >= 0.5
TEST_F(IonsMoveCGTest, SetupCgGradNcggradIsNotMultipleOf10000Case2)
{
    double grad[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double grad0[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double cggrad[6] = {9.0, 8.0, 7.0, 6.0, 5.0, 4.0};
    double cggrad0[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    int ncggrad = 100;
    int flag = 0;

    im_cg.setup_cg_grad(grad, grad0, cggrad, cggrad0, ncggrad, flag);

    EXPECT_DOUBLE_EQ(cggrad[0], grad[0]);
    EXPECT_DOUBLE_EQ(cggrad[1], grad[1]);
    EXPECT_DOUBLE_EQ(cggrad[2], grad[2]);
    EXPECT_DOUBLE_EQ(cggrad[3], grad[3]);
    EXPECT_DOUBLE_EQ(cggrad[4], grad[4]);
    EXPECT_DOUBLE_EQ(cggrad[5], grad[5]);
}

// Test function third_order() case 1
TEST_F(IonsMoveCGTest, ThirdOrderCase1)
{
    double e0 = 1.0;
    double e1 = 1.0;
    double fa = 10.0;
    double fb = -9.99;
    double x = 1.0;
    double bestX = -1.0; // arbitrary initial value

    im_cg.third_order(e0, e1, fa, fb, x, bestX);

    EXPECT_DOUBLE_EQ(bestX, x * fb / (fa - fb));
}

// Test function third_order() case 2
TEST_F(IonsMoveCGTest, ThirdOrderCase2)
{
    double e0 = 1.0;
    double e1 = 1.0;
    double fa = -10.0;
    double fb = 9.9;
    double x = 1.0;
    double bestX = -1.0; // arbitrary initial value

    im_cg.third_order(e0, e1, fa, fb, x, bestX);

    EXPECT_DOUBLE_EQ(bestX, x * fb / (fa - fb));
}

// Test function third_order() case 3
TEST_F(IonsMoveCGTest, ThirdOrderCase3)
{
    double e0 = 1.0;
    double e1 = 1.0;
    double fa = 10.0;
    double fb = -10.1;
    double x = 1.0;
    double bestX = -1.0; // arbitrary initial value

    im_cg.third_order(e0, e1, fa, fb, x, bestX);

    EXPECT_DOUBLE_EQ(bestX, x * fb / (fa - fb));
}

// Test function Brent() case 1
TEST_F(IonsMoveCGTest, BrentCase1)
{
    double fa = 2.0;
    double fb = 1.0;
    double fc = 1.0;
    double xa = -3.0;
    double xb = 2.0;
    double xc = 1.0;
    double best_x = 0.0;
    double xpt = 0.0;

    im_cg.Brent(fa, fb, fc, xa, xb, xc, best_x, xpt);

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
TEST_F(IonsMoveCGTest, BrentCase2)
{
    double fa = -2.0;
    double fb = 3.0;
    double fc = -4.0;
    double xa = 1.0;
    double xb = 2.0;
    double xc = 3.0;
    double best_x = 0.0;
    double xpt = 0.0;

    im_cg.Brent(fa, fb, fc, xa, xb, xc, best_x, xpt);

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
TEST_F(IonsMoveCGTest, BrentCase3)
{
    double fa = 1.0;
    double fb = -3.0;
    double fc = -4.0;
    double xa = 3.0;
    double xb = 2.0;
    double xc = 1.0;
    double best_x = 0.0;
    double xpt = 0.0;

    im_cg.Brent(fa, fb, fc, xa, xb, xc, best_x, xpt);

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
TEST_F(IonsMoveCGTest, BrentCase4)
{
    double fa = 2.0;
    double fb = -3.0;
    double fc = 4.0;
    double xa = 0.0;
    double xb = 2.0;
    double xc = 1.0;
    double best_x = 0.0;
    double xpt = 0.0;

    im_cg.Brent(fa, fb, fc, xa, xb, xc, best_x, xpt);

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
TEST_F(IonsMoveCGTest, Fcal)
{
    Ions_Move_Basic::dim = 9;
    double g0[9] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double g1[9] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double f_value;

    im_cg.f_cal(g0, g1, Ions_Move_Basic::dim, f_value);

    EXPECT_DOUBLE_EQ(f_value, 3.0);
}

// Test function setup_move()
TEST_F(IonsMoveCGTest, SetupMove)
{
    Ions_Move_Basic::dim = 9;
    double trust_radius = 1.0;
    double cg_gradn[9] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double move[9];

    im_cg.setup_move(move, cg_gradn, trust_radius);

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
