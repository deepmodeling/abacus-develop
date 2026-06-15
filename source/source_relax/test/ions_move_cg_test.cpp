#include <regex>
#include "for_test.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#define private public
#include "source_io/module_parameter/parameter.h"
#include "source_relax/ions_move_basic.h"
#include "source_relax/ions_move_cg.h"
#undef private

/************************************************
 *  unit tests of class Ions_Move_CG
 ***********************************************/

class IonsMoveCGTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        // Initialize variables before each test
        Ions_Move_Basic::dim = 6;
        update_iter = 5;
        im_cg.allocate(Ions_Move_Basic::dim);
        PARAM.input.force_thr = 0.001;

        // ban the 'cout' 
        // mohan add 2025-05-02
        std::cout.rdbuf(NULL);
    }

    void TearDown() override
    {
        // Clean up after each test
    }
    void setupucell(UnitCell& ucell)
    {
        for (int it = 0; it < ucell.ntype; it++)
        {
            Atom* atom = &ucell.atoms[it];
            atom->label="test";
            for (int ia = 0; ia < atom->na; ia++)
            {
                atom->mag[ia]= 1;
                for (int ik = 0; ik < 3; ++ik)
                {
                    atom->tau[ia][ik] = 1;
                    atom->mbl[ia][ik] = 1;
                    atom->vel[ia][ik] = 1;
                }
            }
        }
        ucell.lat.GT.Zero();
    }
    Ions_Move_CG im_cg;
    int update_iter;
};

// Test whether the allocate() function can correctly allocate memory space
TEST_F(IonsMoveCGTest, TestAllocate)
{
    const int dim = 4;
    im_cg.allocate(dim);

    // Check if allocated vectors are not empty
    EXPECT_EQ(im_cg.pos0.size(), 4U);
    EXPECT_EQ(im_cg.grad0.size(), 4U);
    EXPECT_EQ(im_cg.cg_grad0.size(), 4U);
    EXPECT_EQ(im_cg.move0.size(), 4U);
}

// Test if a dimension less than or equal to 0 results in an assertion error
TEST_F(IonsMoveCGTest, TestAllocateWithZeroDimension)
{
    const int dim = 0;
    ASSERT_DEATH(im_cg.allocate(dim), "");
}

// Check that the arrays are correctly initialized to 0
TEST_F(IonsMoveCGTest, TestAllocateAndInitialize)
{
    const int dim = 3;
    im_cg.allocate(dim);

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
    const int istep = 1;
    UnitCell ucell;
    setupucell(ucell);
    ModuleBase::matrix force(2, 3);
    double etot = 0.0;
    std::vector<double> etot_info(2, 0.0);
    std::vector<std::string> relax_method = {"cg", "1"};

    // call function
    std::ofstream ofs("TestStartConverged.log");
    im_cg.start(ucell, force, etot, istep, update_iter, ofs, etot_info, relax_method);
    ofs.close();

    // Check output
    std::string expected_output = "\n Largest force is 0 eV/Angstrom while threshold is -1 eV/Angstrom\n"
                                  " largest force is 0, no movement is possible.\n it may converged, otherwise no "
                                  "movement of atom is allowed.\n end of geometry optimization\n                       "
                                  "             istep = 1\n                         update iteration = 5\n";
    std::ifstream ifs("TestStartConverged.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("TestStartConverged.log");


    std::regex pattern(R"(==> .*::.*\t[\d\.]+ GB\t\d+ s\n )");
    output = std::regex_replace(output, pattern, "");
    EXPECT_THAT(output, testing::HasSubstr(expected_output));
    EXPECT_EQ(update_iter, 5);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.0);
}

// Test function start() sd branch
TEST_F(IonsMoveCGTest, TestStartSd)
{
    // setup data
    const int istep = 1;
    std::vector<std::string> relax_method = {"cg_bfgs", "1"};
    Ions_Move_CG::RELAX_CG_THR = 100.0;
    UnitCell ucell;
    setupucell(ucell);
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 0.01;
    double etot = 0.0;
    std::vector<double> etot_info(2, 0.0);

    // call function
    std::ofstream ofs("TestStartSd.log");
    im_cg.start(ucell, force, etot, istep, update_iter, ofs, etot_info, relax_method);
    ofs.close();

    // Check output
    std::string expected_output = "\n Largest force is 0.257111 eV/Angstrom while threshold is -1 eV/Angstrom\n\n"
                                  " Ion relaxation is not converged yet (threshold is 0.0257111)\n";
    std::ifstream ifs("TestStartSd.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("TestStartSd.log"); // mohan

    EXPECT_THAT(output, testing::HasSubstr(expected_output));
    EXPECT_EQ(update_iter, 5);
    EXPECT_EQ(relax_method[0], "bfgs");
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.01);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::best_xxx, -1.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::relax_bfgs_init, 1.0);
}

// Test function start() trial branch with goto
TEST_F(IonsMoveCGTest, TestStartTrialGoto)
{
    // setup data
    const int istep = 1;
    Ions_Move_CG::RELAX_CG_THR = 100.0;
    UnitCell ucell;
    setupucell(ucell);
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 0.1;
    double etot = 0.0;
    std::vector<double> etot_info(2, 0.0);
    std::vector<std::string> relax_method = {"cg_bfgs", "1"};

    // call function
    im_cg.move0[0] = 1.0;
    std::ofstream ofs1("TestStartTrialGoto_temp1.log");
    im_cg.start(ucell, force, etot, istep, update_iter, ofs1, etot_info, relax_method);
    ofs1.close();
    std::remove("TestStartTrialGoto_temp1.log");
    int istep_2 = 2;
    im_cg.move0[0] = 10.0;
    force(0, 0) = 0.001;
    relax_method = {"cg_bfgs", "1"};
    std::ofstream ofs("TestStartTrialGoto.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs, etot_info, relax_method);
    ofs.close();

    // Check output
    std::string expected_output = "\n Largest force is 0.0257111 eV/Angstrom while threshold is -1 eV/Angstrom\n\n"
                                  " Ion relaxation is not converged yet (threshold is 0.0257111)\n";
    std::ifstream ifs("TestStartTrialGoto.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("TestStartTrialGoto.log");

    EXPECT_THAT(output, testing::HasSubstr(expected_output));
    EXPECT_EQ(update_iter, 5);
    EXPECT_EQ(relax_method[0], "bfgs");
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.001);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::best_xxx, 10.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::relax_bfgs_init, 10.0);
}

// Test function start() trial branch without goto
TEST_F(IonsMoveCGTest, TestStartTrial)
{
    // setup data
    const int istep = 1;
    UnitCell ucell;
    setupucell(ucell);
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 0.01;
    double etot = 0.0;
    std::vector<double> etot_info(2, 0.0);
    std::vector<std::string> relax_method = {"cg_bfgs", "1"};

    // call function
    im_cg.move0[0] = 1.0;
    std::ofstream ofs1("TestStartTrial_temp1.log");
    im_cg.start(ucell, force, etot, istep, update_iter, ofs1, etot_info, relax_method);
    ofs1.close();
    std::remove("TestStartTrial_temp1.log");
    int istep_2 = 2;
    im_cg.move0[0] = 10.0;
    std::ofstream ofs("TestStartTrial.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs, etot_info, relax_method);
    ofs.close();

    // Check output
    std::string expected_output = "\n Largest force is 0.257111 eV/Angstrom while threshold is -1 eV/Angstrom\n\n"
                                  " Ion relaxation is not converged yet (threshold is 0.0257111)\n";
    std::ifstream ifs("TestStartTrial.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("TestStartTrial.log");

    EXPECT_THAT(output, testing::HasSubstr(expected_output));
    EXPECT_EQ(update_iter, 5);
    EXPECT_EQ(relax_method[0], "bfgs");
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.01);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::best_xxx, 10.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::relax_bfgs_init, 70.0);
}

// Test function start() no trial branch with goto case 1
TEST_F(IonsMoveCGTest, TestStartNoTrialGotoCase1)
{
    // setup data
    const int istep = 1;
    Ions_Move_CG::RELAX_CG_THR = 100.0;
    UnitCell ucell;
    setupucell(ucell);
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 0.1;
    double etot = 0.0;
    std::vector<double> etot_info(2, 0.0);
    std::vector<std::string> relax_method = {"cg_bfgs", "1"};

    // call function
    im_cg.move0[0] = 1.0;
    std::ofstream ofs1("TestStartNoTrialGotoCase1_temp1.log");
    im_cg.start(ucell, force, etot, istep, update_iter, ofs1, etot_info, relax_method);
    ofs1.close();
    std::remove("TestStartNoTrialGotoCase1_temp1.log");
    int istep_2 = 2;
    std::ofstream ofs2("TestStartNoTrialGotoCase1_temp2.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs2, etot_info, relax_method);
    ofs2.close();
    std::remove("TestStartNoTrialGotoCase1_temp2.log");
    im_cg.move0[0] = 1.0;
    force(0, 0) = 0.001;
    relax_method = {"cg_bfgs", "1"};
    std::ofstream ofs("TestStartNoTrialGotoCase1.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs, etot_info, relax_method);
    ofs.close();

    // Check output
    std::string expected_output = "\n Largest force is 0.0257111 eV/Angstrom while threshold is -1 eV/Angstrom\n\n"
                                  " Ion relaxation is not converged yet (threshold is 0.0257111)\n";
    std::ifstream ifs("TestStartNoTrialGotoCase1.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("TestStartNoTrialGotoCase1.log");

    EXPECT_THAT(output, testing::HasSubstr(expected_output));
    EXPECT_EQ(update_iter, 5);
    EXPECT_EQ(relax_method[0], "bfgs");
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.001);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::best_xxx, 490.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::relax_bfgs_init, 490.0);
}

// Test function start() no trial branch with goto case 2
TEST_F(IonsMoveCGTest, TestStartNoTrialGotoCase2)
{
    // setup data
    const int istep = 1;
    Ions_Move_CG::RELAX_CG_THR = 100.0;
    UnitCell ucell;
    setupucell(ucell);
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 0.01;
    double etot = 0.0;
    std::vector<double> etot_info(2, 0.0);
    std::vector<std::string> relax_method = {"cg_bfgs", "1"};
    Ions_Move_Basic::best_xxx = 1.0;
    Ions_Move_Basic::relax_bfgs_init = 1.0;

    // call function
    im_cg.move0[0] = 1.0;
    std::ofstream ofs1("TestStartNoTrialGotoCase2_temp1.log");
    im_cg.start(ucell, force, etot, istep, update_iter, ofs1, etot_info, relax_method);
    ofs1.close();
    std::remove("TestStartNoTrialGotoCase2_temp1.log");
    int istep_2 = 2;
    im_cg.move0[0] = 10.0;
    std::ofstream ofs2("TestStartNoTrialGotoCase2_temp2.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs2, etot_info, relax_method);
    ofs2.close();
    std::remove("TestStartNoTrialGotoCase2_temp2.log");
    relax_method = {"cg_bfgs", "1"};
    std::ofstream ofs("TestStartNoTrialGotoCase2.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs, etot_info, relax_method);
    ofs.close();

    // Check output
    std::string expected_output = "\n Largest force is 0.257111 eV/Angstrom while threshold is -1 eV/Angstrom\n\n"
                                  " Ion relaxation is not converged yet (threshold is 0.0257111)\n";
    std::ifstream ifs("TestStartNoTrialGotoCase2.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("TestStartNoTrialGotoCase2.log");

    EXPECT_THAT(output, testing::HasSubstr(expected_output));
    EXPECT_EQ(update_iter, 5);
    EXPECT_EQ(relax_method[0], "bfgs");
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.01);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::best_xxx, 70.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::relax_bfgs_init, 70.0);
}

// Test function start() no trial branch without goto
TEST_F(IonsMoveCGTest, TestStartNoTrial)
{
    // setup data
    const int istep = 1;
    Ions_Move_CG::RELAX_CG_THR = 100.0;
    UnitCell ucell;
    setupucell(ucell);
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 0.01;
    double etot = 0.0;
    std::vector<double> etot_info(2, 0.0);
    std::vector<std::string> relax_method = {"cg_bfgs", "1"};
    Ions_Move_Basic::best_xxx = 1.0;
    Ions_Move_Basic::relax_bfgs_init = 1.0;

    // call function
    im_cg.move0[0] = 1.0;
    std::ofstream ofs1("TestStartNoTrial_temp1.log");
    im_cg.start(ucell, force, etot, istep, update_iter, ofs1, etot_info, relax_method);
    ofs1.close();
    std::remove("TestStartNoTrial_temp1.log");
    int istep_2 = 2;
    im_cg.move0[0] = 1.0;
    force(0, 0) = 0.001;
    std::ofstream ofs2("TestStartNoTrial_temp2.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs2, etot_info, relax_method);
    ofs2.close();
    std::remove("TestStartNoTrial_temp2.log");
    std::ofstream ofs("TestStartNoTrial.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs, etot_info, relax_method);
    ofs.close();

    // Check output
    std::string expected_output = "\n Largest force is 0.0257111 eV/Angstrom while threshold is -1 eV/Angstrom\n\n"
                                  " Ion relaxation is not converged yet (threshold is 0.0257111)\n";
    std::ifstream ifs("TestStartNoTrial.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("TestStartNoTrial.log");

    EXPECT_THAT(output, testing::HasSubstr(expected_output));
    EXPECT_EQ(update_iter, 5);
    EXPECT_EQ(relax_method[0], "bfgs");
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::largest_grad, 0.001);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::best_xxx, 1.0);
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

    im_cg.setup_cg_grad(Ions_Move_Basic::dim, grad, grad0, cggrad, cggrad0, ncggrad, flag);

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

    im_cg.setup_cg_grad(Ions_Move_Basic::dim, grad, grad0, cggrad, cggrad0, ncggrad, flag);

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

    im_cg.setup_cg_grad(Ions_Move_Basic::dim, grad, grad0, cggrad, cggrad0, ncggrad, flag);

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
    double f_value = 0.0;

    im_cg.f_cal(Ions_Move_Basic::dim, g0, g1, f_value);

    EXPECT_DOUBLE_EQ(f_value, 3.0);
}

// Test function setup_move()
TEST_F(IonsMoveCGTest, SetupMove)
{
    Ions_Move_Basic::dim = 9;
    double trust_radius = 1.0;
    double cg_gradn[9] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double move[9] = {0.0};

    im_cg.setup_move(Ions_Move_Basic::dim, move, cg_gradn, trust_radius);

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
