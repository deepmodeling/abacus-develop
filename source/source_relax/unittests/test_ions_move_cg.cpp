#include "source_relax/relax_criteria.h"
#include <regex>
#include "for_test.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "source_relax/ions_move_basic.h"
#include "source_relax/ions_move_cg.h"

/************************************************
 *  unit tests of class Ions_Move_CG
 ***********************************************/

class IonsMoveCGTest : public ::testing::Test
{
  public:
    Relax_Criteria criteria;

  protected:
    void SetUp() override
    {
        // Initialize variables before each test
        Ions_Move_Basic::dim = 6;
        update_iter = 5;
        im_cg.allocate(Ions_Move_Basic::dim);
        criteria.force_thr = 0.001;

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

    // Friendship is not inherited: a TEST_F body lives in a class derived from
    // this fixture, so the internals of Ions_Move_CG (which friends the fixture
    // itself) are reached through these forwarders.
    const std::vector<double>& get_pos0() const
    {
        return im_cg.pos0;
    }
    const std::vector<double>& get_grad0() const
    {
        return im_cg.grad0;
    }
    const std::vector<double>& get_cg_grad0() const
    {
        return im_cg.cg_grad0;
    }
    const std::vector<double>& get_move0() const
    {
        return im_cg.move0;
    }
    void set_move0(const int i, const double value)
    {
        im_cg.move0[i] = value;
    }
};

// Test whether the allocate() function can correctly allocate memory space
TEST_F(IonsMoveCGTest, TestAllocate)
{
    const int dim = 4;
    im_cg.allocate(dim);

    // Check if allocated vectors are not empty
    EXPECT_EQ(get_pos0().size(), 4U);
    EXPECT_EQ(get_grad0().size(), 4U);
    EXPECT_EQ(get_cg_grad0().size(), 4U);
    EXPECT_EQ(get_move0().size(), 4U);
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
    EXPECT_DOUBLE_EQ(0.0, get_pos0()[0]);
    EXPECT_DOUBLE_EQ(0.0, get_grad0()[1]);
    EXPECT_DOUBLE_EQ(0.0, get_cg_grad0()[2]);
    EXPECT_DOUBLE_EQ(0.0, get_move0()[0]);
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
    im_cg.start(ucell, force, etot, istep, update_iter, ofs, etot_info, relax_method, criteria);
    ofs.close();

    // Check output
    // Reporting moved to IonCellOptimizer::relax_step; check_converged no longer
    // prints "Largest force is ...". The zero-force branch message and terminate
    // output remain.
    std::string expected_output = " largest force is 0, no movement is possible.\n it may converged, otherwise no "
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
    im_cg.start(ucell, force, etot, istep, update_iter, ofs, etot_info, relax_method, criteria);
    ofs.close();

    // Check output: reporting moved to IonCellOptimizer::relax_step;
    // check_converged prints nothing here anymore.
    std::string expected_output = "";
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
    set_move0(0, 1.0);
    std::ofstream ofs1("TestStartTrialGoto_temp1.log");
    im_cg.start(ucell, force, etot, istep, update_iter, ofs1, etot_info, relax_method, criteria);
    ofs1.close();
    std::remove("TestStartTrialGoto_temp1.log");
    int istep_2 = 2;
    set_move0(0, 10.0);
    force(0, 0) = 0.001;
    relax_method = {"cg_bfgs", "1"};
    std::ofstream ofs("TestStartTrialGoto.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs, etot_info, relax_method, criteria);
    ofs.close();

    // Check output: reporting moved to IonCellOptimizer::relax_step;
    // check_converged prints nothing here anymore.
    std::string expected_output = "";
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
    set_move0(0, 1.0);
    std::ofstream ofs1("TestStartTrial_temp1.log");
    im_cg.start(ucell, force, etot, istep, update_iter, ofs1, etot_info, relax_method, criteria);
    ofs1.close();
    std::remove("TestStartTrial_temp1.log");
    int istep_2 = 2;
    set_move0(0, 10.0);
    std::ofstream ofs("TestStartTrial.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs, etot_info, relax_method, criteria);
    ofs.close();

    // Check output: reporting moved to IonCellOptimizer::relax_step;
    // check_converged prints nothing here anymore.
    std::string expected_output = "";
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
    set_move0(0, 1.0);
    std::ofstream ofs1("TestStartNoTrialGotoCase1_temp1.log");
    im_cg.start(ucell, force, etot, istep, update_iter, ofs1, etot_info, relax_method, criteria);
    ofs1.close();
    std::remove("TestStartNoTrialGotoCase1_temp1.log");
    int istep_2 = 2;
    std::ofstream ofs2("TestStartNoTrialGotoCase1_temp2.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs2, etot_info, relax_method, criteria);
    ofs2.close();
    std::remove("TestStartNoTrialGotoCase1_temp2.log");
    set_move0(0, 1.0);
    force(0, 0) = 0.001;
    relax_method = {"cg_bfgs", "1"};
    std::ofstream ofs("TestStartNoTrialGotoCase1.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs, etot_info, relax_method, criteria);
    ofs.close();

    // Check output: reporting moved to IonCellOptimizer::relax_step;
    // check_converged prints nothing here anymore.
    std::string expected_output = "";
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
    set_move0(0, 1.0);
    std::ofstream ofs1("TestStartNoTrialGotoCase2_temp1.log");
    im_cg.start(ucell, force, etot, istep, update_iter, ofs1, etot_info, relax_method, criteria);
    ofs1.close();
    std::remove("TestStartNoTrialGotoCase2_temp1.log");
    int istep_2 = 2;
    set_move0(0, 10.0);
    std::ofstream ofs2("TestStartNoTrialGotoCase2_temp2.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs2, etot_info, relax_method, criteria);
    ofs2.close();
    std::remove("TestStartNoTrialGotoCase2_temp2.log");
    relax_method = {"cg_bfgs", "1"};
    std::ofstream ofs("TestStartNoTrialGotoCase2.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs, etot_info, relax_method, criteria);
    ofs.close();

    // Check output: reporting moved to IonCellOptimizer::relax_step;
    // check_converged prints nothing here anymore.
    std::string expected_output = "";
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
    set_move0(0, 1.0);
    std::ofstream ofs1("TestStartNoTrial_temp1.log");
    im_cg.start(ucell, force, etot, istep, update_iter, ofs1, etot_info, relax_method, criteria);
    ofs1.close();
    std::remove("TestStartNoTrial_temp1.log");
    int istep_2 = 2;
    set_move0(0, 1.0);
    force(0, 0) = 0.001;
    std::ofstream ofs2("TestStartNoTrial_temp2.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs2, etot_info, relax_method, criteria);
    ofs2.close();
    std::remove("TestStartNoTrial_temp2.log");
    std::ofstream ofs("TestStartNoTrial.log");
    im_cg.start(ucell, force, etot, istep_2, update_iter, ofs, etot_info, relax_method, criteria);
    ofs.close();

    // Check output: reporting moved to IonCellOptimizer::relax_step;
    // check_converged prints nothing here anymore.
    std::string expected_output = "";
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
