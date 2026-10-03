#include "source_relax/relax_criteria.h"
#include "for_test.h"
#include "gtest/gtest.h"
#include "mock_remake_cell.h"
#include "source_relax/lattice_change_basic.h"
#include "source_relax/lattice_change_cg.h"

/************************************************
 *  unit tests of class Lattice_Change_CG
 ***********************************************/

class LatticeChangeCGTest : public ::testing::Test
{
  public:
    Relax_Criteria criteria;

  protected:
    void SetUp() override
    {
        // Initialize variables before each test
        Lattice_Change_Basic::dim = 9;
        Lattice_Change_Basic::stress_step = 1;
        Lattice_Change_Basic::update_iter = 5;
        lc_cg.allocate();
        etot_info.resize(2, 0.0);
    }

    void TearDown() override
    {
        // Clean up after each test
    }

    Lattice_Change_CG lc_cg;
    std::vector<double> etot_info;

    // Friendship is not inherited: a TEST_F body lives in a class derived from
    // this fixture, so the internals of Lattice_Change_CG (which friends the
    // fixture itself) are reached through these forwarders.
    const std::vector<double>& get_lat0() const
    {
        return lc_cg.lat0;
    }
    const std::vector<double>& get_grad0() const
    {
        return lc_cg.grad0;
    }
    const std::vector<double>& get_cg_grad0() const
    {
        return lc_cg.cg_grad0;
    }
    const std::vector<double>& get_move0() const
    {
        return lc_cg.move0;
    }
    void set_move0(const int i, const double value)
    {
        lc_cg.move0[i] = value;
    }
};

// Test whether the allocate() function can correctly allocate memory space
TEST_F(LatticeChangeCGTest, TestAllocate)
{
    Lattice_Change_Basic::dim = 4;
    lc_cg.allocate();

    // Check if allocated vectors are not empty
    EXPECT_EQ(get_lat0().size(), 4U);
    EXPECT_EQ(get_grad0().size(), 4U);
    EXPECT_EQ(get_cg_grad0().size(), 4U);
    EXPECT_EQ(get_move0().size(), 4U);
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
    EXPECT_DOUBLE_EQ(0.0, get_lat0()[0]);
    EXPECT_DOUBLE_EQ(0.0, get_grad0()[1]);
    EXPECT_DOUBLE_EQ(0.0, get_cg_grad0()[2]);
    EXPECT_DOUBLE_EQ(0.0, get_move0()[0]);
}

// Test function start() when converged
TEST_F(LatticeChangeCGTest, TestStartConverged)
{
    // setup data
    UnitCell ucell;
    ucell.lat_axis_free[0] = 1;
    ucell.lat_axis_free[1] = 1;
    ucell.lat_axis_free[2] = 1;
    ModuleBase::matrix stress(3, 3);
    double etot = 0.0;

    // call function
    std::ofstream ofs("test_lc_cg_start_converged.log");
    lc_cg.start(ucell, stress, etot, ofs, etot_info, criteria);
    ofs.close();

    // Check output
    std::string expected_output
        = " Largest stress is 0, movement is impossible.\n end of lattice optimization\n                              stress_step = 1\n       "
          "                  update iteration = 5\n";
    std::ifstream ifs("test_lc_cg_start_converged.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    ifs.close();
    std::remove("test_lc_cg_start_converged.log");
}

// Test function start() sd branch
TEST_F(LatticeChangeCGTest, TestStartSd)
{
    // setup data
    UnitCell ucell;
    ucell.lat_axis_free[0] = 1;
    ucell.lat_axis_free[1] = 1;
    ucell.lat_axis_free[2] = 1;
    ModuleBase::matrix stress(3, 3);
    stress(0, 0) = 0.01;
    double etot = 0.0;

    // call function
    std::ofstream ofs("test_lc_cg_start_sd.log");
    lc_cg.start(ucell, stress, etot, ofs, etot_info, criteria);
    ofs.close();

    // Check output: reporting moved to IonCellOptimizer::relax_step;
    // check_converged prints nothing here anymore.
    std::string expected_output = "";
    std::ifstream ifs("test_lc_cg_start_sd.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    EXPECT_DOUBLE_EQ(Lattice_Change_Basic::lattice_change_ini, 0.01);
    ifs.close();
    std::remove("test_lc_cg_start_sd.log");
}

// Test function start() trial branch with goto
TEST_F(LatticeChangeCGTest, TestStartTrialGoto)
{
    // setup data
    UnitCell ucell;
    ucell.lat_axis_free[0] = 1;
    ucell.lat_axis_free[1] = 1;
    ucell.lat_axis_free[2] = 1;
    ModuleBase::matrix stress(3, 3);
    stress(0, 1) = 0.01;
    double etot = 0.0;

    // call function
    set_move0(0, 1.0);
    std::ofstream ofs1("test_lc_cg_start_trial_goto_temp1.log");
    lc_cg.start(ucell, stress, etot, ofs1, etot_info, criteria);
    ofs1.close();
    std::remove("test_lc_cg_start_trial_goto_temp1.log");
    Lattice_Change_Basic::stress_step = 2;
    set_move0(0, 10.0);
    std::ofstream ofs("test_lc_cg_start_trial_goto.log");
    lc_cg.start(ucell, stress, etot, ofs, etot_info, criteria);
    ofs.close();

    // Check output: reporting moved to IonCellOptimizer::relax_step;
    // check_converged prints nothing here anymore.
    std::string expected_output = "";
    std::ifstream ifs("test_lc_cg_start_trial_goto.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    EXPECT_NEAR(Lattice_Change_Basic::lattice_change_ini, 10.000004999998749, 1e-12);
    ifs.close();
    std::remove("test_lc_cg_start_trial_goto.log");
}

// Test function start() trial branch without goto
TEST_F(LatticeChangeCGTest, TestStartTrial)
{
    // setup data
    UnitCell ucell;
    ucell.lat_axis_free[0] = 1;
    ucell.lat_axis_free[1] = 1;
    ucell.lat_axis_free[2] = 1;
    ModuleBase::matrix stress(3, 3);
    stress(0, 1) = 0.01;
    double etot = 0.0;

    // call function
    std::ofstream ofs1("test_lc_cg_start_trial_temp1.log");
    lc_cg.start(ucell, stress, etot, ofs1, etot_info, criteria);
    ofs1.close();
    std::remove("test_lc_cg_start_trial_temp1.log");
    Lattice_Change_Basic::stress_step = 2;
    std::ofstream ofs("test_lc_cg_start_trial.log");
    lc_cg.start(ucell, stress, etot, ofs, etot_info, criteria);
    ofs.close();

    // Check output: reporting moved to IonCellOptimizer::relax_step;
    // check_converged prints nothing here anymore.
    std::string expected_output = "";
    std::ifstream ifs("test_lc_cg_start_trial.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    EXPECT_NEAR(Lattice_Change_Basic::lattice_change_ini, 70.000034999991243, 1e-12);
    ifs.close();
    std::remove("test_lc_cg_start_trial.log");
}

// Test function start() no trial branch with goto case 1
TEST_F(LatticeChangeCGTest, TestStartNoTrialGotoCase1)
{
    // setup data
    UnitCell ucell;
    ucell.lat_axis_free[0] = 1;
    ucell.lat_axis_free[1] = 1;
    ucell.lat_axis_free[2] = 1;
    ModuleBase::matrix stress(3, 3);
    stress(0, 1) = 0.01;
    double etot = 0.0;

    // call function
    std::ofstream ofs1("test_lc_cg_start_notrial_goto_case1_temp1.log");
    lc_cg.start(ucell, stress, etot, ofs1, etot_info, criteria);
    ofs1.close();
    std::remove("test_lc_cg_start_notrial_goto_case1_temp1.log");
    Lattice_Change_Basic::stress_step = 2;
    std::ofstream ofs2("test_lc_cg_start_notrial_goto_case1_temp2.log");
    lc_cg.start(ucell, stress, etot, ofs2, etot_info, criteria);
    ofs2.close();
    std::remove("test_lc_cg_start_notrial_goto_case1_temp2.log");
    std::ofstream ofs("test_lc_cg_start_notrial_goto_case1.log");
    lc_cg.start(ucell, stress, etot, ofs, etot_info, criteria);
    ofs.close();

    // Check output: reporting moved to IonCellOptimizer::relax_step;
    // check_converged prints nothing here anymore.
    std::string expected_output = "";
    std::ifstream ifs("test_lc_cg_start_notrial_goto_case1.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    EXPECT_NEAR(Lattice_Change_Basic::lattice_change_ini, 490.00024499993867, 1e-12);
    ifs.close();
    std::remove("test_lc_cg_start_notrial_goto_case1.log");
}

// Test function start() no trial branch with goto case 2
TEST_F(LatticeChangeCGTest, TestStartNoTrialGotoCase2)
{
    // setup data
    UnitCell ucell;
    ucell.lat_axis_free[0] = 1;
    ucell.lat_axis_free[1] = 1;
    ucell.lat_axis_free[2] = 1;
    ModuleBase::matrix stress(3, 3);
    stress(0, 1) = 0.01;
    double etot = 0.0;

    // call function
    set_move0(0, 0.1);
    std::ofstream ofs1("test_lc_cg_start_notrial_goto_case2_temp1.log");
    lc_cg.start(ucell, stress, etot, ofs1, etot_info, criteria);
    ofs1.close();
    std::remove("test_lc_cg_start_notrial_goto_case2_temp1.log");
    Lattice_Change_Basic::stress_step = 2;
    std::ofstream ofs2("test_lc_cg_start_notrial_goto_case2_temp2.log");
    lc_cg.start(ucell, stress, etot, ofs2, etot_info, criteria);
    ofs2.close();
    std::remove("test_lc_cg_start_notrial_goto_case2_temp2.log");
    std::ofstream ofs("test_lc_cg_start_notrial_goto_case2.log");
    set_move0(0, 0.1);
    stress(0, 1) = 0.0001;
    lc_cg.start(ucell, stress, etot, ofs, etot_info, criteria);
    ofs.close();

    // Check output: reporting moved to IonCellOptimizer::relax_step;
    // check_converged prints nothing here anymore.
    std::string expected_output = "";
    std::ifstream ifs("test_lc_cg_start_notrial_goto_case2.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    EXPECT_NEAR(Lattice_Change_Basic::lattice_change_ini, 3430.0017149995706, 1e-12);
    ifs.close();
    std::remove("test_lc_cg_start_notrial_goto_case2.log");
}

// Test function start() no trial branch without goto
TEST_F(LatticeChangeCGTest, TestStartNoTrial)
{
    // setup data
    UnitCell ucell;
    ucell.lat_axis_free[0] = 1;
    ucell.lat_axis_free[1] = 1;
    ucell.lat_axis_free[2] = 1;
    ModuleBase::matrix stress(3, 3);
    stress(0, 1) = 0.01;
    double etot = 0.0;

    // call function
    set_move0(0, 1.0);
    std::ofstream ofs1("test_lc_cg_start_notrial_temp1.log");
    lc_cg.start(ucell, stress, etot, ofs1, etot_info, criteria);
    ofs1.close();
    std::remove("test_lc_cg_start_notrial_temp1.log");
    Lattice_Change_Basic::stress_step = 2;
    set_move0(0, 10.0);
    std::ofstream ofs2("test_lc_cg_start_notrial_temp2.log");
    lc_cg.start(ucell, stress, etot, ofs2, etot_info, criteria);
    ofs2.close();
    std::remove("test_lc_cg_start_notrial_temp2.log");
    std::ofstream ofs("test_lc_cg_start_notrial.log");
    lc_cg.start(ucell, stress, etot, ofs, etot_info, criteria);
    ofs.close();

    // Check output: reporting moved to IonCellOptimizer::relax_step;
    // check_converged prints nothing here anymore.
    std::string expected_output = "";
    std::ifstream ifs("test_lc_cg_start_notrial.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected_output, output);
    EXPECT_NEAR(Lattice_Change_Basic::lattice_change_ini, 96040.106328872833, 1e-12);
    ifs.close();
    std::remove("test_lc_cg_start_notrial.log");
}

