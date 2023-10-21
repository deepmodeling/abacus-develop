#include "for_test.h"
#define private public
#define protected public
#include "../esolver_dp.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
/************************************************
 *  unit tests of class ESolver_DP
 ***********************************************/

/**
 * - Tested Functions:
 *   - ESolver_DP::Init()
 *   - ESolver_DP::Run()
 *   - ESolver_DP::cal_Energy()
 *   - ESolver_DP::cal_Force()
 *   - ESolver_DP::cal_Stress()
 *   - ESolver_DP::postprocess()
 *   - ESolver_DP::type_map()
 */
namespace ModuleIO
{
void print_force(std::ofstream& ofs_running,
                 const UnitCell& cell,
                 const std::string& name,
                 const ModuleBase::matrix& force,
                 bool ry = true)
{
}
void print_stress(const std::string& name, const ModuleBase::matrix& scs, const bool screen, const bool ry)
{
}
} // namespace ModuleIO

class ESolverDPTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        // Initialize variables before each test
        esolver = new ModuleESolver::ESolver_DP("./support/case_1.pb");
        esolver->Init(inp, ucell);
    }

    void TearDown() override
    {
        // Clean up after each test
    }

    ModuleESolver::ESolver_DP* esolver;
    Input inp;
    UnitCell ucell;
};

// Test the Init() funciton case 1
TEST_F(ESolverDPTest, InitCase1)
{
    // Check the initialized variables
    EXPECT_DOUBLE_EQ(esolver->dp_potential, 0.0);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(esolver->dp_virial(i, j), 0.0);
            EXPECT_DOUBLE_EQ(esolver->cell[3 * i + j], 0.0);
        }
    }
    for (int i = 0; i < ucell.nat; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(esolver->dp_force(i, j), 0.0);
            EXPECT_DOUBLE_EQ(esolver->coord[3 * i + j], 0.0);
        }
    }
    EXPECT_EQ(esolver->dp_type[0], 1);
    EXPECT_EQ(esolver->dp_type[1], 0);
    EXPECT_EQ(esolver->atype[0], 1);
    EXPECT_EQ(esolver->atype[1], 0);
}

// Test the Init() funciton case 2
TEST_F(ESolverDPTest, InitCase2)
{
    esolver->dp_type[0] = 0;
    esolver->dp_type[1] = 0;
    esolver->dp_file = "./support/case_2.pb";
    esolver->Init(inp, ucell);

    // Check the initialized variables
    EXPECT_EQ(esolver->dp_type[0], 0);
    EXPECT_EQ(esolver->dp_type[1], 0);
    EXPECT_EQ(esolver->atype[0], 0);
    EXPECT_EQ(esolver->atype[1], 1);
}

// Test the Run() funciton WARNING_QUIT
TEST_F(ESolverDPTest, RunWarningQuit)
{
    int istep = 0;

    testing::internal::CaptureStdout();
    EXPECT_EXIT(esolver->Run(istep, ucell), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Please recompile with -D__DPMD"));
}

// Test the cal_Energy() funciton
TEST_F(ESolverDPTest, CalEnergy)
{
    double etot = 0.0;
    esolver->dp_potential = 9.8;
    etot = esolver->cal_Energy();

    // Check the results
    EXPECT_DOUBLE_EQ(etot, 9.8);
}

// Test the cal_Force() funciton
TEST_F(ESolverDPTest, CalForce)
{
    ModuleBase::matrix force(ucell.nat, 3);
    for (int i = 0; i < ucell.nat; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            esolver->dp_force(i, j) = 3.0 * i + j;
        }
    }

    esolver->cal_Force(force);

    // Check the results
    for (int i = 0; i < ucell.nat; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(force(i, j), 3.0 * i + j);
        }
    }
}

// Test the cal_Stress() funciton
TEST_F(ESolverDPTest, CalStress)
{
    ModuleBase::matrix stress(3, 3);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            esolver->dp_virial(i, j) = 3.0 * i + j;
        }
    }

    esolver->cal_Stress(stress);

    // Check the results
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(stress(i, j), 3.0 * i + j);
        }
    }
}

// Test the postprocess() funciton
TEST_F(ESolverDPTest, Postprocess)
{
    esolver->dp_potential = 9.8;

    // Check the results
    GlobalV::ofs_running.open("log");
    esolver->postprocess();
    GlobalV::ofs_running.close();

    std::string expected_output = "\n\n --------------------------------------------\n !FINAL_ETOT_IS 133.3358404 eV\n "
                                  "--------------------------------------------\n\n\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_EQ(expected_output, output);
}

// Test the type_map() funciton when find type_map
TEST_F(ESolverDPTest, TypeMapCase1)
{
    bool find = false;
    esolver->dp_file = "./support/case_1.pb";
    GlobalV::ofs_running.open("log");
    find = esolver->type_map(ucell);
    GlobalV::ofs_running.close();

    std::string expected_output = "\nDetermine the type map from DP model\nntype read from DP model: 2\n  Cu  Al\n\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    // Check the results
    EXPECT_TRUE(find);
    EXPECT_EQ(expected_output, output);
    EXPECT_EQ(esolver->dp_type[0], 1);
    EXPECT_EQ(esolver->dp_type[1], 0);
}

// Test the type_map() funciton when find not type_map
TEST_F(ESolverDPTest, TypeMapCase2)
{
    bool find = false;
    esolver->dp_type[0] = 0;
    esolver->dp_type[1] = 0;
    esolver->dp_file = "./support/case_2.pb";
    find = esolver->type_map(ucell);

    // Check the results
    EXPECT_FALSE(find);
    EXPECT_EQ(esolver->dp_type[0], 0);
    EXPECT_EQ(esolver->dp_type[1], 0);
}

// Test the type_map() funciton WARNING_QUIT
TEST_F(ESolverDPTest, TypeMapWarningQuit)
{
    esolver->dp_file = "./support/case_3.pb";

    testing::internal::CaptureStdout();
    EXPECT_EXIT(esolver->type_map(ucell), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("can not find the DP model"));
}