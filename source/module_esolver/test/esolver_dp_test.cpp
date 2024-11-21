#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "for_test.h"

#define private public
#define protected public
#include "../esolver_dp.h"
#undef private
/************************************************
 *  unit tests of class ESolver_DP
 ***********************************************/

/**
 * - Tested Functions:
 *   - ESolver_DP::init()
 *   - ESolver_DP::run()
 *   - ESolver_DP::cal_energy()
 *   - ESolver_DP::cal_force()
 *   - ESolver_DP::cal_stress()
 *   - ESolver_DP::post_process()
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
        ucell.iat2it = new int[2];
        ucell.iat2it[0] = 0;
        ucell.iat2it[1] = 1;
        ucell.iat2ia = new int[2];
        ucell.iat2ia[0] = 0;
        ucell.iat2ia[1] = 1;
        ucell.nat = 2;
        ucell.ntype = 2;
        ucell.atoms = new Atom[2];
        ucell.atoms[0].na = 1;
        ucell.atoms[1].na = 1;
        ucell.atoms[0].taud.resize(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
        ucell.atoms[1].taud.resize(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));

        ucell.atom_label = new std::string[2];
        ucell.atom_label[0] = "Cu";
        ucell.atom_label[1] = "Al";
        esolver->before_all_runners(inp, ucell);
    }

    void TearDown() override
    {
        // Clean up after each test
        delete esolver;
        delete[] ucell.iat2it;
        delete[] ucell.iat2ia;
        delete[] ucell.atoms;
        delete[] ucell.atom_label;
    }

    ModuleESolver::ESolver_DP* esolver;
    Input_para inp;
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
    EXPECT_EQ(esolver->atype[0], 0);
    EXPECT_EQ(esolver->atype[1], 0);
}

// Test the Run() funciton WARNING_QUIT
TEST_F(ESolverDPTest, RunWarningQuit)
{
    int istep = 0;

    testing::internal::CaptureStdout();
    EXPECT_EXIT(esolver->runner(istep, ucell), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Please recompile with -D__DPMD"));
}

// Test the cal_energy() funciton
TEST_F(ESolverDPTest, CalEnergy)
{
    double etot = 0.0;
    esolver->dp_potential = 9.8;
    etot = esolver->cal_energy();

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

    esolver->cal_force(force);

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

    esolver->cal_stress(stress);

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
    esolver->after_all_runners();
    GlobalV::ofs_running.close();

    std::string expected_output = "\n\n --------------------------------------------\n !FINAL_ETOT_IS 133.3358404 eV\n "
                                  "--------------------------------------------\n\n\n";
    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_EQ(expected_output, output);
}
