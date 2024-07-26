#include <unistd.h>

#include <cstdio>
#include <iostream>
#include <sstream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_base/constants.h"
#include "module_base/global_variable.h"
#include "module_io/output_log.h"

#ifdef __MPI
#include "module_basis/module_pw/test/test_tool.h"
#include "mpi.h"
#endif
/**
 * - Tested Functions:
 *  - output_convergence_after_scf()
 *  - output_efermi()
*/

// Test the output_convergence_after_scf function
TEST(OutputConvergenceAfterSCFTest, TestConvergence) {
    bool convergence = true;
    double energy = 2.0;
    std::ofstream ofs_running("test_output_convergence.txt");
    ModuleIO::output_convergence_after_scf(convergence, energy, ofs_running);
    ofs_running.close();

    std::ifstream ifs_running("test_output_convergence.txt");
    std::stringstream ss;
    ss << ifs_running.rdbuf();
    std::string file_content = ss.str();
    ifs_running.close();

    std::string expected_content = "\n charge density convergence is achieved\n"
                                   " final etot is 27.211396 eV\n";

    EXPECT_EQ(file_content, expected_content);
     std::remove("test_output_convergence.txt");
}

TEST(OutputConvergenceAfterSCFTest, TestNotConvergence) {
    bool convergence = false;
    double energy = 2.0;
    std::ofstream ofs_running("test_output_convergence_noconvergence.txt");
    testing::internal::CaptureStdout();
    ModuleIO::output_convergence_after_scf(convergence, energy, ofs_running);
    std::string screen_output = testing::internal::GetCapturedStdout();
    ofs_running.close();

    std::ifstream ifs_running("test_output_convergence_noconvergence.txt");
    std::stringstream ss;
    ss << ifs_running.rdbuf();
    std::string file_content = ss.str();
    ifs_running.close();

    std::string expected_content = " !! convergence has not been achieved @_@\n";
    std::string expected_content_screen = " !! CONVERGENCE HAS NOT BEEN ACHIEVED !!\n";

    EXPECT_EQ(file_content, expected_content);
    EXPECT_EQ(screen_output, expected_content_screen);
    std::remove("test_output_convergence_noconvergence.txt");
}

// Test the output_efermi function
TEST(OutputEfermiTest, TestConvergence) {
    bool convergence = true;
    double efermi = 1.0;
    std::ofstream ofs_running("test_output_efermi.txt");
    ModuleIO::output_efermi(convergence, efermi, ofs_running);
    ofs_running.close();

    std::ifstream ifs_running("test_output_efermi.txt");
    std::stringstream ss;
    ss << ifs_running.rdbuf();
    std::string file_content = ss.str();
    ifs_running.close();

    std::string expected_content = " EFERMI = 13.605698 eV\n";

    EXPECT_EQ(file_content, expected_content);
    std::remove("test_output_efermi.txt");
}

TEST(OutputEfermiTest, TestNotConvergence) {
    bool convergence = false;
    double efermi = 1.0;
    std::ofstream ofs_running("test_output_efermi_noconvergence.txt");
    ModuleIO::output_efermi(convergence, efermi, ofs_running);
    ofs_running.close();

    std::ifstream ifs_running("test_output_efermi_noconvergence.txt");
    std::stringstream ss;
    ss << ifs_running.rdbuf();
    std::string file_content = ss.str();
    ifs_running.close();

    std::string expected_content = ""; // No output expected if convergence is false

    EXPECT_EQ(file_content, expected_content);
    std::remove("test_output_efermi_noconvergence.txt");
}

TEST(OutputEfermiTest, TestMOutputLevel) {
    bool convergence = true;
    double efermi = 1.0;
    GlobalV::OUT_LEVEL = "m"; // Setting output level to "m"
    std::ofstream ofs_running("test_output_efermi_m_outputlevel.txt");
    ModuleIO::output_efermi(convergence, efermi, ofs_running);
    ofs_running.close();

    std::ifstream ifs_running("test_output_efermi_m_outputlevel.txt");
    std::stringstream ss;
    ss << ifs_running.rdbuf();
    std::string file_content = ss.str();
    ifs_running.close();

    std::string expected_content = ""; // No output expected if OUT_LEVEL is "m"

    EXPECT_EQ(file_content, expected_content);
    std::remove("test_output_efermi_m_outputlevel.txt");
}

UnitCell::UnitCell()
{
    ntype = 1;
    nat = 2;
    atoms = new Atom[ntype];

    latvec.e11 = latvec.e22 = latvec.e33 = 10;
    latvec.e12 = latvec.e13 = latvec.e23 = 0;
    latvec.e21 = latvec.e31 = latvec.e32 = 0;

    atoms[0].taud = new ModuleBase::Vector3<double>[2];
    atoms[0].taud[0].set(0.5456, 0, 0.54275);
    atoms[0].taud[1].set(0.54, 0.8495, 0.34175);
}
UnitCell::~UnitCell()
{
    if (atoms != nullptr)
        delete[] atoms;
}
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
Atom::Atom()
{
    na = 2;
    label = "Al";
}
Atom::~Atom()
{
    if (taud != nullptr)
        delete[] taud;
}
Atom_pseudo::Atom_pseudo()
{
}
Atom_pseudo::~Atom_pseudo()
{
}
pseudo::pseudo()
{
}
pseudo::~pseudo()
{
}

class OutputVacuumLevelTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis* rhopw = nullptr;

    virtual void SetUp()
    {
        rhopw = new ModulePW::PW_Basis;
    }
    virtual void TearDown()
    {
        if (rhopw != nullptr)
            delete rhopw;
    }
};

TEST_F(OutputVacuumLevelTest, OutputVacuumLevel)
{
    GlobalV::NSPIN = 1;
    std::cout << GlobalV::NSPIN << std::endl;
    UnitCell ucell;
    std::cout << ucell.latvec.e11 << std::endl;
    std::cout << ucell.latvec.e22 << std::endl;
    std::cout << ucell.latvec.e33 << std::endl;

#ifdef __MPI
    rhopw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, MPI_COMM_WORLD);
#endif
    rhopw->initgrids(1.5, ucell.latvec, 60);
    rhopw->initparameters(false, 60);
    rhopw->setuptransform();
    rhopw->collect_local_pw();

    double** rho = new double*[1];
    rho[0] = new double[rhopw->nrxx];
    double* v_elecstat = new double[rhopw->nrxx];
    for (int ir = 0; ir < rhopw->nrxx; ++ir)
    {
        rho[0][ir] = 0.01 * ir;
        v_elecstat[ir] = 0.02 * ir;
    }

    std::ofstream ofs_running("test_output_vacuum_level.txt");
    ModuleIO::output_vacuum_level(rhopw, &ucell, rho, v_elecstat, ofs_running);
    ofs_running.close();

    std::ifstream ifs_running("test_output_vacuum_level.txt");
    std::stringstream ss;
    ss << ifs_running.rdbuf();
    std::string file_content = ss.str();
    ifs_running.close();

    std::string expected_content = "The vacuum level is 1523.7 eV\n";

    EXPECT_EQ(file_content, expected_content);
    std::remove("test_output_vacuum_level.txt");

    delete[] rho[0];
    delete[] rho;
    delete[] v_elecstat;
}

TEST(PrintForce, PrintForce)
{
    UnitCell ucell;
    GlobalV::TEST_FORCE = 1;
    std::string name = "test";
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 1.0;
    force(0, 1) = 2.0;
    force(0, 2) = 3.0;
    force(1, 0) = 0.0;
    force(1, 1) = 0.0;
    force(1, 2) = 0.0;

    GlobalV::ofs_running.open("test.txt");
    ModuleIO::print_force(GlobalV::ofs_running, ucell, name, force, false);
    GlobalV::ofs_running.close();

    std::ifstream ifs("test.txt");
    std::string output_str;
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("---------------------------------------------------------------------------"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("test"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("---------------------------------------------------------------------------"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("Al1        25.7110532015        51.4221064030        77.1331596044"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("Al2         0.0000000000         0.0000000000         0.0000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("---------------------------------------------------------------------------"));
    ifs.close();
    std::remove("test.txt");
}

TEST(PrintStress, PrintStress)
{
    ModuleBase::matrix stress(3, 3);
    stress(0, 0) = 1.0;
    stress(0, 1) = 2.0;
    stress(0, 2) = 3.0;
    stress(1, 0) = 0.0;
    stress(1, 1) = 0.0;
    stress(1, 2) = 0.0;
    stress(2, 0) = 0.0;
    stress(2, 1) = 0.0;
    stress(2, 2) = 0.0;

    GlobalV::ofs_running.open("test.txt");
    ModuleIO::print_stress("TOTAL-STRESS", stress, true, false);
    GlobalV::ofs_running.close();

    std::ifstream ifs("test.txt");
    std::string output_str;
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("----------------------------------------------------------------"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("TOTAL-STRESS (KBAR)                                           "));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("----------------------------------------------------------------"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("   147105.2279754489    294210.4559508978    441315.6839263467"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("        0.0000000000         0.0000000000         0.0000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("        0.0000000000         0.0000000000         0.0000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("----------------------------------------------------------------"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr(" TOTAL-PRESSURE: 49035.075992 KBAR"));
    ifs.close();
    std::remove("test.txt");
}

int main(int argc, char** argv)
{
#ifdef __MPI
    setupmpi(argc, argv, GlobalV::NPROC, GlobalV::MY_RANK);
    divide_pools(GlobalV::NPROC,
                 GlobalV::MY_RANK,
                 GlobalV::NPROC_IN_POOL,
                 GlobalV::KPAR,
                 GlobalV::MY_POOL,
                 GlobalV::RANK_IN_POOL);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    finishmpi();
#endif
    return result;
}