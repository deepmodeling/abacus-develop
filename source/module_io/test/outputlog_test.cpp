#include "gtest/gtest.h"
#include <sstream>
#include "module_io/output_log.h"
#include "module_base/global_variable.h"
#include "module_base/constants.h"
#include <iostream>
#include <cstdio>
#include <unistd.h>
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