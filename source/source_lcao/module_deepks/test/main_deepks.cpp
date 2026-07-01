#include "deepks_test_runner.h"

#include <gtest/gtest.h>

#ifdef __MPI
#include <mpi.h>
#endif

#include <cerrno>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>

#ifndef DEEPKS_UT_CHECK_NAME
#error "DEEPKS_UT_CHECK_NAME must be defined by CMake."
#endif

#ifndef DEEPKS_UT_CASE_DIR
#error "DEEPKS_UT_CASE_DIR must be defined by CMake."
#endif

#ifndef DEEPKS_UT_RUNNER
#error "DEEPKS_UT_RUNNER must be defined by CMake."
#endif

template <typename T>
void DEEPKS_UT_RUNNER(test_deepks<T>& test);

namespace
{
std::string shell_quote(const std::string& value)
{
    std::string quoted = "'";
    for (std::string::const_iterator it = value.begin(); it != value.end(); ++it)
    {
        if (*it == '\'')
        {
            quoted += "'\\''";
        }
        else
        {
            quoted += *it;
        }
    }
    quoted += "'";
    return quoted;
}

void prepare_workdir()
{
    const std::string case_dir = DEEPKS_UT_CASE_DIR;
    const std::string check_name = DEEPKS_UT_CHECK_NAME;
    const std::string run_root = "deepks_unit_run_" + case_dir + "_" + check_name;

    std::ostringstream command;
    command << "rm -rf " << shell_quote(run_root) << " && "
            << "mkdir -p " << shell_quote(run_root) << " && "
            << "cp -R support/PP_ORB " << shell_quote(run_root + "/PP_ORB") << " && "
            << "cp -R " << shell_quote("support/" + case_dir) << " " << shell_quote(run_root + "/" + case_dir);

    ASSERT_EQ(std::system(command.str().c_str()), 0) << "Failed to prepare DeePKS unit-test work directory";

    const std::string workdir = run_root + "/" + case_dir;
    ASSERT_EQ(chdir(workdir.c_str()), 0) << "Failed to chdir to " << workdir << ": " << std::strerror(errno);
}

template <typename T>
void run_typed_check()
{
    test_deepks<T> test;
    test.preparation();
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    DEEPKS_UT_RUNNER(test);
}

bool gamma_only_case()
{
    std::ifstream ifs("INPUT");
    std::string key;
    bool gamma_only_local = false;
    ifs >> key >> gamma_only_local;
    return gamma_only_local;
}
} // namespace

TEST(DeePKSUnitTest, ConfiguredCheck)
{
    prepare_workdir();

    if (gamma_only_case())
    {
        run_typed_check<double>();
    }
    else
    {
        run_typed_check<std::complex<double>>();
    }
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
