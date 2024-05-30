#include "gtest/gtest.h"
#include "../write_wfc_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"

//Can not test the 2d-to-grid convertion function now, because the refactor of related functions is not finished.
//So just mock the function globalIndex() here for compiling.
int Local_Orbital_wfc::globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    return 0;
}


TEST(GenWfcLcaoFnameTest, OutType1GammaOnlyOutAppFlagTrue)
{
    int out_type = 1;
    bool gamma_only = true;
    bool out_app_flag = true;
    int ik = 0;
    int istep = 0;
    GlobalV::global_out_dir = "/path/to/global_out_dir/";

    std::string expected_output = "/path/to/global_out_dir/WFC_LCAO_GAMMA1.txt";
    std::string result = ModuleIO::wfc_lcao_gen_fname(out_type, gamma_only, out_app_flag, ik, istep);

    EXPECT_EQ(result, expected_output);
}

TEST(GenWfcLcaoFnameTest, OutType2GammaOnlyOutAppFlagFalse)
{
    int out_type = 2;
    bool gamma_only = true;
    bool out_app_flag = false;
    int ik = 1;
    int istep = 2;
    GlobalV::global_out_dir = "/path/to/global_out_dir/";

    std::string expected_output = "/path/to/global_out_dir/WFC_LCAO_GAMMA2_ION3.dat";
    std::string result = ModuleIO::wfc_lcao_gen_fname(out_type, gamma_only, out_app_flag, ik, istep);

    EXPECT_EQ(result, expected_output);
}

TEST(GenWfcLcaoFnameTest, OutTypeInvalid)
{
    int out_type = 3;
    bool gamma_only = false;
    bool out_app_flag = true;
    int ik = 2;
    int istep = 3;
    GlobalV::global_out_dir = "/path/to/global_out_dir/";

    std::string expected_output = "/path/to/global_out_dir/WFC_LCAO_K3.txt";
    // catch the screen output
    testing::internal::CaptureStdout();
    std::string result = ModuleIO::wfc_lcao_gen_fname(out_type, gamma_only, out_app_flag, ik, istep);
    std::string output = testing::internal::GetCapturedStdout();


    EXPECT_EQ(result, expected_output);
}
