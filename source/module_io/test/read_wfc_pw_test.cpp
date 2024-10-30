#include "module_io/read_wfc_pw.h"

#define private public
#include "module_parameter/parameter.h"
#undef private
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#ifdef __MPI
#include "module_base/parallel_global.h"
#include "module_basis/module_pw/test/test_tool.h"
#include "mpi.h"
#endif

/**
 * - Tested Functions:
 *  - read_wfc_pw()
 */

class ReadWfcPwTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis_K* wfcpw = nullptr;
    ModuleBase::Vector3<double>* kvec_d = nullptr;
    int nkstot = 8;

    virtual void SetUp()
    {
        wfcpw = new ModulePW::PW_Basis_K;
        kvec_d = new ModuleBase::Vector3<double>[nkstot];
    }
    virtual void TearDown()
    {
        if (wfcpw != nullptr) {
            delete wfcpw;
}
        if (kvec_d != nullptr) {
            delete[] kvec_d;
}
    }
};

// Test the read_wfc_pw function
TEST_F(ReadWfcPwTest, ReadWfcPw)
{
    std::string filename = "./support/WAVEFUNC1.dat";

#ifdef __MPI
    wfcpw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
    wfcpw->initgrids(5.3233, ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 12, 12, 12);
    wfcpw->initparameters(false, 20, 8, kvec_d);
    wfcpw->setuptransform();
    wfcpw->collect_local_pw();

    PARAM.input.nbands = 8;
    const int nbasis = wfcpw->npwk[0];
    ModuleBase::ComplexMatrix wfcatom(PARAM.input.nbands, nbasis);
    ModuleIO::read_wfc_pw(filename, wfcpw, 0, nkstot, wfcatom);

    if (GlobalV::NPROC_IN_POOL == 1)
    {
        EXPECT_DOUBLE_EQ(wfcatom(0, 0).real(), -0.017953720885562179);
        EXPECT_DOUBLE_EQ(wfcatom(0, 0).imag(), 0.035959236666204548);
        EXPECT_DOUBLE_EQ(wfcatom(0, wfcpw->npwk[0] - 1).real(), -0.021041484787624309);
        EXPECT_DOUBLE_EQ(wfcatom(0, wfcpw->npwk[0] - 1).imag(), 0.042143574695220835);
        EXPECT_DOUBLE_EQ(wfcatom(7, wfcpw->npwk[0] - 1).real(), -0.011075023130363163);
        EXPECT_DOUBLE_EQ(wfcatom(7, wfcpw->npwk[0] - 1).imag(), 0.017342817352703006);
    }
    else if (GlobalV::NPROC_IN_POOL == 4)
    {
        if (GlobalV::RANK_IN_POOL == 0)
        {
            EXPECT_DOUBLE_EQ(wfcatom(0, 0).real(), -0.017953720885562179);
            EXPECT_DOUBLE_EQ(wfcatom(0, 0).imag(), 0.035959236666204548);
            EXPECT_DOUBLE_EQ(wfcatom(0, wfcpw->npwk[0] - 1).real(), -0.021041489893031052);
            EXPECT_DOUBLE_EQ(wfcatom(0, wfcpw->npwk[0] - 1).imag(), 0.04214358371778136);
            EXPECT_DOUBLE_EQ(wfcatom(7, wfcpw->npwk[0] - 1).real(), -0.0048838415065336213);
            EXPECT_DOUBLE_EQ(wfcatom(7, wfcpw->npwk[0] - 1).imag(), 0.0078610803827715778);
        }
        else if (GlobalV::RANK_IN_POOL == 1)
        {
            EXPECT_DOUBLE_EQ(wfcatom(0, 0).real(), -0.021041442735212794);
            EXPECT_DOUBLE_EQ(wfcatom(0, 0).imag(), 0.042143594413604955);
            EXPECT_DOUBLE_EQ(wfcatom(0, wfcpw->npwk[0] - 1).real(), -0.021041483700533322);
            EXPECT_DOUBLE_EQ(wfcatom(0, wfcpw->npwk[0] - 1).imag(), 0.042143578780482846);
            EXPECT_DOUBLE_EQ(wfcatom(7, wfcpw->npwk[0] - 1).real(), -0.0052306421970293327);
            EXPECT_DOUBLE_EQ(wfcatom(7, wfcpw->npwk[0] - 1).imag(), 0.008388410016516171);
        }
        else if (GlobalV::RANK_IN_POOL == 2)
        {
            EXPECT_DOUBLE_EQ(wfcatom(0, 0).real(), -0.021041446966127354);
            EXPECT_DOUBLE_EQ(wfcatom(0, 0).imag(), 0.042143576374759073);
            EXPECT_DOUBLE_EQ(wfcatom(0, wfcpw->npwk[0] - 1).real(), -0.021041484787624309);
            EXPECT_DOUBLE_EQ(wfcatom(0, wfcpw->npwk[0] - 1).imag(), 0.042143574695220835);
            EXPECT_DOUBLE_EQ(wfcatom(7, wfcpw->npwk[0] - 1).real(), -0.011075023130363163);
            EXPECT_DOUBLE_EQ(wfcatom(7, wfcpw->npwk[0] - 1).imag(), 0.017342817352703006);
        }
        else if (GlobalV::RANK_IN_POOL == 3)
        {
            EXPECT_DOUBLE_EQ(wfcatom(0, 0).real(), -0.035800521528771376);
            EXPECT_DOUBLE_EQ(wfcatom(0, 0).imag(), 0.071704302066073339);
            EXPECT_DOUBLE_EQ(wfcatom(0, wfcpw->npwk[0] - 1).real(), -0.035800589849852141);
            EXPECT_DOUBLE_EQ(wfcatom(0, wfcpw->npwk[0] - 1).imag(), 0.071704271988304383);
            EXPECT_DOUBLE_EQ(wfcatom(7, wfcpw->npwk[0] - 1).real(), -0.02862895801312244);
            EXPECT_DOUBLE_EQ(wfcatom(7, wfcpw->npwk[0] - 1).imag(), 0.045177780166697649);
        }
    }
}

// Test the read_wfc_pw function when the file is not found or wrong type
TEST_F(ReadWfcPwTest, NotFoundFile)
{

#ifdef __MPI
    wfcpw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
    wfcpw->initgrids(5.3233, ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 12, 12, 12);
    wfcpw->initparameters(false, 20, 8, kvec_d);
    wfcpw->setuptransform();
    wfcpw->collect_local_pw();

    ModuleBase::ComplexMatrix wfcatom(PARAM.input.nbands, wfcpw->npwk[0]);

    if(GlobalV::RANK_IN_POOL == 0)
    {
    // dat file
    std::string filename = "notfound.dat";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleIO::read_wfc_pw(filename, wfcpw, 0, nkstot, wfcatom), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,testing::HasSubstr("Can't open file notfound.dat"));


    // txt file
    filename = "notfound.txt";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleIO::read_wfc_pw(filename, wfcpw, 0, nkstot, wfcatom), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,testing::HasSubstr("Can't open file notfound.txt"));

    // other file
    filename = "notfound";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleIO::read_wfc_pw(filename, wfcpw, 0, nkstot, wfcatom), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,testing::HasSubstr("Unknown file type und"));
    }
}

// Test the read_wfc_pw function when nbands is inconsistent
TEST_F(ReadWfcPwTest, InconsistentBands)
{
    if (GlobalV::NPROC_IN_POOL == 1)
    {
        std::string filename = "./support/WAVEFUNC1.dat";

#ifdef __MPI
        wfcpw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
        wfcpw->initgrids(5.3233, ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 12, 12, 12);
        wfcpw->initparameters(false, 20, 8, kvec_d);
        wfcpw->setuptransform();
        wfcpw->collect_local_pw();

        PARAM.input.nbands = 4;
        const int nbasis = wfcpw->npwk[0];
        ModuleBase::ComplexMatrix wfcatom(PARAM.input.nbands, nbasis);
        testing::internal::CaptureStdout();
        EXPECT_EXIT(ModuleIO::read_wfc_pw(filename, wfcpw, 0, nkstot, wfcatom), ::testing::ExitedWithCode(0), "");
        std::string output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("nbands_in = 8"));
        EXPECT_THAT(output, testing::HasSubstr("nbands = 4"));
        EXPECT_THAT(
            output,
            testing::HasSubstr(
                "ikstot_in != ikstot || nkstot_in != nkstot || npwtot_in != npwtot || nbands_in != PARAM.inp.nbands"));
    }
}

// Test the read_wfc_pw function when kevc is inconsistent
TEST_F(ReadWfcPwTest, InconsistentKvec)
{
    if (GlobalV::NPROC_IN_POOL == 1)
    {
        std::string filename = "./support/WAVEFUNC1.dat";

        kvec_d[0] = ModuleBase::Vector3<double>(0.0, 0.0, 1.0);

#ifdef __MPI
        wfcpw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
        wfcpw->initgrids(5.3233, ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 12, 12, 12);
        wfcpw->initparameters(false, 20, 8, kvec_d);
        wfcpw->setuptransform();
        wfcpw->collect_local_pw();

        PARAM.input.nbands = 8;
        const int nbasis = wfcpw->npwk[0];
        ModuleBase::ComplexMatrix wfcatom(PARAM.input.nbands, nbasis);
        testing::internal::CaptureStdout();
        EXPECT_EXIT(ModuleIO::read_wfc_pw(filename, wfcpw, 0, nkstot, wfcatom), ::testing::ExitedWithCode(0), "");
        std::string output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("kvec_in[0] = 0 0 0"));
        EXPECT_THAT(output, testing::HasSubstr("kvec[0] = -1 1 -1"));
        EXPECT_THAT(output, testing::HasSubstr("k vector in file is not the same as the one in memory"));
    }
}

// Test the read_wfc_pw function when lat0 is inconsistent
TEST_F(ReadWfcPwTest, InconsistentLat0)
{
    if (GlobalV::NPROC_IN_POOL == 1)
    {
        std::string filename = "./support/WAVEFUNC1.dat";
        kvec_d[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);

#ifdef __MPI
        wfcpw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
        wfcpw->initgrids(5, ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 12, 12, 12);
        wfcpw->initparameters(false, 20, 8, kvec_d);
        wfcpw->setuptransform();
        wfcpw->collect_local_pw();

        PARAM.input.nbands = 8;
        const int nbasis = wfcpw->npwk[0];
        ModuleBase::ComplexMatrix wfcatom(PARAM.input.nbands, nbasis);
        testing::internal::CaptureStdout();
        EXPECT_EXIT(ModuleIO::read_wfc_pw(filename, wfcpw, 0, nkstot, wfcatom), ::testing::ExitedWithCode(0), "");
        std::string output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("lat0_in = 5.3233"));
        EXPECT_THAT(output, testing::HasSubstr("lat0 = 5"));
        EXPECT_THAT(output, testing::HasSubstr("lat0_in != pw_wfc->lat0 || tpiba_in != pw_wfc->tpiba"));
    }
}

// Test the read_wfc_pw function when G is inconsistent
TEST_F(ReadWfcPwTest, InconsistentG)
{
    if (GlobalV::NPROC_IN_POOL == 1)
    {
        std::string filename = "./support/WAVEFUNC1.dat";
        kvec_d[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);

#ifdef __MPI
        wfcpw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
        wfcpw->initgrids(5.3233, ModuleBase::Matrix3(-0.49, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 12, 12, 12);
        wfcpw->initparameters(false, 20, 8, kvec_d);
        wfcpw->setuptransform();
        wfcpw->collect_local_pw();

        PARAM.input.nbands = 8;
        const int nbasis = wfcpw->npwk[0];
        ModuleBase::ComplexMatrix wfcatom(PARAM.input.nbands, nbasis);
        testing::internal::CaptureStdout();
        EXPECT_EXIT(ModuleIO::read_wfc_pw(filename, wfcpw, 0, nkstot, wfcatom), ::testing::ExitedWithCode(0), "");
        std::string output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("G_in[0] = -1 -1 1\nG_in[1] = 1 1 1\nG_in[2] = -1 1 -1\n"));
        EXPECT_THAT(
            output,
            testing::HasSubstr(
                "G[0] = -1.0101 -1.0101 1.0101\nG[1] = 1.0101 1.0101 0.989899\nG[2] = -1.0101 0.989899 -0.989899\n"));
        EXPECT_THAT(output, testing::HasSubstr("G_in != G"));
    }
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