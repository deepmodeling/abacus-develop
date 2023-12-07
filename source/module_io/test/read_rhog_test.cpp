#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_io/rhog_io.h"
#ifdef __MPI
#include "module_basis/module_pw/test/test_tool.h"
#include "mpi.h"
#endif

/**
 * - Tested Functions:
 *  - read_rhog()
 */

class ReadRhogTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis* rhopw = nullptr;
    std::complex<double>** rhog = nullptr;

    virtual void SetUp()
    {
        rhopw = new ModulePW::PW_Basis;
        rhog = new std::complex<double>*[1];
        rhog[0] = new std::complex<double>[1471];
    }
    virtual void TearDown()
    {
        if (rhopw != nullptr)
            delete rhopw;
        if (rhog[0] != nullptr)
            delete[] rhog[0];
        if (rhog != nullptr)
            delete[] rhog;
    }
};

// Test the read_rhog function
TEST_F(ReadRhogTest, ReadRhog)
{
    std::string filename = "./support/charge-density.dat";
    GlobalV::NSPIN = 1;
#ifdef __MPI
    rhopw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, MPI_COMM_WORLD);
#endif
    rhopw->initgrids(6.5, ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 120);
    rhopw->initparameters(false, 120);
    rhopw->setuptransform();
    rhopw->collect_local_pw();

    bool result = ModuleIO::read_rhog(filename, rhopw, rhog);

    EXPECT_TRUE(result);
    EXPECT_DOUBLE_EQ(rhog[0][0].real(), -1.0304462993299456e-05);
    EXPECT_DOUBLE_EQ(rhog[0][0].imag(), -1.2701788626185278e-13);
    EXPECT_DOUBLE_EQ(rhog[0][1].real(), -0.0003875762482855959);
    EXPECT_DOUBLE_EQ(rhog[0][1].imag(), -4.2556814316812048e-12);
    EXPECT_DOUBLE_EQ(rhog[0][1470].real(), -3.5683133614445107e-05);
    EXPECT_DOUBLE_EQ(rhog[0][1470].imag(), 1.6176615686863767e-12);
}

// Test the read_rhog function when the file is not found
TEST_F(ReadRhogTest, NotFoundFile)
{
    std::string filename = "notfound.txt";

    GlobalV::ofs_warning.open("test_read_rhog.txt");
    bool result = ModuleIO::read_rhog(filename, rhopw, rhog);
    GlobalV::ofs_warning.close();

    std::ifstream ifs_running("test_read_rhog.txt");
    std::stringstream ss;
    ss << ifs_running.rdbuf();
    std::string file_content = ss.str();
    ifs_running.close();

    std::string expected_content = " ModuleIO::read_rhog  warning : Can't open file notfound.txt\n";

    EXPECT_FALSE(result);
    EXPECT_EQ(file_content, expected_content);
    std::remove("test_read_rhog.txt");
}

// Test the read_rhog function when tgamma_only is inconsistent
TEST_F(ReadRhogTest, InconsistentGammaOnly)
{
    std::string filename = "./support/charge-density.dat";
    GlobalV::NSPIN = 2;
    rhopw->gamma_only = true;

    GlobalV::ofs_warning.open("test_read_rhog.txt");
    bool result = ModuleIO::read_rhog(filename, rhopw, rhog);
    GlobalV::ofs_warning.close();

    std::ifstream ifs_running("test_read_rhog.txt");
    std::stringstream ss;
    ss << ifs_running.rdbuf();
    std::string file_content = ss.str();
    ifs_running.close();

    std::string expected_content
        = " ModuleIO::read_rhog  warning : some planewaves in file are not used\n ModuleIO::read_rhog  warning : some "
          "spin channels in file are missing\n ModuleIO::read_rhog  warning : gamma_only read from file is "
          "inconsistent with INPUT\n";

    EXPECT_FALSE(result);
    EXPECT_EQ(file_content, expected_content);
    std::remove("test_read_rhog.txt");
}

// Test the read_rhog function when some planewaves in file are missing
TEST_F(ReadRhogTest, SomePWMissing)
{
    std::string filename = "./support/charge-density.dat";
    GlobalV::NSPIN = 1;
    rhopw->npwtot = 2000;

    GlobalV::ofs_warning.open("test_read_rhog.txt");
    bool result = ModuleIO::read_rhog(filename, rhopw, rhog);
    GlobalV::ofs_warning.close();

    std::ifstream ifs_running("test_read_rhog.txt");
    std::stringstream ss;
    ss << ifs_running.rdbuf();
    std::string file_content = ss.str();
    ifs_running.close();

    std::string expected_content = " ModuleIO::read_rhog  warning : some planewaves in file are missing\n";

    EXPECT_TRUE(result);
    EXPECT_EQ(file_content, expected_content);
    std::remove("test_read_rhog.txt");
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