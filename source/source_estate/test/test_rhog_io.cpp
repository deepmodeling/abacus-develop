#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_estate/rhog_io.h"
#include "source_base/module_parallel/para_world.h"
#include "source_base/module_parallel/para_tag.h"
#include "source_base/module_parallel/para_bridge.h"
#ifdef __MPI
#include "source_basis/module_pw/test/test_tool.h"
#include "mpi.h"
#endif
#include <fstream>
#include <sstream>

/**
 * - Tested Functions:
 *  - read_rhog()
 *  - write_rhog()
 */

class ReadRhogTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis rhopw;
    std::vector<std::vector<std::complex<double>>> rhog_data;
    std::vector<std::complex<double>*> rhog;
    Parallel::ParaWorld pw_world = Parallel::make_pw_world();
    std::ofstream warning_stream;

    void setup_pw_basis()
    {
#ifdef __MPI
        rhopw.initmpi(pw_world.size(), pw_world.rank(), pw_world.comm());
#endif
        rhopw.initgrids(6.5, ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 120);
        rhopw.initparameters(false, 120);
        rhopw.setuptransform();
        rhopw.collect_local_pw();
    }

    void open_warning(const std::string& path)
    {
        warning_stream.open(path);
    }

    void close_warning()
    {
        if (warning_stream.is_open())
        {
            warning_stream.close();
        }
    }

    std::string read_warning_file(const std::string& path)
    {
        std::ifstream ifs(path);
        std::stringstream ss;
        ss << ifs.rdbuf();
        ifs.close();
        return ss.str();
    }

    virtual void SetUp()
    {
        rhog_data.resize(1, std::vector<std::complex<double>>(1471));
        rhog.push_back(rhog_data[0].data());
    }

    virtual void TearDown()
    {
        close_warning();
    }
};

// Test the read_rhog function with normal file
TEST_F(ReadRhogTest, ReadRhog)
{
    std::string filename = "./support/charge-density.dat";
    setup_pw_basis();

    bool result = elecstate::read_rhog(filename, &rhopw, 1, rhog.data(), pw_world, nullptr);

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
    setup_pw_basis();
    std::string filename = "notfound.txt";

    open_warning("test_read_rhog.txt");
    bool result = elecstate::read_rhog(filename, &rhopw, 1, rhog.data(), pw_world, &warning_stream);
    close_warning();

    std::string expected_content = " elecstate::read_rhog  warning : Can't open file notfound.txt\n";
    EXPECT_FALSE(result);
    EXPECT_EQ(read_warning_file("test_read_rhog.txt"), expected_content);
    std::remove("test_read_rhog.txt");
}

// Test the read_rhog function when gamma_only is inconsistent
TEST_F(ReadRhogTest, InconsistentGammaOnly)
{
    setup_pw_basis();
    std::string filename = "./support/charge-density.dat";
    rhopw.gamma_only = true;
    // Fewer planewaves than the file holds (1471) triggers the
    // "some planewaves in file are not used" warning.
    rhopw.npwtot = 1000;

    open_warning("test_read_rhog.txt");
    bool result = elecstate::read_rhog(filename, &rhopw, 2, rhog.data(), pw_world, &warning_stream);
    close_warning();

    std::string expected_content
        = " elecstate::read_rhog  warning : some planewaves in file are not used\n elecstate::read_rhog  warning : some "
          "spin channels in file are missing\n elecstate::read_rhog  warning : gamma_only read from file is "
          "inconsistent with INPUT\n";

    EXPECT_FALSE(result);
    EXPECT_EQ(read_warning_file("test_read_rhog.txt"), expected_content);
    std::remove("test_read_rhog.txt");
}

// Test the read_rhog function when some planewaves in file are missing
TEST_F(ReadRhogTest, SomePWMissing)
{
    setup_pw_basis();
    std::string filename = "./support/charge-density.dat";
    rhopw.npwtot = 2000;

    open_warning("test_read_rhog.txt");
    bool result = elecstate::read_rhog(filename, &rhopw, 1, rhog.data(), pw_world, &warning_stream);
    close_warning();

    std::string expected_content = " elecstate::read_rhog  warning : some planewaves in file are missing\n";
    EXPECT_TRUE(result);
    EXPECT_EQ(read_warning_file("test_read_rhog.txt"), expected_content);
    std::remove("test_read_rhog.txt");
}

// Test read_rhog with os_warning=nullptr (silent mode, must not crash)
TEST_F(ReadRhogTest, OsNullptrSilent)
{
    std::string filename = "notfound.txt";
    bool result = elecstate::read_rhog(filename, &rhopw, 1, rhog.data(), pw_world, nullptr);
    EXPECT_FALSE(result);
}

// Test write_rhog round-trip: write then read back, verify data consistency
TEST_F(ReadRhogTest, WriteRoundTrip)
{
    setup_pw_basis();

    // initialize some rhog data
    rhog_data[0].assign(rhopw.npw, std::complex<double>(1.5, 2.5));

    std::string tmpfile = "test_rhog_roundtrip.dat";

    // write
    bool write_result = elecstate::write_rhog(
        tmpfile, rhopw.gamma_only, &rhopw, 1,
        ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0),
        rhog.data(), pw_world, nullptr);
    EXPECT_TRUE(write_result);

    // read back into a fresh buffer
    std::vector<std::vector<std::complex<double>>> rhog_read_data(
        1, std::vector<std::complex<double>>(rhopw.npw));
    std::vector<std::complex<double>*> rhog_read;
    rhog_read.push_back(rhog_read_data[0].data());

    bool read_result = elecstate::read_rhog(tmpfile, &rhopw, 1, rhog_read.data(), pw_world, nullptr);
    EXPECT_TRUE(read_result);

    // compare: within MPI precision tolerance
    int diff_count = 0;
    for (int ig = 0; ig < rhopw.npw; ++ig)
    {
        if (std::abs(rhog[0][ig] - rhog_read[0][ig]) > 1e-10)
        {
            ++diff_count;
        }
    }
    EXPECT_EQ(diff_count, 0) << diff_count << " planewave values differ after round-trip";

    std::remove(tmpfile.c_str());
}

// Test write_rhog when the output path is not writable
TEST_F(ReadRhogTest, WriteFileFail)
{
    setup_pw_basis();
    rhog_data[0].assign(rhopw.npw, std::complex<double>(1.0, 0.0));

    // try to write to a directory path (not a file) — should fail
    bool result = elecstate::write_rhog(
        "/tmp", rhopw.gamma_only, &rhopw, 1,
        ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0),
        rhog.data(), pw_world, nullptr);
    EXPECT_FALSE(result);
}

// Test write_rhog with nspin=2, round-trip both channels
TEST_F(ReadRhogTest, WriteRoundTripNspin2)
{
    setup_pw_basis();

    // expand to nspin=2
    rhog_data.resize(2, std::vector<std::complex<double>>(rhopw.npw));
    rhog.clear();
    rhog.push_back(rhog_data[0].data());
    rhog.push_back(rhog_data[1].data());

    // initialize distinct values for each spin channel
    for (int ig = 0; ig < rhopw.npw; ++ig)
    {
        rhog_data[0][ig] = std::complex<double>(1.0 * ig, 0.1 * ig);
        rhog_data[1][ig] = std::complex<double>(2.0 * ig, 0.2 * ig);
    }

    std::string tmpfile = "test_rhog_roundtrip_nspin2.dat";

    // write nspin=2
    bool write_result = elecstate::write_rhog(
        tmpfile, rhopw.gamma_only, &rhopw, 2,
        ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0),
        rhog.data(), pw_world, nullptr);
    EXPECT_TRUE(write_result);

    // read back
    std::vector<std::vector<std::complex<double>>> rhog_read_data(
        2, std::vector<std::complex<double>>(rhopw.npw));
    std::vector<std::complex<double>*> rhog_read;
    rhog_read.push_back(rhog_read_data[0].data());
    rhog_read.push_back(rhog_read_data[1].data());

    bool read_result = elecstate::read_rhog(tmpfile, &rhopw, 2, rhog_read.data(), pw_world, nullptr);
    EXPECT_TRUE(read_result);

    int diff_count = 0;
    for (int is = 0; is < 2; ++is)
    {
        for (int ig = 0; ig < rhopw.npw; ++ig)
        {
            if (std::abs(rhog[is][ig] - rhog_read[is][ig]) > 1e-10)
            {
                ++diff_count;
            }
        }
    }
    EXPECT_EQ(diff_count, 0) << diff_count << " planewave values differ after nspin=2 round-trip";

    std::remove(tmpfile.c_str());
}

// Test write_rhog with nspin=4, round-trip all 4 channels
TEST_F(ReadRhogTest, WriteRoundTripNspin4)
{
    setup_pw_basis();

    rhog_data.resize(4, std::vector<std::complex<double>>(rhopw.npw));
    rhog.clear();
    for (int is = 0; is < 4; ++is)
    {
        rhog.push_back(rhog_data[is].data());
    }

    // initialize distinct values for each spin channel
    for (int is = 0; is < 4; ++is)
    {
        for (int ig = 0; ig < rhopw.npw; ++ig)
        {
            rhog_data[is][ig] = std::complex<double>((is + 1) * 1.0 * ig, (is + 1) * 0.1 * ig);
        }
    }

    std::string tmpfile = "test_rhog_roundtrip_nspin4.dat";

    bool write_result = elecstate::write_rhog(
        tmpfile, rhopw.gamma_only, &rhopw, 4,
        ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0),
        rhog.data(), pw_world, nullptr);
    EXPECT_TRUE(write_result);

    // read back as nspin=4
    std::vector<std::vector<std::complex<double>>> rhog_read_data(
        4, std::vector<std::complex<double>>(rhopw.npw));
    std::vector<std::complex<double>*> rhog_read;
    for (int is = 0; is < 4; ++is)
    {
        rhog_read.push_back(rhog_read_data[is].data());
    }

    bool read_result = elecstate::read_rhog(tmpfile, &rhopw, 4, rhog_read.data(), pw_world, nullptr);
    EXPECT_TRUE(read_result);

    int diff_count = 0;
    for (int is = 0; is < 4; ++is)
    {
        for (int ig = 0; ig < rhopw.npw; ++ig)
        {
            if (std::abs(rhog[is][ig] - rhog_read[is][ig]) > 1e-10)
            {
                ++diff_count;
            }
        }
    }
    EXPECT_EQ(diff_count, 0) << diff_count << " planewave values differ after nspin=4 round-trip";

    std::remove(tmpfile.c_str());
}

// Test the special path L173-181: file nspin=2 read as input nspin=4
// Expected behavior: rhog[0] preserved, rhog[1] and rhog[2] zeroed,
// rhog[3] <- old rhog[1]
TEST_F(ReadRhogTest, ReadRhogNspin2To4SpecialPath)
{
    setup_pw_basis();

    // Step 1: write a nspin=2 binary with known values
    rhog_data.resize(2, std::vector<std::complex<double>>(rhopw.npw));
    rhog.clear();
    rhog.push_back(rhog_data[0].data());
    rhog.push_back(rhog_data[1].data());

    for (int ig = 0; ig < rhopw.npw; ++ig)
    {
        rhog_data[0][ig] = std::complex<double>(10.0 + ig, 0.0);
        rhog_data[1][ig] = std::complex<double>(20.0 + ig, 0.0);
    }

    std::string tmpfile = "test_rhog_nspin2_to_4.dat";

    bool write_result = elecstate::write_rhog(
        tmpfile, rhopw.gamma_only, &rhopw, 2,
        ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0),
        rhog.data(), pw_world, nullptr);
    EXPECT_TRUE(write_result);

    // Step 2: read back as nspin=4 — triggers the L173-181 special path
    std::vector<std::vector<std::complex<double>>> rhog_read_data(
        4, std::vector<std::complex<double>>(rhopw.npw));
    std::vector<std::complex<double>*> rhog_read;
    for (int is = 0; is < 4; ++is)
    {
        rhog_read.push_back(rhog_read_data[is].data());
    }

    bool read_result = elecstate::read_rhog(tmpfile, &rhopw, 4, rhog_read.data(), pw_world, nullptr);
    EXPECT_TRUE(read_result);

    // Verify the special transformation at L173-181:
    // rhog[0]  <- file spin 0
    // rhog[1]  <- ZEROED (was file spin 1, then ZEROS)
    // rhog[2]  <- ZEROED
    // rhog[3]  <- file spin 1 (copied before ZEROS)
    for (int ig = 0; ig < rhopw.npw; ++ig)
    {
        // rhog[0] should match original spin 0
        EXPECT_NEAR(rhog_read_data[0][ig].real(), 10.0 + ig, 1e-10);
        EXPECT_NEAR(rhog_read_data[0][ig].imag(), 0.0, 1e-10);

        // rhog[1] should be zeroed
        EXPECT_NEAR(rhog_read_data[1][ig].real(), 0.0, 1e-10);
        EXPECT_NEAR(rhog_read_data[1][ig].imag(), 0.0, 1e-10);

        // rhog[2] should be zeroed
        EXPECT_NEAR(rhog_read_data[2][ig].real(), 0.0, 1e-10);
        EXPECT_NEAR(rhog_read_data[2][ig].imag(), 0.0, 1e-10);

        // rhog[3] should equal original spin 1 (copied before zero)
        EXPECT_NEAR(rhog_read_data[3][ig].real(), 20.0 + ig, 1e-10);
        EXPECT_NEAR(rhog_read_data[3][ig].imag(), 0.0, 1e-10);
    }

    std::remove(tmpfile.c_str());
}

int main(int argc, char** argv)
{
#ifdef __MPI
    int nproc = 1;
    int myrank = 0;
    int nproc_in_pool = 1;
    int kpar = 1;
    int mypool = 0;
    int rank_in_pool = 0;
    setupmpi(argc, argv, nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    finishmpi();
#endif
    return result;
}
