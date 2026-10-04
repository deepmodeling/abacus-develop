#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_estate/module_charge/chg_rhog_io.h"
#include "source_base/module_parallel/para_world.h"
#include "source_base/module_parallel/para_tag.h"
#include "source_base/module_parallel/para_bridge.h"
#ifdef __MPI
#include "source_basis/module_pw/test/test_tool.h"
#include "mpi.h"
#endif
#include <array>
#include <complex>
#include <cstring>
#include <fstream>
#include <sstream>

/**
 * - Tested Functions:
 *  - read_rhog()
 *  - write_rhog()
 *
 * All binary inputs are generated inside the tests (either via write_rhog
 * round-trip or by hand-crafting the byte layout), so no support/*.dat file
 * is required. The hand-crafted case (ReadLegacyBinaryFormat) pins the
 * on-disk format against silent drift of the writer.
 */

class ReadRhogTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis rhopw;
    std::vector<std::vector<std::complex<double>>> rhog_data;
    std::vector<std::complex<double>*> rhog;
    Parallel::ParaWorld pw_world = Parallel::make_pw_world();
    std::ofstream warning_stream;

    static ModuleBase::Matrix3 latvec()
    {
        return ModuleBase::Matrix3(-0.5, 0.0, 0.5,
                                   0.0, 0.5, 0.5,
                                   -0.5, 0.5, 0.0);
    }

    void setup_pw_basis()
    {
#ifdef __MPI
        rhopw.initmpi(pw_world.size(), pw_world.rank(), pw_world.comm());
#endif
        // Small, explicitly-sized FFT grid: avoids the expensive automatic
        // grid search in initgrids(lat0, latvec, ecut) while still providing
        // a few dozen planewaves for meaningful IO tests.
        rhopw.initgrids(6.5, latvec(), 8, 8, 8);
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

    // Fill rhog_data with distinct, deterministic values for nspin channels.
    void fill_rhog(const int nspin)
    {
        rhog_data.assign(nspin, std::vector<std::complex<double>>(rhopw.npw));
        rhog.clear();
        for (int is = 0; is < nspin; ++is)
        {
            for (int ig = 0; ig < rhopw.npw; ++ig)
            {
                rhog_data[is][ig] = std::complex<double>((is + 1) * 1.0 * ig, (is + 1) * 0.1 * ig);
            }
            rhog.push_back(rhog_data[is].data());
        }
    }

    // Write a binary rhog file by hand (no write_rhog), keeping the exact
    // field order documented in chg_rhog_io.cpp:
    //   [3][gamma_only][npwtot][nspin][3]
    //   [9][b1..b3 (9 doubles)][9]
    //   [3*npwtot][miller ints][3*npwtot]
    //   per spin: [npwtot][complex doubles][npwtot]
    // Only rank 0 writes; call sites must ensure the data lives on rank 0.
    void write_binary_by_hand(const std::string& filename,
                              const int gamma_only_flag,
                              const int nspin_file,
                              const std::vector<std::array<int, 3>>& miller,
                              const std::vector<std::vector<std::complex<double>>>& values) const
    {
        if (pw_world.rank() != 0)
        {
            return;
        }
        const int npw_file = static_cast<int>(miller.size());
        std::ofstream ofs(filename, std::ios::binary);

        int size = 3;
        ofs.write(reinterpret_cast<const char*>(&size), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&gamma_only_flag), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&npw_file), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&nspin_file), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&size), sizeof(int));

        size = 9;
        const ModuleBase::Matrix3 GT = latvec().Inverse();
        const double b[9] = {GT.e11, GT.e12, GT.e13, GT.e21, GT.e22, GT.e23, GT.e31, GT.e32, GT.e33};
        ofs.write(reinterpret_cast<const char*>(&size), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(b), 9 * sizeof(double));
        ofs.write(reinterpret_cast<const char*>(&size), sizeof(int));

        size = 3 * npw_file;
        ofs.write(reinterpret_cast<const char*>(&size), sizeof(int));
        for (const auto& m : miller)
        {
            ofs.write(reinterpret_cast<const char*>(m.data()), 3 * sizeof(int));
        }
        ofs.write(reinterpret_cast<const char*>(&size), sizeof(int));

        size = npw_file;
        for (int is = 0; is < nspin_file; ++is)
        {
            ofs.write(reinterpret_cast<const char*>(&size), sizeof(int));
            ofs.write(reinterpret_cast<const char*>(values[is].data()), npw_file * sizeof(std::complex<double>));
            ofs.write(reinterpret_cast<const char*>(&size), sizeof(int));
        }
        ofs.close();
    }

    virtual void SetUp()
    {
        // Buffers are sized to the actual basis in each test via fill_rhog
        // after setup_pw_basis(); allocate a minimal placeholder here.
        rhog_data.resize(1, std::vector<std::complex<double>>(1));
        rhog.push_back(rhog_data[0].data());
    }

    virtual void TearDown()
    {
        close_warning();
    }
};

// Round-trip: write known data with write_rhog, read back, verify values.
// Replaces the old support/charge-density.dat based test.
TEST_F(ReadRhogTest, ReadRhog)
{
    setup_pw_basis();
    fill_rhog(1);

    const std::string tmpfile = "test_rhog_read.dat";
    bool write_result = module_charge::write_rhog(tmpfile, rhopw.gamma_only, &rhopw, 1, latvec(), rhog.data(), pw_world, nullptr);
    ASSERT_TRUE(write_result);

    std::vector<std::complex<double>> read_back(rhopw.npw);
    std::complex<double>* read_ptr = read_back.data();
    bool read_result = module_charge::read_rhog(tmpfile, &rhopw, 1, &read_ptr, pw_world, nullptr);
    ASSERT_TRUE(read_result);

    for (int ig = 0; ig < rhopw.npw; ++ig)
    {
        EXPECT_DOUBLE_EQ(read_back[ig].real(), rhog_data[0][ig].real());
        EXPECT_DOUBLE_EQ(read_back[ig].imag(), rhog_data[0][ig].imag());
    }
    std::remove(tmpfile.c_str());
}

// Pin the on-disk binary format: hand-craft a file (without write_rhog)
// and check read_rhog maps values to the correct G-vectors.
TEST_F(ReadRhogTest, ReadLegacyBinaryFormat)
{
    setup_pw_basis();
    fill_rhog(1);

    // Use two well-separated G-vectors: the Gamma point and the first
    // non-zero planewave, so the Miller-index -> ig mapping is exercised.
    std::vector<std::array<int, 3>> miller(2);
    miller[0] = {0, 0, 0};
    const ModuleBase::Vector3<double> g1 = rhopw.gdirect[1];
    miller[1] = {static_cast<int>(g1.x), static_cast<int>(g1.y), static_cast<int>(g1.z)};

    std::vector<std::vector<std::complex<double>>> values(1, std::vector<std::complex<double>>(2));
    values[0][0] = std::complex<double>(1.5, -0.5);
    values[0][1] = std::complex<double>(2.5, 3.5);

    const std::string tmpfile = "test_rhog_legacy.dat";
    write_binary_by_hand(tmpfile, 0, 1, miller, values);

    bool read_result = module_charge::read_rhog(tmpfile, &rhopw, 1, rhog.data(), pw_world, nullptr);
    ASSERT_TRUE(read_result);

    // ig_gge0 is the index of the Gamma point in the local basis.
    EXPECT_DOUBLE_EQ(rhog_data[0][rhopw.ig_gge0].real(), 1.5);
    EXPECT_DOUBLE_EQ(rhog_data[0][rhopw.ig_gge0].imag(), -0.5);
    // The second entry lands on gdirect[1]; find its ig via the Miller index.
    bool found = false;
    for (int ig = 0; ig < rhopw.npw; ++ig)
    {
        if (static_cast<int>(rhopw.gdirect[ig].x) == miller[1][0] &&
            static_cast<int>(rhopw.gdirect[ig].y) == miller[1][1] &&
            static_cast<int>(rhopw.gdirect[ig].z) == miller[1][2])
        {
            EXPECT_DOUBLE_EQ(rhog_data[0][ig].real(), 2.5);
            EXPECT_DOUBLE_EQ(rhog_data[0][ig].imag(), 3.5);
            found = true;
        }
    }
    EXPECT_TRUE(found);
    std::remove(tmpfile.c_str());
}

// Test the read_rhog function when the file is not found
TEST_F(ReadRhogTest, NotFoundFile)
{
    setup_pw_basis();
    fill_rhog(1);
    std::string filename = "notfound.txt";

    open_warning("test_read_rhog.txt");
    bool result = module_charge::read_rhog(filename, &rhopw, 1, rhog.data(), pw_world, &warning_stream);
    close_warning();

    std::string expected_content = " module_charge::read_rhog  warning : Can't open file notfound.txt\n";
    EXPECT_FALSE(result);
    EXPECT_EQ(read_warning_file("test_read_rhog.txt"), expected_content);
    std::remove("test_read_rhog.txt");
}

// Test the read_rhog function when gamma_only is inconsistent
TEST_F(ReadRhogTest, InconsistentGammaOnly)
{
    setup_pw_basis();
    fill_rhog(1);

    // Self-generate a file with gamma_only=0, nspin=1.
    const std::string tmpfile = "test_rhog_gamma.dat";
    bool write_result = module_charge::write_rhog(tmpfile, rhopw.gamma_only, &rhopw, 1, latvec(), rhog.data(), pw_world, nullptr);
    ASSERT_TRUE(write_result);

    // Flip gamma_only and shrink npwtot to trigger the warning branches.
    rhopw.gamma_only = true;
    rhopw.npwtot -= 1;

    open_warning("test_read_rhog.txt");
    bool result = module_charge::read_rhog(tmpfile, &rhopw, 2, rhog.data(), pw_world, &warning_stream);
    close_warning();

    std::string expected_content
        = " module_charge::read_rhog  warning : some planewaves in file are not used\n module_charge::read_rhog  warning : some "
          "spin channels in file are missing\n module_charge::read_rhog  warning : gamma_only read from file is "
          "inconsistent with INPUT\n";

    EXPECT_FALSE(result);
    EXPECT_EQ(read_warning_file("test_read_rhog.txt"), expected_content);
    std::remove(tmpfile.c_str());
    std::remove("test_read_rhog.txt");
}

// Test the read_rhog function when some planewaves in file are missing
TEST_F(ReadRhogTest, SomePWMissing)
{
    setup_pw_basis();
    fill_rhog(1);

    const std::string tmpfile = "test_rhog_missing.dat";
    bool write_result = module_charge::write_rhog(tmpfile, rhopw.gamma_only, &rhopw, 1, latvec(), rhog.data(), pw_world, nullptr);
    ASSERT_TRUE(write_result);

    // Pretend the basis holds more planewaves than the file.
    rhopw.npwtot += 1;

    open_warning("test_read_rhog.txt");
    bool result = module_charge::read_rhog(tmpfile, &rhopw, 1, rhog.data(), pw_world, &warning_stream);
    close_warning();

    std::string expected_content = " module_charge::read_rhog  warning : some planewaves in file are missing\n";
    EXPECT_TRUE(result);
    EXPECT_EQ(read_warning_file("test_read_rhog.txt"), expected_content);
    std::remove(tmpfile.c_str());
    std::remove("test_read_rhog.txt");
}

// Test read_rhog with os_warning=nullptr (silent mode, must not crash)
TEST_F(ReadRhogTest, OsNullptrSilent)
{
    std::string filename = "notfound.txt";
    bool result = module_charge::read_rhog(filename, &rhopw, 1, rhog.data(), pw_world, nullptr);
    EXPECT_FALSE(result);
}

// Test write_rhog round-trip: write then read back, verify data consistency
TEST_F(ReadRhogTest, WriteRoundTrip)
{
    setup_pw_basis();
    fill_rhog(1);

    std::string tmpfile = "test_rhog_roundtrip.dat";

    bool write_result = module_charge::write_rhog(tmpfile, rhopw.gamma_only, &rhopw, 1, latvec(), rhog.data(), pw_world, nullptr);
    EXPECT_TRUE(write_result);

    std::vector<std::vector<std::complex<double>>> rhog_read_data(1, std::vector<std::complex<double>>(rhopw.npw));
    std::vector<std::complex<double>*> rhog_read;
    rhog_read.push_back(rhog_read_data[0].data());

    bool read_result = module_charge::read_rhog(tmpfile, &rhopw, 1, rhog_read.data(), pw_world, nullptr);
    EXPECT_TRUE(read_result);

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
    fill_rhog(1);

    // try to write to a directory path (not a file) — should fail
    bool result = module_charge::write_rhog("/tmp", rhopw.gamma_only, &rhopw, 1, latvec(), rhog.data(), pw_world, nullptr);
    EXPECT_FALSE(result);
}

// Test write_rhog with nspin=2, round-trip both channels
TEST_F(ReadRhogTest, WriteRoundTripNspin2)
{
    setup_pw_basis();
    fill_rhog(2);

    std::string tmpfile = "test_rhog_roundtrip_nspin2.dat";

    bool write_result = module_charge::write_rhog(tmpfile, rhopw.gamma_only, &rhopw, 2, latvec(), rhog.data(), pw_world, nullptr);
    EXPECT_TRUE(write_result);

    std::vector<std::vector<std::complex<double>>> rhog_read_data(2, std::vector<std::complex<double>>(rhopw.npw));
    std::vector<std::complex<double>*> rhog_read;
    rhog_read.push_back(rhog_read_data[0].data());
    rhog_read.push_back(rhog_read_data[1].data());

    bool read_result = module_charge::read_rhog(tmpfile, &rhopw, 2, rhog_read.data(), pw_world, nullptr);
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
    fill_rhog(4);

    std::string tmpfile = "test_rhog_roundtrip_nspin4.dat";

    bool write_result = module_charge::write_rhog(tmpfile, rhopw.gamma_only, &rhopw, 4, latvec(), rhog.data(), pw_world, nullptr);
    EXPECT_TRUE(write_result);

    std::vector<std::vector<std::complex<double>>> rhog_read_data(4, std::vector<std::complex<double>>(rhopw.npw));
    std::vector<std::complex<double>*> rhog_read;
    for (int is = 0; is < 4; ++is)
    {
        rhog_read.push_back(rhog_read_data[is].data());
    }

    bool read_result = module_charge::read_rhog(tmpfile, &rhopw, 4, rhog_read.data(), pw_world, nullptr);
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
    fill_rhog(2);

    // Override with simple, distinct values for this scenario.
    for (int ig = 0; ig < rhopw.npw; ++ig)
    {
        rhog_data[0][ig] = std::complex<double>(10.0 + ig, 0.0);
        rhog_data[1][ig] = std::complex<double>(20.0 + ig, 0.0);
    }

    std::string tmpfile = "test_rhog_nspin2_to_4.dat";

    bool write_result = module_charge::write_rhog(tmpfile, rhopw.gamma_only, &rhopw, 2, latvec(), rhog.data(), pw_world, nullptr);
    EXPECT_TRUE(write_result);

    std::vector<std::vector<std::complex<double>>> rhog_read_data(4, std::vector<std::complex<double>>(rhopw.npw));
    std::vector<std::complex<double>*> rhog_read;
    for (int is = 0; is < 4; ++is)
    {
        rhog_read.push_back(rhog_read_data[is].data());
    }

    bool read_result = module_charge::read_rhog(tmpfile, &rhopw, 4, rhog_read.data(), pw_world, nullptr);
    EXPECT_TRUE(read_result);

    for (int ig = 0; ig < rhopw.npw; ++ig)
    {
        EXPECT_NEAR(rhog_read_data[0][ig].real(), 10.0 + ig, 1e-10);
        EXPECT_NEAR(rhog_read_data[0][ig].imag(), 0.0, 1e-10);
        EXPECT_NEAR(rhog_read_data[1][ig].real(), 0.0, 1e-10);
        EXPECT_NEAR(rhog_read_data[1][ig].imag(), 0.0, 1e-10);
        EXPECT_NEAR(rhog_read_data[2][ig].real(), 0.0, 1e-10);
        EXPECT_NEAR(rhog_read_data[2][ig].imag(), 0.0, 1e-10);
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
