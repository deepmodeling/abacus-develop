#include "gmock/gmock.h"
#include "gtest/gtest.h"

#define private public
#include "source_io/module_parameter/parameter.h"
#undef private

#include "source_io/module_output/cube_io.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>

#ifdef __MPI
#include "source_basis/module_pw/test/test_tool.h"
#include "mpi.h"
#endif

/**
 * Tested functions:
 *   - ModuleIO::write_cube()
 *   - ModuleIO::read_cube()
 *
 * Correctness tests validate write→readback roundtrip fidelity.
 * Performance tests record serial baseline timing for pre/post optimization comparison.
 *
 * Extension points (after MPI-IO / binary format implemented):
 *   - MPIWriteAllRanksConsistent
 *   - TextVsBinaryConsistency
 *   - Bench_WriteCube_MPIIO_np{N}_{size}
 *   - Bench_WriteCube_BinaryFormat_{size}
 */

// --- helper: build Cube metadata & data ---
struct CubeMeta
{
    std::vector<std::string> comment;
    int natom = 0;
    std::vector<double> origin;
    int nx = 0, ny = 0, nz = 0;
    std::vector<double> dx, dy, dz;
    std::vector<int> atom_type;
    std::vector<double> atom_charge;
    std::vector<std::vector<double>> atom_pos;
    std::vector<double> data;
};

static CubeMeta make_cube_meta(int nx, int ny, int nz)
{
    CubeMeta m;
    m.comment = {"Test cube from ABACUS", "Inner loop is z, followed by y and x"};
    m.natom = 2;
    m.origin = {0.0, 0.0, 0.0};
    m.nx = nx;
    m.ny = ny;
    m.nz = nz;

    double hx = 1.0 / nx;
    double hy = 1.0 / ny;
    double hz = 1.0 / nz;
    m.dx = {hx, 0.0, 0.0};
    m.dy = {0.0, hy, 0.0};
    m.dz = {0.0, 0.0, hz};

    m.atom_type = {1, 8};
    m.atom_charge = {1.0, 6.0};
    m.atom_pos = {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}};

    int nxyz = nx * ny * nz;
    m.data.resize(nxyz);
    int nxy = nx * ny;
    for (int ix = 0; ix < nx; ++ix)
    {
        for (int iy = 0; iy < ny; ++iy)
        {
            for (int iz = 0; iz < nz; ++iz)
            {
                double x = static_cast<double>(ix) / nx;
                double y = static_cast<double>(iy) / ny;
                double z = static_cast<double>(iz) / nz;
                double r2 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) + (z - 0.5) * (z - 0.5);
                m.data[ix * nxy + iy * nz + iz] = std::exp(-r2 * 20.0);
            }
        }
    }
    return m;
}

class WriteCubeTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        meta_ = make_cube_meta(4, 4, 4);
#ifdef __MPI
        tmp_file_ = "test_write_cube_tmp_" + std::to_string(GlobalV::MY_RANK) + ".cube";
#else
        tmp_file_ = "test_write_cube_tmp.cube";
#endif
    }

    void TearDown() override
    {
        std::remove(tmp_file_.c_str());
    }

    CubeMeta meta_;
    std::string tmp_file_;
};

// ===================================================================
// Correctness tests
// ===================================================================

TEST_F(WriteCubeTest, WriteTextCubeAndReadBack)
{
    ModuleIO::write_cube(tmp_file_, meta_.comment, meta_.natom, meta_.origin,
                         meta_.nx, meta_.ny, meta_.nz,
                         meta_.dx, meta_.dy, meta_.dz,
                         meta_.atom_type, meta_.atom_charge, meta_.atom_pos,
                         meta_.data, 6, 6);

    // read back
    std::vector<std::string> cmt;
    int natom = 0;
    std::vector<double> org;
    int nx = 0, ny = 0, nz = 0;
    std::vector<double> dx(3), dy(3), dz(3);
    std::vector<int> atype;
    std::vector<double> acharge;
    std::vector<std::vector<double>> apos;
    std::vector<double> rdata;

    bool ok = ModuleIO::read_cube(tmp_file_, cmt, natom, org, nx, ny, nz,
                                  dx, dy, dz, atype, acharge, apos, rdata);
    ASSERT_TRUE(ok);

    EXPECT_EQ(meta_.natom, natom);
    EXPECT_EQ(meta_.nx, nx);
    EXPECT_EQ(meta_.ny, ny);
    EXPECT_EQ(meta_.nz, nz);
    ASSERT_EQ(meta_.data.size(), rdata.size());

    for (size_t i = 0; i < meta_.data.size(); ++i)
    {
        EXPECT_NEAR(meta_.data[i], rdata[i], 1e-6) << "mismatch at index " << i;
    }
}

TEST_F(WriteCubeTest, ReadbackHeaderCorrect)
{
    ModuleIO::write_cube(tmp_file_, meta_.comment, meta_.natom, meta_.origin,
                         meta_.nx, meta_.ny, meta_.nz,
                         meta_.dx, meta_.dy, meta_.dz,
                         meta_.atom_type, meta_.atom_charge, meta_.atom_pos,
                         meta_.data, 6, 6);

    std::vector<std::string> cmt;
    int natom = 0;
    std::vector<double> org;
    int nx = 0, ny = 0, nz = 0;
    std::vector<double> dx(3), dy(3), dz(3);
    std::vector<int> atype;
    std::vector<double> acharge;
    std::vector<std::vector<double>> apos;
    std::vector<double> rdata;

    ModuleIO::read_cube(tmp_file_, cmt, natom, org, nx, ny, nz,
                        dx, dy, dz, atype, acharge, apos, rdata);

    EXPECT_EQ(cmt[0], meta_.comment[0]);
    EXPECT_EQ(cmt[1], meta_.comment[1]);
    EXPECT_EQ(natom, meta_.natom);
    EXPECT_DOUBLE_EQ(org[0], meta_.origin[0]);
    EXPECT_EQ(nx, meta_.nx);
    EXPECT_EQ(ny, meta_.ny);
    EXPECT_EQ(nz, meta_.nz);
    EXPECT_NEAR(dx[0], meta_.dx[0], 1e-8);
    EXPECT_NEAR(dy[1], meta_.dy[1], 1e-8);
    EXPECT_NEAR(dz[2], meta_.dz[2], 1e-8);
    EXPECT_EQ(atype[0], meta_.atom_type[0]);
    EXPECT_EQ(atype[1], meta_.atom_type[1]);
}

TEST_F(WriteCubeTest, DataLayoutZFastest)
{
    // verify z-fastest indexing: data[ix*nxy + iy*nz + iz]
    int n = 4;
    auto m = make_cube_meta(n, n, n);
    int nxy = n * n;
    for (int ix = 0; ix < n; ++ix)
        for (int iy = 0; iy < n; ++iy)
            for (int iz = 0; iz < n; ++iz)
                m.data[ix * nxy + iy * n + iz] = ix * 100.0 + iy * 10.0 + iz;

    ModuleIO::write_cube(tmp_file_, m.comment, m.natom, m.origin,
                         m.nx, m.ny, m.nz, m.dx, m.dy, m.dz,
                         m.atom_type, m.atom_charge, m.atom_pos, m.data, 6, 6);

    std::vector<double> rdata;
    {
        std::vector<std::string> cmt;
        int natom = 0;
        std::vector<double> org;
        int nx_r = 0, ny_r = 0, nz_r = 0;
        std::vector<double> dx_r(3), dy_r(3), dz_r(3);
        std::vector<int> at;
        std::vector<double> ac;
        std::vector<std::vector<double>> ap;
        ModuleIO::read_cube(tmp_file_, cmt, natom, org, nx_r, ny_r, nz_r,
                            dx_r, dy_r, dz_r, at, ac, ap, rdata);
    }

    for (int ix = 0; ix < n; ++ix)
        for (int iy = 0; iy < n; ++iy)
            for (int iz = 0; iz < n; ++iz)
            {
                int idx = ix * nxy + iy * n + iz;
                double expected = ix * 100.0 + iy * 10.0 + iz;
                EXPECT_NEAR(rdata[idx], expected, 1e-8)
                    << "layout mismatch at (" << ix << "," << iy << "," << iz << ")";
            }
}

TEST_F(WriteCubeTest, ReadCubeFileNotFound)
{
    std::vector<std::string> cmt;
    int natom = 0;
    std::vector<double> org;
    int nx = 0, ny = 0, nz = 0;
    std::vector<double> dx(3), dy(3), dz(3);
    std::vector<int> atype;
    std::vector<double> acharge;
    std::vector<std::vector<double>> apos;
    std::vector<double> rdata;

    bool ok = ModuleIO::read_cube("nonexistent_file.cube", cmt, natom, org,
                                  nx, ny, nz, dx, dy, dz, atype, acharge, apos, rdata);
    EXPECT_FALSE(ok);
}

// ===================================================================
// Serial performance benchmarks (baseline for pre/post comparison)
// ===================================================================

static long long bench_write_cube(const CubeMeta& m, const std::string& fn, int repeat)
{
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int r = 0; r < repeat; ++r)
    {
        ModuleIO::write_cube(fn, m.comment, m.natom, m.origin,
                             m.nx, m.ny, m.nz,
                             m.dx, m.dy, m.dz,
                             m.atom_type, m.atom_charge, m.atom_pos,
                             m.data, 6, 6);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
}

static void bench_report(const std::string& name, long long time_us, int repeat,
                         double data_mb, int nproc)
{
    double ms = static_cast<double>(time_us) / 1000.0 / repeat;
    double mbps = data_mb / (ms / 1000.0);
    printf("[BENCH] %-36s  time=%8.2f ms  throughput=%8.2f MB/s  size=%6.1f MB  np=%d\n",
           name.c_str(), ms, mbps, data_mb, nproc);
}

TEST_F(WriteCubeTest, Bench_WriteCube_SerialText_64)
{
    auto m = make_cube_meta(64, 64, 64);
    std::string fn = "bench_wc_64.cube";
    double data_mb = static_cast<double>(m.data.size() * sizeof(double)) / 1048576.0 + 0.05;
    int repeat = 5;
    long long t = bench_write_cube(m, fn, repeat);
    bench_report("WriteCube_SerialText_64", t, repeat, data_mb, 1);
    std::remove(fn.c_str());
}

TEST_F(WriteCubeTest, Bench_WriteCube_SerialText_128)
{
    auto m = make_cube_meta(128, 128, 128);
    std::string fn = "bench_wc_128.cube";
    double data_mb = static_cast<double>(m.data.size() * sizeof(double)) / 1048576.0 + 0.2;
    int repeat = 3;
    long long t = bench_write_cube(m, fn, repeat);
    bench_report("WriteCube_SerialText_128", t, repeat, data_mb, 1);
    std::remove(fn.c_str());
}

TEST_F(WriteCubeTest, Bench_WriteCube_SerialText_256)
{
    auto m = make_cube_meta(256, 256, 256);
    std::string fn = "bench_wc_256.cube";
    double data_mb = static_cast<double>(m.data.size() * sizeof(double)) / 1048576.0 + 0.8;
    int repeat = 1;
    long long t = bench_write_cube(m, fn, repeat);
    bench_report("WriteCube_SerialText_256", t, repeat, data_mb, 1);
    std::remove(fn.c_str());
}

// ===================================================================
// MPI-IO parallel write tests
// ===================================================================

#ifdef __MPI

TEST_F(WriteCubeTest, MPIWriteBinaryCubeAndReadBack)
{
    // All ranks do MPI-IO write, then rank 0 reads back and verifies
    auto m = make_cube_meta(16, 16, 16);
    std::string fn = "test_mpi_binary.cube";

    ModuleIO::write_cube_mpi(fn, m.comment, m.natom, m.origin,
                             m.nx, m.ny, m.nz,
                             m.dx, m.dy, m.dz,
                             m.atom_type, m.atom_charge, m.atom_pos,
                             m.data, 6, MPI_COMM_WORLD);

    // Rank 0 reads back
    if (GlobalV::MY_RANK == 0)
    {
        std::vector<std::string> cmt;
        int nr = 0;
        std::vector<double> org(3);
        int nx_r = 0, ny_r = 0, nz_r = 0;
        std::vector<double> dxr(3), dyr(3), dzr(3);
        std::vector<int> at;
        std::vector<double> ac;
        std::vector<std::vector<double>> ap;
        std::vector<double> rdata;

        bool ok = ModuleIO::read_cube(fn, cmt, nr, org, nx_r, ny_r, nz_r,
                                      dxr, dyr, dzr, at, ac, ap, rdata);
        ASSERT_TRUE(ok);
        EXPECT_EQ(m.natom, nr);
        EXPECT_EQ(m.nx, nx_r);
        EXPECT_EQ(m.ny, ny_r);
        EXPECT_EQ(m.nz, nz_r);
        ASSERT_EQ(m.data.size(), rdata.size());

        for (size_t i = 0; i < m.data.size(); ++i)
            EXPECT_NEAR(m.data[i], rdata[i], 1e-10) << "mismatch at " << i;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (GlobalV::MY_RANK == 0)
        std::remove(fn.c_str());
}

TEST_F(WriteCubeTest, TextVsBinaryConsistency)
{
    // Write the same data with text and MPI binary, verify both read identically
    auto m = make_cube_meta(10, 10, 10);
    std::string fn_txt = "test_consist_txt.cube";
    std::string fn_bin = "test_consist_bin.cube";

    // Write text (rank 0 only)
    if (GlobalV::MY_RANK == 0)
    {
        ModuleIO::write_cube(fn_txt, m.comment, m.natom, m.origin,
                             m.nx, m.ny, m.nz, m.dx, m.dy, m.dz,
                             m.atom_type, m.atom_charge, m.atom_pos,
                             m.data, 6);
    }

    // Write binary (all ranks via MPI-IO)
    ModuleIO::write_cube_mpi(fn_bin, m.comment, m.natom, m.origin,
                             m.nx, m.ny, m.nz, m.dx, m.dy, m.dz,
                             m.atom_type, m.atom_charge, m.atom_pos,
                             m.data, 6, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    // Rank 0 reads both and compares
    if (GlobalV::MY_RANK == 0)
    {
        auto read_one = [](const std::string& f) {
            std::vector<std::string> cmt;
            int nr = 0;
            std::vector<double> org(3);
            int nx_r = 0, ny_r = 0, nz_r = 0;
            std::vector<double> dxr(3), dyr(3), dzr(3);
            std::vector<int> at;
            std::vector<double> ac;
            std::vector<std::vector<double>> ap;
            std::vector<double> d;
            ModuleIO::read_cube(f, cmt, nr, org, nx_r, ny_r, nz_r,
                                dxr, dyr, dzr, at, ac, ap, d);
            return d;
        };

        auto d_txt = read_one(fn_txt);
        auto d_bin = read_one(fn_bin);

        ASSERT_EQ(d_txt.size(), d_bin.size());
        // Text format uses precision=6 (6 sig digits), so 1e-7 tolerance is appropriate
        for (size_t i = 0; i < d_txt.size(); ++i)
            EXPECT_NEAR(d_txt[i], d_bin[i], 1e-6) << "txt vs bin mismatch at " << i;

        std::remove(fn_txt.c_str());
        std::remove(fn_bin.c_str());
    }
}

TEST_F(WriteCubeTest, Bench_WriteCube_MPIIO_Binary_256)
{
    auto m = make_cube_meta(256, 256, 256);
    std::string fn = "bench_mpi_bin_256.cube";
    double data_mb = static_cast<double>(m.data.size() * sizeof(double)) / 1048576.0 + 0.8;
    int repeat = 1;
    int np = GlobalV::NPROC;

    MPI_Barrier(MPI_COMM_WORLD);
    auto t0 = std::chrono::high_resolution_clock::now();

    for (int r = 0; r < repeat; ++r)
    {
        ModuleIO::write_cube_mpi(fn, m.comment, m.natom, m.origin,
                                 m.nx, m.ny, m.nz,
                                 m.dx, m.dy, m.dz,
                                 m.atom_type, m.atom_charge, m.atom_pos,
                                 m.data, 6, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    auto t1 = std::chrono::high_resolution_clock::now();
    long long t = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

    bench_report("WriteCube_MPIIO_Binary_256", t, repeat, data_mb, np);

    if (GlobalV::MY_RANK == 0)
        std::remove(fn.c_str());
}

#endif // __MPI

int main(int argc, char** argv)
{
#ifdef __MPI
    setupmpi(argc, argv, GlobalV::NPROC, GlobalV::MY_RANK);
    divide_pools(GlobalV::NPROC, GlobalV::MY_RANK, GlobalV::NPROC_IN_POOL,
                 GlobalV::KPAR, GlobalV::MY_POOL, GlobalV::RANK_IN_POOL);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    finishmpi();
#endif
    return result;
}
