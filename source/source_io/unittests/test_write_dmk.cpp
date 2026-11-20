#include <complex>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "source_io/module_dm/write_dmk.h"
#include "source_io/module_hs/dense_text_body.h"
#include "source_base/global_variable.h"
#include "prepare_unitcell.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_base/module_external/scalapack_connector.h"

#ifdef __MPI
#include "mpi.h"
#include "source_base/module_parallel/para_mpi_func.h"
#endif

namespace
{
#ifdef __MPI
/**
 * @brief RAII guard for a test-scoped MPI communicator created with
 *        MPI_Comm_split; frees the communicator on destruction.
 *
 * The project parallel helpers deliberately do not own communicator
 * lifetimes (see ParaWorld), so this local guard encapsulates the
 * split/free pair inside the test that needs it.
 */
class SplitCommGuard
{
public:
    /**
     * @brief Split MPI_COMM_WORLD by color.
     * @param color Processes sharing a color join the same new communicator.
     * @param key   Rank ordering key inside the new communicator.
     */
    SplitCommGuard(const int color, const int key) : comm_(MPI_COMM_NULL)
    {
        MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm_);
    }

    ~SplitCommGuard()
    {
        if (comm_ != MPI_COMM_NULL)
        {
            MPI_Comm_free(&comm_);
        }
    }

    SplitCommGuard(const SplitCommGuard&) = delete;
    SplitCommGuard& operator=(const SplitCommGuard&) = delete;

    /// @brief Wrapped communicator (owned by this guard).
    MPI_Comm get() const
    {
        return comm_;
    }

private:
    MPI_Comm comm_;
};

/// @brief Synchronize all ranks on MPI_COMM_WORLD via the parallel helpers.
void world_barrier()
{
    const Parallel::ParaWorld world = Parallel::ParaWorld::make_mpi("test_write_dmk", MPI_COMM_WORLD);
    Parallel::barrier(world);
}
#endif

/// @brief Rank that owns the written files: MPI rank 0, or the only process
///        in a serial build.
bool on_io_rank()
{
#ifdef __MPI
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank == 0;
#else
    return true;
#endif
}
} // namespace

Magnetism::Magnetism() {
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism() { }

/************************************************
 *  unit test of read_dmk and write_dmk
 ***********************************************/

/**
 * - Tested Functions:
 *   - read_dmk()
 *     - the function to read density matrix K from file
 *     - the serial version without MPI
 *   - write_dmk()
 *     - the function to write density matrix K to file
 *     - the serial version without MPI
 */

/// @brief Fixture owning the Si UnitCell shared by most tests.
class DMKTest : public testing::Test
{
protected:
    void SetUp() override
    {
        UcellTestPrepare utp = UcellTestLib["Si"];
        ucell = utp.SetUcellInfo();
    }

    void TearDown() override
    {
        delete ucell;
    }

    UnitCell* ucell = nullptr;
};

/// @brief Element-wise near-equality check for two identically shaped buffers.
template <typename T>
void expect_dmk_eq(const std::vector<std::vector<T>>& actual,
                   const std::vector<std::vector<T>>& expected)
{
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t ib = 0; ib < actual.size(); ++ib)
    {
        ASSERT_EQ(actual[ib].size(), expected[ib].size());
        for (std::size_t ie = 0; ie < actual[ib].size(); ++ie)
        {
            EXPECT_NEAR(actual[ib][ie], expected[ib][ie], 1e-8);
        }
    }
}

/// @brief Complex overload: checks real and imaginary parts separately.
void expect_dmk_eq(const std::vector<std::vector<std::complex<double>>>& actual,
                   const std::vector<std::vector<std::complex<double>>>& expected)
{
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t ib = 0; ib < actual.size(); ++ib)
    {
        ASSERT_EQ(actual[ib].size(), expected[ib].size());
        for (std::size_t ie = 0; ie < actual[ib].size(); ++ie)
        {
            EXPECT_NEAR(actual[ib][ie].real(), expected[ib][ie].real(), 1e-8);
            EXPECT_NEAR(actual[ib][ie].imag(), expected[ib][ie].imag(), 1e-8);
        }
    }
}

void init_pv(int nlocal, Parallel_Blkcyclic& pv)
{
#ifdef __MPI
        pv.init(nlocal, nlocal, 1, MPI_COMM_WORLD, false);
#else
        pv.set_serial(nlocal, nlocal);
#endif             
}

template <typename T>
void gen_dmk(std::vector<std::vector<T>>& dmk, std::vector<double>& efs,  int nspin, int nk, int nlocal, Parallel_Blkcyclic& pv)
{
    int myrank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
    std::vector<std::vector<T>> dmk_global(nspin * nk, std::vector<T>(nlocal * nlocal, T(0.0)));
    if (myrank == 0)
    {
        for (int i = 0; i < nspin * nk; i++)
        {
            for (int j = 0; j < nlocal * nlocal; j++)
            {
                std::complex<double> value = std::complex<double>(1.0 * i + 0.1 * j, 0.1 * i + 1.0 * j);
                dmk_global[i][j] = reinterpret_cast<T*>(&value)[0];
            }
        }
    }
#ifdef __MPI
    Parallel_Blkcyclic pv_global;
    pv_global.init(nlocal, nlocal, nlocal, MPI_COMM_WORLD, false);
    dmk.resize(nspin * nk, std::vector<T>(pv.size(), 0.0));
    for (int i = 0; i < nspin * nk; i++)
    {
        Cpxgemr2d(nlocal,
                  nlocal,
                  dmk_global[i].data(),
                  1,
                  1,
                  const_cast<int*>(pv_global.desc()),
                  dmk[i].data(),
                  1,
                  1,
                  const_cast<int*>(pv.desc()),
                  pv.ictxt());
    }
#else
    dmk = dmk_global;
#endif

    efs.resize(nspin, 0.0);
    for (int i = 0; i < nspin; i++)
    {
        efs[i] = 0.1 * i;
    }
}

/// @brief Fill the single-k-point K_Vectors metadata used by the writers.
void init_kv(K_Vectors& kv, const double weight)
{
    kv.set_nkstot(1);
    kv.set_nkstot_nospin(1);
    kv.set_nks(1);
    kv.set_spin_mult(2);
    kv.kvec_c.resize(1);
    kv.kvec_c[0].x = 0.0;
    kv.kvec_c[0].y = 0.0;
    kv.kvec_c[0].z = 0.0;
    kv.kvec_d.resize(1);
    kv.kvec_d[0].x = 0.0;
    kv.kvec_d[0].y = 0.0;
    kv.kvec_d[0].z = 0.0;
    kv.wk.resize(1);
    kv.wk[0] = weight;
    kv.isk.resize(1);
    kv.isk[0] = 0;
    kv.kc_done = true;
    kv.kd_done = true;
}

/// @brief Tokenize the numbers of the matrix body written at the end of a DMK
///        file.
/// @details The body is the trailing run of lines carrying complex pairs, so it
///          is collected from the end of the file; taking the first "(" instead
///          would hit the UnitCell header, and parsing the whole file would
///          trip over the non-numeric metadata. Parentheses and commas become
///          whitespace, so the value order (real, imaginary, real, ...) is
///          preserved.
std::vector<double> read_numbers(const std::string& fname)
{
    std::ifstream ifs(fname.c_str());
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(ifs, line))
    {
        lines.push_back(line);
    }
    std::size_t first_body_line = lines.size();
    while (first_body_line > 0)
    {
        const std::string& candidate = lines[first_body_line - 1];
        if (candidate.find('(') == std::string::npos || candidate.find(')') == std::string::npos)
        {
            break;
        }
        --first_body_line;
    }
    std::string text;
    for (std::size_t i = first_body_line; i < lines.size(); ++i)
    {
        text += lines[i];
        text += " ";
    }
    for (std::size_t i = 0; i < text.size(); ++i)
    {
        const char c = text[i];
        if (c == '(' || c == ')' || c == ',')
        {
            text[i] = ' ';
        }
    }
    std::istringstream iss(text);
    std::vector<double> values;
    double value = 0.0;
    while (iss >> value)
    {
        values.push_back(value);
    }
    return values;
}

/// @brief Fill the local blocks of a single-k-point DM with a row-major probe,
///        stored column-major as a column-major eigensolver would leave it.
template <typename T>
void fill_local_dmk(const std::vector<T>& probe, const int nlocal,
                    const Parallel_Blkcyclic& pv, std::vector<std::vector<T>>& dmk)
{
    dmk.assign(1, std::vector<T>(pv.size(), T(0)));
    for (int i = 0; i < nlocal; ++i)
    {
        const int ir = pv.irloc(i);
        if (ir < 0)
        {
            continue;
        }
        for (int j = 0; j < nlocal; ++j)
        {
            const int ic = pv.icloc(j);
            if (ic >= 0)
            {
                dmk[0][ir + ic * pv.nrloc()] = probe[i * nlocal + j];
            }
        }
    }
}

/// @brief Fill the K_Vectors metadata for nspin spin channels and nk k-points.
/// @details write_dmk reads one k-point coordinate, weight and index per file,
///          so a multi-k round trip needs real k-point arrays instead of the
///          single-point metadata set up by init_kv.
void init_kv_multi(K_Vectors& kv, const int nspin, const int nk)
{
    kv.set_nkstot(nspin * nk);
    kv.set_nkstot_nospin(nk);
    kv.set_nks(nspin * nk);
    kv.set_spin_mult(2);
    kv.kvec_c.resize(nk);
    kv.kvec_d.resize(nk);
    kv.wk.resize(nk);
    kv.isk.resize(nk);
    for (int ik = 0; ik < nk; ++ik)
    {
        kv.kvec_c[ik].x = 0.1 * ik;
        kv.kvec_c[ik].y = 0.0;
        kv.kvec_c[ik].z = 0.0;
        kv.kvec_d[ik].x = 0.1 * ik;
        kv.kvec_d[ik].y = 0.0;
        kv.kvec_d[ik].z = 0.0;
        kv.wk[ik] = 1.0 / nk;
        kv.isk[ik] = 0;
    }
    kv.kc_done = true;
    kv.kd_done = true;
}

/// @brief Copy a file inside the scratch directory.
/// @details Used to build a file whose header contradicts the name the reader
///          looks for, which cannot be produced through write_dmk itself.
bool copy_file(const std::string& src, const std::string& dst)
{
    std::ifstream ifs(src.c_str());
    if (!ifs.good())
    {
        return false;
    }
    std::ofstream ofs(dst.c_str());
    if (!ofs.good())
    {
        return false;
    }
    ofs << ifs.rdbuf();
    return ofs.good();
}

/// @brief Drop the final line of a text file to simulate a truncated write.
bool drop_last_line(const std::string& fname)
{
    std::ifstream ifs(fname.c_str());
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(ifs, line))
    {
        lines.push_back(line);
    }
    ifs.close();
    if (lines.empty())
    {
        return false;
    }
    std::ofstream ofs(fname.c_str());
    for (std::size_t i = 0; i + 1 < lines.size(); ++i)
    {
        ofs << lines[i] << "\n";
    }
    ofs.close();
    return true;
}

TEST_F(DMKTest, GenFileName) {
    bool gamma_only = true;
    int ispin = 0;
    int nspin = 2;
    int ik = 0;
    int istep = 0;
    std::string fname = ModuleIO::dmk_gen_fname(gamma_only, ispin, nspin, ik, istep);
    EXPECT_EQ(fname, "dms1g1_nao.txt");

    ispin = 1;

    fname = ModuleIO::dmk_gen_fname(gamma_only, ispin, nspin, ik, istep);
    EXPECT_EQ(fname, "dms2g1_nao.txt");

    ispin = 0;
    gamma_only = false;    

    fname = ModuleIO::dmk_gen_fname(gamma_only, ispin, nspin, ik, istep);
    EXPECT_EQ(fname, "dmk1s1g1_nao.txt");

    ispin = 1;
    ik = 1;

    fname = ModuleIO::dmk_gen_fname(gamma_only, ispin, nspin, ik, istep);
    EXPECT_EQ(fname, "dmk2s2g1_nao.txt");
};


TEST_F(DMKTest,WriteDMK) {
    int nspin = 2;
    int nk = 1;
    int nk_multik = 2;
    int nlocal = 20;
    std::vector<std::vector<double>> dmk;
    std::vector<std::vector<std::complex<double>>> dmk_multik;
    Parallel_Blkcyclic pv;
    std::vector<double> efs;
    init_pv(nlocal, pv);

    gen_dmk(dmk, efs, nspin, nk, nlocal, pv);
    gen_dmk(dmk_multik, efs, nspin, nk_multik, nlocal, pv);
    const std::string out_dir = "./";

    const int istep = -1;
    K_Vectors kv;
    init_kv(kv, 1.0);

    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, 3, efs, ucell, pv, /*column_major=*/true,
                        out_dir, istep, 0, ofs_log);
    ModuleIO::write_dmk(dmk_multik, kv, 3, efs, ucell, pv, /*column_major=*/true,
                        out_dir, istep, 0, ofs_log);

    if (on_io_rank())
    {
        // {filename, expected Fermi energy, check trailing newline}
        const struct { const char* name; const char* fermi; bool check_newline; } cases[] = {
            {"dms1_nao.txt",  "0 # Fermi energy in Ry",   false},
            {"dms2_nao.txt",  "0.1 # Fermi energy in Ry", false},
            {"dmk1s1_nao.txt","0 # Fermi energy in Ry",   true},
            {"dmk2s1_nao.txt","0 # Fermi energy in Ry",   false},
            {"dmk1s2_nao.txt","0.1 # Fermi energy in Ry", false},
            {"dmk2s2_nao.txt","0.1 # Fermi energy in Ry", false},
        };
        for (const auto& c : cases)
        {
            std::ifstream ifs(c.name);
            const std::string str((std::istreambuf_iterator<char>(ifs)),
                                  std::istreambuf_iterator<char>());
            ifs.close();
            EXPECT_THAT(str, testing::HasSubstr(c.fermi)) << c.name;
            EXPECT_THAT(str, testing::HasSubstr("20 # number of localized basis")) << c.name;
            // the matrix body must end with a newline so appended records are not glued
            if (c.check_newline)
            {
                EXPECT_EQ(str.back(), '\n') << c.name;
            }
            remove(c.name);
        }
    }
};

// The DM(k) text body is row-major even when the local blocks are stored
// column-major (as left behind by genelpa/elpa/scalapack_gvx). A non-symmetric
// complex probe makes an accidental transpose visible.
TEST_F(DMKTest, BodyIsRowMajorForColumnMajorStorage) {
    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    // Row-major probe: (0,0)=1, (0,1)=2+3i, (1,0)=4+5i, (1,1)=6.
    const std::vector<std::complex<double>> probe
        = {std::complex<double>(1.0, 0.0), std::complex<double>(2.0, 3.0),
           std::complex<double>(4.0, 5.0), std::complex<double>(6.0, 0.0)};

    // Column-major storage of the same matrix: local[ir + ic * nrloc].
    std::vector<std::vector<std::complex<double>>> dmk;
    fill_local_dmk(probe, nlocal, pv, dmk);

    std::vector<double> efs(nspin, 0.0);
    K_Vectors kv;
    init_kv(kv, 1.0);

    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, 8, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    if (on_io_rank())
    {
        // nk = 1 with a multi-k (complex) DM -> dmk1_nao.txt.
        const std::string fn = out_dir + "dmk1_nao.txt";
        const std::vector<double> values = read_numbers(fn);
        const std::vector<double> expected = {1.0, 0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.0};
        ASSERT_EQ(values.size(), expected.size());
        for (std::size_t i = 0; i < expected.size(); ++i)
        {
            EXPECT_NEAR(values[i], expected[i], 1e-8);
        }
        remove(fn.c_str());
    }
}

// out_app_flag = 0: one fixed filename, truncated and overwritten each call.
TEST_F(DMKTest, Mode0FixedNameOverwrite) {
    int nspin = 2;
    int nk = 1;
    int nlocal = 20;
    std::vector<std::vector<double>> dmk;
    std::vector<double> efs;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);
    gen_dmk(dmk, efs, nspin, nk, nlocal, pv);

    K_Vectors kv;
    init_kv(kv, 1.0);

    const std::string out_dir = "./";
    const int out_app_flag = 0;
    std::ofstream ofs_log;

    // Two different isteps -> the same fixed filename, overwritten each call.
    ModuleIO::write_dmk(dmk, kv, /*precision=*/3, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, out_app_flag, ofs_log);
    ModuleIO::write_dmk(dmk, kv, /*precision=*/3, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/1, out_app_flag, ofs_log);

    if (on_io_rank())
    {
        std::ifstream ifs("dms1_nao.txt");
        std::string str((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>());
        ifs.close();
        // The second write replaced the first record: only step 2 remains.
        EXPECT_THAT(str, testing::HasSubstr("Ionic_Step 2"));
        EXPECT_THAT(str, testing::Not(testing::HasSubstr("Ionic_Step 1")));
        // g-suffixed files must NOT exist.
        ifs.open("dms1g1_nao.txt");
        EXPECT_FALSE(ifs.good());
        ifs.open("dms1g2_nao.txt");
        EXPECT_FALSE(ifs.good());
        remove("dms1_nao.txt");
        remove("dms2_nao.txt");
    }
}

// out_app_flag = 2: per-step g{step} filename; directory routing is handled by
// set_global_dir and is covered by io_sys_var_test. Mode 1 (append) is covered
// by ReadDMKMode1FrameSelection, which would fail on a wrong write.
TEST_F(DMKTest, Mode2PerStepFile) {
    int nspin = 2;
    int nk = 1;
    int nlocal = 20;
    std::vector<std::vector<double>> dmk;
    std::vector<double> efs;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);
    gen_dmk(dmk, efs, nspin, nk, nlocal, pv);

    K_Vectors kv;
    init_kv(kv, 1.0);

    const std::string out_dir = "./";
    const int out_app_flag = 2;
    std::ofstream ofs_log;

    ModuleIO::write_dmk(dmk, kv, /*precision=*/3, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/2, out_app_flag, ofs_log);

    if (on_io_rank())
    {
        std::ifstream ifs("dms1g3_nao.txt");
        EXPECT_TRUE(ifs.good());
        ifs.close();
        remove("dms1g3_nao.txt");
    }
}

// Round trip: write a column-major DM(k) and read it back with the same
// storage order. A transposed read would return the complex conjugate of the
// Hermitian probe, i.e. the imaginary off-diagonal elements would flip sign.
TEST_F(DMKTest, ReadDMKRoundTripComplex) {
    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    // Row-major probe: (0,0)=1, (0,1)=2+3i, (1,0)=4+5i, (1,1)=6.
    const std::vector<std::complex<double>> probe
        = {std::complex<double>(1.0, 0.0), std::complex<double>(2.0, 3.0),
           std::complex<double>(4.0, 5.0), std::complex<double>(6.0, 0.0)};

    std::vector<std::vector<std::complex<double>>> dmk;
    fill_local_dmk(probe, nlocal, pv, dmk);

    std::vector<double> efs(nspin, 0.0);
    K_Vectors kv;
    init_kv(kv, 1.0);

    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    std::vector<std::vector<std::complex<double>>> dmk_read;
    EXPECT_TRUE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                   /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));
    expect_dmk_eq(dmk_read, dmk);

    if (on_io_rank())
    {
        remove((out_dir + "dmk1_nao.txt").c_str());
    }
}

TEST_F(DMKTest, ReadDMKRoundTripGamma) {
    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    // Row-major probe: 1, 2 / 3, 4.
    const std::vector<double> probe = {1.0, 2.0, 3.0, 4.0};
    std::vector<std::vector<double>> dmk;
    fill_local_dmk(probe, nlocal, pv, dmk);

    std::vector<double> efs(nspin, 0.0);
    K_Vectors kv;
    init_kv(kv, 1.0);

    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    std::vector<std::vector<double>> dmk_read;
    EXPECT_TRUE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                   /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));
    expect_dmk_eq(dmk_read, dmk);

    if (on_io_rank())
    {
        remove((out_dir + "dm_nao.txt").c_str());
    }
}

TEST_F(DMKTest, ReadDMKRoundTripMultipleMatrices) {
    const int nspin = 2;
    const int nk = 2;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    K_Vectors kv;
    init_kv_multi(kv, nspin, nk);

    std::vector<std::vector<std::complex<double>>> dmk(nspin * nk);
    for (int index = 0; index < nspin * nk; ++index)
    {
        std::vector<std::complex<double>> probe(nlocal * nlocal);
        for (int element = 0; element < nlocal * nlocal; ++element)
        {
            probe[element] = std::complex<double>(0.25 * (index + 1) * (element + 1),
                                                  0.1 * index - 0.2 * element);
        }
        std::vector<std::vector<std::complex<double>>> one_dmk;
        fill_local_dmk(probe, nlocal, pv, one_dmk);
        dmk[index] = one_dmk[0];
    }

    const std::string out_dir = testing::TempDir();
    const std::vector<double> efs(nspin, 0.0);
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    std::vector<std::vector<std::complex<double>>> dmk_read;
    EXPECT_TRUE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                   /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));
    expect_dmk_eq(dmk_read, dmk);

    if (on_io_rank())
    {
        remove((out_dir + "dmk1s1_nao.txt").c_str());
        remove((out_dir + "dmk2s1_nao.txt").c_str());
        remove((out_dir + "dmk1s2_nao.txt").c_str());
        remove((out_dir + "dmk2s2_nao.txt").c_str());
    }
}

// A missing file must be reported cooperatively instead of aborting a single
// rank, so the caller can decide how to continue.
TEST_F(DMKTest, ReadDMKMissingFileReturnsFalse) {
    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    K_Vectors kv;
    init_kv(kv, 1.0);

    std::vector<std::vector<double>> dmk(1, std::vector<double>(nlocal * nlocal, 7.0));
    const std::vector<std::vector<double>> original = dmk;
    std::ofstream ofs_log;
    const std::string missing_dir = testing::TempDir() + "dmk_no_such_dir/";
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, missing_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk, ofs_log));
    EXPECT_EQ(dmk, original);
}

#ifdef __MPI
/// @brief Verify that split descriptor communicators keep their read status independent.
TEST_F(DMKTest, ReadDMKSplitCommunicatorStatus)
{
    const Parallel::ParaWorld world = Parallel::ParaWorld::make_mpi("test_write_dmk", MPI_COMM_WORLD);
    const int world_rank = world.rank();
    const int world_size = world.size();
    if (world_size < 2)
    {
        GTEST_SKIP() << "This test requires at least two MPI ranks";
    }

    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;

    Parallel_Blkcyclic pv_world;
    init_pv(nlocal, pv_world);
    const std::vector<double> probe = {1.0, 2.0, 3.0, 4.0};
    std::vector<std::vector<double>> dmk;
    fill_local_dmk(probe, nlocal, pv_world, dmk);

    K_Vectors kv;
    init_kv(kv, 1.0);
    const std::vector<double> efs(nspin, 0.0);
    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv_world, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);
    world_barrier();

    const int color = world_rank % 2;
    const SplitCommGuard split_guard(color, world_rank);

    int local_pass = 1;
    {
        Parallel_Blkcyclic pv_split;
        pv_split.init(nlocal, nlocal, 1, split_guard.get(), false);
        const std::vector<double> sentinel(pv_split.size(), 7.0);
        std::vector<std::vector<double>> dmk_read(1, sentinel);
        const bool expect_success = (color == 0);
        const std::string read_dir = expect_success ? out_dir : out_dir + "dmk_split_missing_dir/";
        const bool read_success = ModuleIO::read_dmk(nspin, nk, kv, pv_split, /*column_major=*/true,
                                                     read_dir, /*istep=*/0, /*out_app_flag=*/0,
                                                     dmk_read, ofs_log);
        if (read_success != expect_success)
        {
            local_pass = 0;
        }

        if (expect_success)
        {
            std::vector<std::vector<double>> expected;
            fill_local_dmk(probe, nlocal, pv_split, expected);
            if (dmk_read != expected)
            {
                local_pass = 0;
            }
        }
        else if (dmk_read.size() != 1 || dmk_read[0] != sentinel)
        {
            local_pass = 0;
        }
    }

    int global_pass = local_pass;
    Parallel::reduce_min(global_pass, world);
    EXPECT_EQ(global_pass, 1) << "Split communicators exchanged read status or returned incorrect data";

    if (world_rank == 0)
    {
        std::remove((out_dir + "dm_nao.txt").c_str());
    }
    world_barrier();
}
#endif

// out_app_flag = 2: the ionic step is part of the filename, so the reader has
// to request the matching g{istep+1} file and reject a frame whose marker
// contradicts the name.
TEST_F(DMKTest, ReadDMKMode2PerStepFile) {
    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    const std::vector<std::complex<double>> probe
        = {std::complex<double>(1.0, 0.0), std::complex<double>(2.0, 3.0),
           std::complex<double>(4.0, 5.0), std::complex<double>(6.0, 0.0)};
    std::vector<std::vector<std::complex<double>>> dmk;
    fill_local_dmk(probe, nlocal, pv, dmk);

    std::vector<double> efs(nspin, 0.0);
    K_Vectors kv;
    init_kv(kv, 1.0);

    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/2, /*out_app_flag=*/2, ofs_log);

    std::vector<std::vector<std::complex<double>>> dmk_read;
    EXPECT_TRUE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                   /*istep=*/2, /*out_app_flag=*/2, dmk_read, ofs_log));
    expect_dmk_eq(dmk_read, dmk);

    // A different step maps to another file name, which does not exist.
    std::vector<std::vector<std::complex<double>>> dmk_other_step;
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/2, dmk_other_step, ofs_log));

    if (on_io_rank())
    {
        remove((out_dir + "dmk1g3_nao.txt").c_str());
    }
}

// out_app_flag = 1: a fixed-name file accumulates one record per ionic step,
// so the reader has to select the requested frame instead of always taking the
// first one.
TEST_F(DMKTest, ReadDMKMode1FrameSelection) {
    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    const std::vector<std::complex<double>> probe_first
        = {std::complex<double>(1.0, 0.0), std::complex<double>(2.0, 3.0),
           std::complex<double>(4.0, 5.0), std::complex<double>(6.0, 0.0)};
    const std::vector<std::complex<double>> probe_second
        = {std::complex<double>(7.0, 0.0), std::complex<double>(8.0, 1.0),
           std::complex<double>(9.0, 2.0), std::complex<double>(10.0, 0.0)};
    std::vector<std::vector<std::complex<double>>> dmk_first;
    std::vector<std::vector<std::complex<double>>> dmk_second;
    fill_local_dmk(probe_first, nlocal, pv, dmk_first);
    fill_local_dmk(probe_second, nlocal, pv, dmk_second);

    std::vector<double> efs(nspin, 0.0);
    K_Vectors kv;
    init_kv(kv, 1.0);

    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    // Step 0 truncates the file and step 1 appends the second record.
    ModuleIO::write_dmk(dmk_first, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/1, ofs_log);
    ModuleIO::write_dmk(dmk_second, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/1, /*out_app_flag=*/1, ofs_log);

    std::vector<std::vector<std::complex<double>>> dmk_read;
    EXPECT_TRUE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                   /*istep=*/0, /*out_app_flag=*/1, dmk_read, ofs_log));
    expect_dmk_eq(dmk_read, dmk_first);

    dmk_read.clear();
    EXPECT_TRUE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                   /*istep=*/1, /*out_app_flag=*/1, dmk_read, ofs_log));
    expect_dmk_eq(dmk_read, dmk_second);

    // A frame that was never appended must be reported, not silently replaced
    // by the last available record.
    std::vector<std::vector<std::complex<double>>> dmk_missing_frame;
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/5, /*out_app_flag=*/1, dmk_missing_frame, ofs_log));

    if (on_io_rank())
    {
        remove((out_dir + "dmk1_nao.txt").c_str());
    }
}

/// @brief Check that every block of a DM(k) buffer still holds the sentinel.
bool all_equal_to(const std::vector<std::vector<std::complex<double>>>& dmk,
                  const std::complex<double>& sentinel)
{
    for (std::size_t ib = 0; ib < dmk.size(); ++ib)
    {
        for (std::size_t ie = 0; ie < dmk[ib].size(); ++ie)
        {
            if (dmk[ib][ie] != sentinel)
            {
                return false;
            }
        }
    }
    return true;
}

// A failure on a later (spin, k) point must discard the matrices that were
// already read; otherwise the caller would keep a half-updated DM(k) and run
// the next SCF from a mixture of two restarts.
TEST_F(DMKTest, ReadDMKPartialFailureKeepsCallerData) {
    const int nspin = 1;
    const int nk = 2;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    K_Vectors kv;
    init_kv_multi(kv, nspin, nk);

    std::vector<std::vector<std::complex<double>>> dmk;
    for (int ik = 0; ik < nk; ++ik)
    {
        std::vector<std::complex<double>> probe(nlocal * nlocal,
                                                std::complex<double>(ik + 1.0, 0.5 * ik));
        std::vector<std::vector<std::complex<double>>> one_dmk;
        fill_local_dmk(probe, nlocal, pv, one_dmk);
        dmk.push_back(one_dmk[0]);
    }

    const std::vector<double> efs(nspin, 0.0);
    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    // Removing the last matrix makes the first one read successfully and the
    // second one fail, so the failure happens after partial progress.
    if (on_io_rank())
    {
        std::remove((out_dir + "dmk2_nao.txt").c_str());
    }

    const std::complex<double> sentinel(7.0, 7.0);
    std::vector<std::vector<std::complex<double>>> dmk_read(
        nspin * nk, std::vector<std::complex<double>>(pv.nrloc() * pv.ncloc(), sentinel));

    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));
    EXPECT_TRUE(all_equal_to(dmk_read, sentinel));

    if (on_io_rank())
    {
        std::remove((out_dir + "dmk1_nao.txt").c_str());
    }
}

// A spin-down file whose header claims spin 1 must be rejected instead of being
// loaded as the down channel.
TEST_F(DMKTest, ReadDMKRejectsSpinMismatch) {
    const int nspin = 2;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    K_Vectors kv;
    init_kv_multi(kv, nspin, nk);

    const std::vector<std::complex<double>> probe_up
        = {std::complex<double>(1.0, 0.0), std::complex<double>(2.0, 3.0),
           std::complex<double>(4.0, 5.0), std::complex<double>(6.0, 0.0)};
    const std::vector<std::complex<double>> probe_down
        = {std::complex<double>(7.0, 0.0), std::complex<double>(8.0, 1.0),
           std::complex<double>(9.0, 2.0), std::complex<double>(10.0, 0.0)};
    std::vector<std::vector<std::complex<double>>> dmk;
    std::vector<std::vector<std::complex<double>>> probe_block;
    fill_local_dmk(probe_up, nlocal, pv, probe_block);
    dmk.push_back(probe_block[0]);
    fill_local_dmk(probe_down, nlocal, pv, probe_block);
    dmk.push_back(probe_block[0]);

    std::vector<double> efs(nspin, 0.0);
    efs[1] = 0.1;
    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    // Overwrite the spin-down file with the spin-up one: the content now
    // contradicts both the name and the requested spin channel.
    if (on_io_rank())
    {
        EXPECT_TRUE(copy_file(out_dir + "dmk1s1_nao.txt", out_dir + "dmk1s2_nao.txt"));
    }

    const std::complex<double> sentinel(7.0, 7.0);
    std::vector<std::vector<std::complex<double>>> dmk_read(
        nspin * nk, std::vector<std::complex<double>>(pv.nrloc() * pv.ncloc(), sentinel));
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));
    EXPECT_TRUE(all_equal_to(dmk_read, sentinel));

    if (on_io_rank())
    {
        std::remove((out_dir + "dmk1s1_nao.txt").c_str());
        std::remove((out_dir + "dmk1s2_nao.txt").c_str());
    }
}

// k-point metadata: a different k-mesh and a file holding the wrong k-point
// must both be rejected instead of being used for the requested index.
TEST_F(DMKTest, ReadDMKRejectsKPointMismatch) {
    const int nspin = 1;
    const int nk = 2;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    K_Vectors kv;
    init_kv_multi(kv, nspin, nk);

    std::vector<std::vector<std::complex<double>>> dmk;
    for (int ik = 0; ik < nk; ++ik)
    {
        std::vector<std::complex<double>> probe(nlocal * nlocal,
                                                std::complex<double>(ik + 1.0, 0.5 * ik));
        std::vector<std::vector<std::complex<double>>> one_dmk;
        fill_local_dmk(probe, nlocal, pv, one_dmk);
        dmk.push_back(one_dmk[0]);
    }

    const std::vector<double> efs(nspin, 0.0);
    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    // A restart with a different total number of k-points must not reuse these
    // files even though the per-k filename happens to match.
    K_Vectors kv_other_mesh;
    init_kv_multi(kv_other_mesh, nspin, nk);
    kv_other_mesh.set_nkstot(nspin * nk + 1);
    std::vector<std::vector<std::complex<double>>> dmk_read;
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv_other_mesh, pv, /*column_major=*/true,
                                    out_dir, /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    // The second k-point's record copied over the first name: the header index
    // no longer matches the requested k-point.
    if (on_io_rank())
    {
        EXPECT_TRUE(copy_file(out_dir + "dmk2_nao.txt", out_dir + "dmk1_nao.txt"));
    }
    dmk_read.clear();
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    if (on_io_rank())
    {
        std::remove((out_dir + "dmk1_nao.txt").c_str());
        std::remove((out_dir + "dmk2_nao.txt").c_str());
    }
}

// k-point coordinates: a mesh of the same size but with reordered or shifted
// k-vectors must be rejected, because cal_dmr pairs the read DM(k) with the
// current KPT phases and would silently build a wrong density.
TEST_F(DMKTest, ReadDMKRejectsKVecMismatch) {
    const int nspin = 1;
    const int nk = 2;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    K_Vectors kv;
    init_kv_multi(kv, nspin, nk);

    std::vector<std::vector<std::complex<double>>> dmk;
    for (int ik = 0; ik < nk; ++ik)
    {
        std::vector<std::complex<double>> probe(nlocal * nlocal,
                                                std::complex<double>(ik + 1.0, 0.5 * ik));
        std::vector<std::vector<std::complex<double>>> one_dmk;
        fill_local_dmk(probe, nlocal, pv, one_dmk);
        dmk.push_back(one_dmk[0]);
    }

    const std::vector<double> efs(nspin, 0.0);
    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    // Same nkstot and same indices, but the k-vectors are swapped between the
    // two files: the stored coordinate no longer matches the requested ik.
    K_Vectors kv_swapped;
    init_kv_multi(kv_swapped, nspin, nk);
    kv_swapped.kvec_d[0].x = 0.1;
    kv_swapped.kvec_d[1].x = 0.0;
    std::vector<std::vector<std::complex<double>>> dmk_read;
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv_swapped, pv, /*column_major=*/true,
                                    out_dir, /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    // A mesh shifted by a fractional amount with the same count and indices
    // must also be rejected.
    K_Vectors kv_shifted;
    init_kv_multi(kv_shifted, nspin, nk);
    kv_shifted.kvec_d[0].x += 0.05;
    kv_shifted.kvec_d[1].x += 0.05;
    dmk_read.clear();
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv_shifted, pv, /*column_major=*/true,
                                    out_dir, /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    // A shift by an integer reciprocal-lattice vector leaves e^{ik·R}
    // unchanged, so the same DM(k) is still valid and must be accepted.
    K_Vectors kv_wrapped;
    init_kv_multi(kv_wrapped, nspin, nk);
    kv_wrapped.kvec_d[0].x += 1.0;
    kv_wrapped.kvec_d[1].x -= 2.0;
    dmk_read.clear();
    EXPECT_TRUE(ModuleIO::read_dmk(nspin, nk, kv_wrapped, pv, /*column_major=*/true,
                                   out_dir, /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    if (on_io_rank())
    {
        std::remove((out_dir + "dmk1_nao.txt").c_str());
        std::remove((out_dir + "dmk2_nao.txt").c_str());
    }
}

// A basis size that differs from the current calculation must be rejected, not
// silently truncated or over-read.
TEST_F(DMKTest, ReadDMKRejectsBasisMismatch) {
    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    K_Vectors kv;
    init_kv_multi(kv, nspin, nk);

    const std::vector<std::complex<double>> probe
        = {std::complex<double>(1.0, 0.0), std::complex<double>(2.0, 3.0),
           std::complex<double>(4.0, 5.0), std::complex<double>(6.0, 0.0)};
    std::vector<std::vector<std::complex<double>>> dmk;
    fill_local_dmk(probe, nlocal, pv, dmk);

    const std::vector<double> efs(nspin, 0.0);
    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    // Read the 2x2 record with a 3-orbital distribution.
    Parallel_Blkcyclic pv_wide;
    init_pv(nlocal + 1, pv_wide);
    std::vector<std::vector<std::complex<double>>> dmk_read(
        nspin * nk, std::vector<std::complex<double>>(pv_wide.nrloc() * pv_wide.ncloc(),
                                                      std::complex<double>(7.0, 7.0)));
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv_wide, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    if (on_io_rank())
    {
        std::remove((out_dir + "dmk1_nao.txt").c_str());
    }
}

// Truncated input: a file cut inside the header and a complete file whose last
// body line was lost must both fail instead of yielding a partial matrix.
TEST_F(DMKTest, ReadDMKRejectsTruncatedFile) {
    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    K_Vectors kv;
    init_kv_multi(kv, nspin, nk);

    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;

    // Cut inside the header: only the spin count survives.
    if (on_io_rank())
    {
        std::ofstream ofs((out_dir + "dm_nao.txt").c_str());
        ofs << "1\n";
        ofs.close();
    }
    std::vector<std::vector<double>> dmk_read;
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    // A complete record whose final line was lost: the body is short.
    const std::vector<double> probe = {1.0, 2.0, 3.0, 4.0};
    std::vector<std::vector<double>> dmk;
    fill_local_dmk(probe, nlocal, pv, dmk);
    const std::vector<double> efs(nspin, 0.0);
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    if (on_io_rank())
    {
        EXPECT_TRUE(drop_last_line(out_dir + "dm_nao.txt"));
    }
    dmk_read.clear();
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    if (on_io_rank())
    {
        std::remove((out_dir + "dm_nao.txt").c_str());
    }
}

/// @brief Replace a specific line of a text file with new content.
/// @param fname Path of the file to rewrite (must exist).
/// @param target Zero-based line index to replace; out-of-range returns false.
/// @param new_content The replacement text, written with a trailing newline.
bool rewrite_line(const std::string& fname, std::size_t target, const std::string& new_content)
{
    std::ifstream ifs(fname.c_str());
    if (!ifs.good())
    {
        return false;
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(ifs, line))
    {
        lines.push_back(line);
    }
    if (target >= lines.size())
    {
        return false;
    }
    lines[target] = new_content;
    std::ofstream ofs(fname.c_str());
    if (!ofs.good())
    {
        return false;
    }
    for (std::size_t i = 0; i < lines.size(); ++i)
    {
        ofs << lines[i] << "\n";
    }
    return ofs.good();
}

// UnitCell block integrity: a file whose per-type atom-count line cannot be
// parsed as integers used to leave the void skipper at the subsequent
// coordinate rows. A small (e.g. nlocal = 2) matrix read would then consume
// those coordinates and the reader would hand back a bogus density matrix.
// The skipper must fail fast instead, so read_dmk reports parse failure.
TEST_F(DMKTest, ReadDMKRejectsMalformedAtomCount) {
    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    // Row-major probe: 1, 2 / 3, 4.
    const std::vector<double> probe = {1.0, 2.0, 3.0, 4.0};
    std::vector<std::vector<double>> dmk;
    fill_local_dmk(probe, nlocal, pv, dmk);

    const std::vector<double> efs(nspin, 0.0);
    K_Vectors kv;
    init_kv(kv, 1.0);

    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    if (on_io_rank())
    {
        // write_dmk emits a 14-line header (Ionic_Step through the trailing
        // blank separator) followed by the UnitCell block written by write_ucell.
        // For a canonical Si single cell (ntype = 1, na = 2) the per-type atom
        // counts sit at global line index 20:
        //   14: latName, 15: lat0, 16-18: latvec rows, 19: atom labels,
        //   20: per-type counts " 2", 21: "Direct", 22-23: coordinates.
        // Replace the counts with a non-numeric token so read_ucell has to
        // bail out with total_natom == 0. With the old void skipper that
        // would leave the stream at the coordinate rows, whose six doubles
        // are enough to satisfy a 2x2 gamma matrix read.
        const std::string fn = out_dir + "dm_nao.txt";
        ASSERT_TRUE(rewrite_line(fn, /*target=*/20, "X # malformed count\n"));
    }

    const std::vector<std::vector<double>> sentinel(
        1, std::vector<double>(pv.nrloc() * pv.ncloc(), 7.0));
    std::vector<std::vector<double>> dmk_read = sentinel;
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));
    // The caller buffer must not be touched on failure.
    EXPECT_EQ(dmk_read, sentinel);

    if (on_io_rank())
    {
        std::remove((out_dir + "dm_nao.txt").c_str());
    }
}

// Mode 0 uses a fixed filename for every step, so an interrupted write can
// leave files from different ionic steps on disk. The reader must reject such
// a mixture instead of building a restart density from two different steps.
TEST_F(DMKTest, ReadDMKMode0RejectsMixedSteps) {
    const int nspin = 1;
    const int nk = 2;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    K_Vectors kv;
    init_kv_multi(kv, nspin, nk);

    std::vector<std::vector<std::complex<double>>> dmk;
    for (int ik = 0; ik < nk; ++ik)
    {
        std::vector<std::complex<double>> probe(nlocal * nlocal,
                                                std::complex<double>(ik + 1.0, 0.5 * ik));
        std::vector<std::vector<std::complex<double>>> one_dmk;
        fill_local_dmk(probe, nlocal, pv, one_dmk);
        dmk.push_back(one_dmk[0]);
    }

    const std::vector<double> efs(nspin, 0.0);
    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    // Bump the step marker of the second file only, mimicking an interrupted
    // mode-0 write where dmk2 was refreshed but dmk1 was not.
    if (on_io_rank())
    {
        ASSERT_TRUE(rewrite_line(out_dir + "dmk2_nao.txt", /*target=*/0, "Ionic_Step 2"));
    }

    const std::complex<double> sentinel(7.0, 7.0);
    std::vector<std::vector<std::complex<double>>> dmk_read(
        nspin * nk, std::vector<std::complex<double>>(pv.nrloc() * pv.ncloc(), sentinel));
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));
    // The caller buffer must stay untouched on failure.
    EXPECT_TRUE(all_equal_to(dmk_read, sentinel));

    if (on_io_rank())
    {
        std::remove((out_dir + "dmk1_nao.txt").c_str());
        std::remove((out_dir + "dmk2_nao.txt").c_str());
    }
}

// Legacy files written before the Ionic_Step marker existed must stay
// readable in mode 0: without a marker there is nothing to cross-check.
TEST_F(DMKTest, ReadDMKMode0AcceptsFilesWithoutStepMarker) {
    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    const std::vector<double> probe = {1.0, 2.0, 3.0, 4.0};
    std::vector<std::vector<double>> dmk;
    fill_local_dmk(probe, nlocal, pv, dmk);

    const std::vector<double> efs(nspin, 0.0);
    K_Vectors kv;
    init_kv(kv, 1.0);

    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    // Strip the "Ionic_Step" marker line to mimic a pre-marker writer.
    if (on_io_rank())
    {
        const std::string fn = out_dir + "dm_nao.txt";
        std::ifstream ifs(fn.c_str());
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(ifs, line))
        {
            if (line.find("Ionic_Step") == std::string::npos)
            {
                lines.push_back(line);
            }
        }
        ifs.close();
        std::ofstream ofs(fn.c_str());
        for (const std::string& out_line : lines)
        {
            ofs << out_line << "\n";
        }
    }

    std::vector<std::vector<double>> dmk_read;
    EXPECT_TRUE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                   /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));
    expect_dmk_eq(dmk_read, dmk);

    if (on_io_rank())
    {
        std::remove((out_dir + "dm_nao.txt").c_str());
    }
}

// The flat index shared by the readers and the writers: a row-major buffer
// stores (i, j) at i * nlocal + j, a column-major one at j * nlocal + i. All
// LCAO ks_solver choices are column-major, so this mapping is what keeps the
// file row-major for every supported solver.
TEST_F(DMKTest, DenseBufferIndexSelectsLayout) {
    EXPECT_EQ(ModuleIO::dense_buffer_index(0, 0, 3, /*column_major=*/false), 0u);
    EXPECT_EQ(ModuleIO::dense_buffer_index(0, 1, 3, /*column_major=*/false), 1u);
    EXPECT_EQ(ModuleIO::dense_buffer_index(1, 0, 3, /*column_major=*/false), 3u);
    EXPECT_EQ(ModuleIO::dense_buffer_index(2, 1, 3, /*column_major=*/false), 7u);
    EXPECT_EQ(ModuleIO::dense_buffer_index(0, 1, 3, /*column_major=*/true), 3u);
    EXPECT_EQ(ModuleIO::dense_buffer_index(1, 0, 3, /*column_major=*/true), 1u);
    EXPECT_EQ(ModuleIO::dense_buffer_index(1, 2, 3, /*column_major=*/true), 7u);
}

// The header must state the conj-first element convention, and read_dmk must
// transparently skip that comment line when the file is read back.
TEST_F(DMKTest, HeaderStatesConjFirstConvention) {
    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    // Non-symmetric complex probe: a conj-last (i.e. transposed) read of the
    // Hermitian-interpreted body would flip the off-diagonal imaginary signs.
    const std::vector<std::complex<double>> probe
        = {std::complex<double>(1.0, 0.0), std::complex<double>(2.0, 3.0),
           std::complex<double>(4.0, 5.0), std::complex<double>(6.0, 0.0)};
    std::vector<std::vector<std::complex<double>>> dmk;
    fill_local_dmk(probe, nlocal, pv, dmk);

    std::vector<double> efs(nspin, 0.0);
    K_Vectors kv;
    init_kv(kv, 1.0);

    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    if (on_io_rank())
    {
        const std::string fn = out_dir + "dmk1_nao.txt";
        std::ifstream ifs(fn.c_str());
        const std::string str((std::istreambuf_iterator<char>(ifs)),
                              std::istreambuf_iterator<char>());
        ifs.close();
        EXPECT_THAT(str, testing::HasSubstr(
            "# matrix element (i,j) = sum_n w_n * conj(c_n,i) * c_n,j"));
    }

    // The '#' comment line inside the header must not disturb the round trip.
    std::vector<std::vector<std::complex<double>>> dmk_read;
    EXPECT_TRUE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                   /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));
    expect_dmk_eq(dmk_read, dmk);

    if (on_io_rank())
    {
        remove((out_dir + "dmk1_nao.txt").c_str());
    }
}

// With column_major = false the reader and the writer must still agree on the
// same layout, otherwise a restart would read a transposed matrix. No current
// LCAO solver takes this branch, so this only pins their mutual consistency.
TEST_F(DMKTest, ReadDMKRoundTripRowMajorFlag) {
    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    const std::vector<std::complex<double>> probe
        = {std::complex<double>(1.0, 0.0), std::complex<double>(2.0, 3.0),
           std::complex<double>(4.0, 5.0), std::complex<double>(6.0, 0.0)};
    std::vector<std::vector<std::complex<double>>> dmk;
    fill_local_dmk(probe, nlocal, pv, dmk);

    const std::vector<double> efs(nspin, 0.0);
    K_Vectors kv;
    init_kv(kv, 1.0);

    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/false,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    std::vector<std::vector<std::complex<double>>> dmk_read;
    EXPECT_TRUE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/false, out_dir,
                                   /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));
    expect_dmk_eq(dmk_read, dmk);

    if (on_io_rank())
    {
        std::remove((out_dir + "dmk1_nao.txt").c_str());
    }
}

// The reader must reject non-finite matrix components: strtod reports an
// overflowing literal such as "1e999" only through errno while consuming the
// full token, and both strtod and the stream extractor accept the literal
// spellings "inf"/"nan" without flagging a failure. Such a component would
// otherwise be installed as the restart density.

/// @brief Overwrite the first matrix-body token of a written DM(k) file.
/// @details The complex body is the only section whose tokens contain '(';
///          for the real body of a 2x2 matrix the whole body is the final
///          line (8 values per line). Both fixtures used below keep the body
///          free of comment markers, which the header and the UnitCell block
///          never are on their numeric lines.
bool patch_first_body_token(const std::string& fname,
                            const std::string& token,
                            const bool complex_body)
{
    std::ifstream ifs(fname.c_str());
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(ifs, line))
    {
        lines.push_back(line);
    }
    ifs.close();
    if (lines.empty())
    {
        return false;
    }

    std::size_t target = lines.size();
    if (complex_body)
    {
        for (std::size_t i = 0; i < lines.size(); ++i)
        {
            if (lines[i].find('(') != std::string::npos)
            {
                target = i;
                break;
            }
        }
    }
    else
    {
        target = lines.size() - 1;
    }
    if (target == lines.size())
    {
        return false;
    }

    std::string& current = lines[target];
    const std::size_t first = current.find_first_not_of(" \t");
    if (first == std::string::npos)
    {
        return false;
    }
    const std::size_t end = current.find_first_of(" \t", first);
    current.replace(first, end == std::string::npos ? std::string::npos : end - first, token);

    std::ofstream ofs(fname.c_str());
    for (const std::string& out_line : lines)
    {
        ofs << out_line << "\n";
    }
    return true;
}

TEST_F(DMKTest, ReadDMKRejectsOverflowingComplexComponent) {
    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    const std::vector<std::complex<double>> probe
        = {std::complex<double>(1.0, 0.0), std::complex<double>(2.0, 3.0),
           std::complex<double>(4.0, 5.0), std::complex<double>(6.0, 0.0)};
    std::vector<std::vector<std::complex<double>>> dmk;
    fill_local_dmk(probe, nlocal, pv, dmk);

    const std::vector<double> efs(nspin, 0.0);
    K_Vectors kv;
    init_kv(kv, 1.0);

    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    // strtod consumes "1e999" fully and only reports the overflow via ERANGE.
    if (on_io_rank())
    {
        EXPECT_TRUE(patch_first_body_token(out_dir + "dmk1_nao.txt", "(1e999,0)", true));
    }
    std::vector<std::vector<std::complex<double>>> dmk_read;
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    // The literal "nan" parses without any errno flag and must be rejected too.
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);
    if (on_io_rank())
    {
        EXPECT_TRUE(patch_first_body_token(out_dir + "dmk1_nao.txt", "(nan,0)", true));
    }
    dmk_read.clear();
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    if (on_io_rank())
    {
        std::remove((out_dir + "dmk1_nao.txt").c_str());
    }
}

TEST_F(DMKTest, ReadDMKRejectsNonFiniteRealComponent) {
    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    const std::vector<double> probe = {1.0, 2.0, 3.0, 4.0};
    std::vector<std::vector<double>> dmk;
    fill_local_dmk(probe, nlocal, pv, dmk);

    const std::vector<double> efs(nspin, 0.0);
    K_Vectors kv;
    init_kv(kv, 1.0);

    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    // The stream extractor accepts the literal "inf" without setting failbit.
    if (on_io_rank())
    {
        EXPECT_TRUE(patch_first_body_token(out_dir + "dm_nao.txt", "inf", false));
    }
    std::vector<std::vector<double>> dmk_read;
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    if (on_io_rank())
    {
        std::remove((out_dir + "dm_nao.txt").c_str());
    }
}


#ifdef __MPI
int main(int argc, char** argv)
{
    GlobalV::MY_RANK = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);

    testing::InitGoogleTest(&argc, argv);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

    if (GlobalV::MY_RANK != 0) {
        delete listeners.Release(listeners.default_result_printer());
    }

    int result = RUN_ALL_TESTS();
    MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (GlobalV::MY_RANK == 0 && result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
    }

    MPI_Finalize();
    return result;
}
#endif
