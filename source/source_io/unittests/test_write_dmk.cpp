#include <complex>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "source_io/module_dm/write_dmk.h"
#include "source_base/global_variable.h"
#include "prepare_unitcell.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_base/module_external/scalapack_connector.h"

#ifdef __MPI
#include "mpi.h"
#endif

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

TEST(DMKTest, GenFileName) {
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


TEST(DMKTest,WriteDMK) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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
    
    std::ifstream ifs;

    int pass = 0;
    if (GlobalV::MY_RANK == 0)
    {
        std::string fn = "dms1_nao.txt";
        ifs.open(fn);
        std::string str((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>());
        EXPECT_THAT(str, testing::HasSubstr("0 # Fermi energy in Ry"));
        EXPECT_THAT(str, testing::HasSubstr("20 # number of localized basis"));
        ifs.close();

        fn = "dms2_nao.txt";
        ifs.open(fn);
        str = std::string((std::istreambuf_iterator<char>(ifs)),
                          std::istreambuf_iterator<char>());
        EXPECT_THAT(str, testing::HasSubstr("0.1 # Fermi energy in Ry"));
        EXPECT_THAT(str, testing::HasSubstr("20 # number of localized basis"));
        ifs.close();

        fn = "dmk1s1_nao.txt";
        ifs.open(fn);
        str = std::string((std::istreambuf_iterator<char>(ifs)),
                          std::istreambuf_iterator<char>());
        EXPECT_THAT(str, testing::HasSubstr("0 # Fermi energy in Ry"));
        EXPECT_THAT(str, testing::HasSubstr("20 # number of localized basis"));
        // the matrix body must end with a newline so appended records are not glued
        EXPECT_EQ(str.back(), '\n');
        ifs.close();

        fn = "dmk2s1_nao.txt";
        ifs.open(fn);
        str = std::string((std::istreambuf_iterator<char>(ifs)),
                          std::istreambuf_iterator<char>());
        EXPECT_THAT(str, testing::HasSubstr("0 # Fermi energy in Ry"));
        EXPECT_THAT(str, testing::HasSubstr("20 # number of localized basis"));
        ifs.close();

        fn = "dmk1s2_nao.txt";
        ifs.open(fn);
        str = std::string((std::istreambuf_iterator<char>(ifs)),
                          std::istreambuf_iterator<char>());
        EXPECT_THAT(str, testing::HasSubstr("0.1 # Fermi energy in Ry"));
        EXPECT_THAT(str, testing::HasSubstr("20 # number of localized basis"));
        ifs.close();

        fn = "dmk2s2_nao.txt";
        ifs.open(fn);
        str = std::string((std::istreambuf_iterator<char>(ifs)),
                          std::istreambuf_iterator<char>());
        EXPECT_THAT(str, testing::HasSubstr("0.1 # Fermi energy in Ry"));
        EXPECT_THAT(str, testing::HasSubstr("20 # number of localized basis"));
        ifs.close();
        remove("dms1_nao.txt");
        remove("dms2_nao.txt");
        remove("dmk1s1_nao.txt");
        remove("dmk2s1_nao.txt");
        remove("dmk1s2_nao.txt");
        remove("dmk2s2_nao.txt");
    }

    delete ucell;
    // remove the generated files
    
};

// ---------------------------------------------------------------------------
// The DM(k) text body is row-major even when the local blocks are stored
// column-major, which is what genelpa/elpa/scalapack_gvx and the other
// column-major solvers leave behind. A non-symmetric complex probe makes an
// accidental transpose (the complex conjugate of a Hermitian matrix) visible.
// ---------------------------------------------------------------------------
TEST(DMKTest, BodyIsRowMajorForColumnMajorStorage) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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

    if (GlobalV::MY_RANK == 0)
    {
        // nk = 1 with a multi-k (complex) DM -> dmk1_nao.txt.
        const std::string fn = out_dir + "dmk1_nao.txt";
        const std::vector<double> values = read_numbers(fn);
        ASSERT_EQ(values.size(), 8u);
        EXPECT_NEAR(values[0], 1.0, 1e-8);
        EXPECT_NEAR(values[1], 0.0, 1e-8);
        EXPECT_NEAR(values[2], 2.0, 1e-8);
        EXPECT_NEAR(values[3], 3.0, 1e-8);
        EXPECT_NEAR(values[4], 4.0, 1e-8);
        EXPECT_NEAR(values[5], 5.0, 1e-8);
        EXPECT_NEAR(values[6], 6.0, 1e-8);
        EXPECT_NEAR(values[7], 0.0, 1e-8);
        remove(fn.c_str());
    }

    delete ucell;
}

// ---------------------------------------------------------------------------
// out_app_flag = 0: one fixed filename, truncated and overwritten each call.
// ---------------------------------------------------------------------------
TEST(DMKTest, Mode0FixedNameOverwrite) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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

    if (GlobalV::MY_RANK == 0)
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
    delete ucell;
}

// ---------------------------------------------------------------------------
// out_app_flag = 1: fixed filename, multiple calls append successive records.
// ---------------------------------------------------------------------------
TEST(DMKTest, Mode1FixedNameAppend) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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
    const int out_app_flag = 1;
    std::ofstream ofs_log;

    // Two ionic steps -> both write to dms1_nao.txt.
    ModuleIO::write_dmk(dmk, kv, /*precision=*/3, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, out_app_flag, ofs_log);
    ModuleIO::write_dmk(dmk, kv, /*precision=*/3, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/1, out_app_flag, ofs_log);

    if (GlobalV::MY_RANK == 0)
    {
        std::ifstream ifs("dms1_nao.txt");
        std::string str((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>());
        ifs.close();
        // Two "Ionic_Step" headers present -> two appended records.
        EXPECT_THAT(str, testing::HasSubstr("Ionic_Step 1"));
        EXPECT_THAT(str, testing::HasSubstr("Ionic_Step 2"));
        // g-suffixed file must NOT exist.
        ifs.open("dms1g1_nao.txt");
        EXPECT_FALSE(ifs.good());
        remove("dms1_nao.txt");
    }
    delete ucell;
}

// ---------------------------------------------------------------------------
// out_app_flag = 2: per-step g{step} filename; directory routing is handled by
// set_global_dir and is covered by io_sys_var_test.
// ---------------------------------------------------------------------------
TEST(DMKTest, Mode2PerStepFile) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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

    if (GlobalV::MY_RANK == 0)
    {
        std::ifstream ifs("dms1g3_nao.txt");
        EXPECT_TRUE(ifs.good());
        ifs.close();
        remove("dms1g3_nao.txt");
    }
    delete ucell;
}

// ---------------------------------------------------------------------------
// Round trip: write a column-major DM(k) and read it back with the same
// storage order. A transposed read would return the complex conjugate of the
// Hermitian probe, i.e. the imaginary off-diagonal elements would flip sign.
// ---------------------------------------------------------------------------
TEST(DMKTest, ReadDMKRoundTripComplex) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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
    ASSERT_EQ(dmk_read.size(), 1u);
    ASSERT_EQ(dmk_read[0].size(), dmk[0].size());
    for (std::size_t idx = 0; idx < dmk[0].size(); ++idx)
    {
        EXPECT_NEAR(dmk_read[0][idx].real(), dmk[0][idx].real(), 1e-8);
        EXPECT_NEAR(dmk_read[0][idx].imag(), dmk[0][idx].imag(), 1e-8);
    }

    if (GlobalV::MY_RANK == 0)
    {
        remove((out_dir + "dmk1_nao.txt").c_str());
    }
    delete ucell;
}

TEST(DMKTest, ReadDMKRoundTripGamma) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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
    ASSERT_EQ(dmk_read.size(), 1u);
    ASSERT_EQ(dmk_read[0].size(), dmk[0].size());
    for (std::size_t idx = 0; idx < dmk[0].size(); ++idx)
    {
        EXPECT_NEAR(dmk_read[0][idx], dmk[0][idx], 1e-8);
    }

    if (GlobalV::MY_RANK == 0)
    {
        remove((out_dir + "dm_nao.txt").c_str());
    }
    delete ucell;
}

// A missing file must be reported cooperatively instead of aborting a single
// rank, so the caller can decide how to continue.
TEST(DMKTest, ReadDMKMissingFileReturnsFalse) {
    const int nspin = 1;
    const int nk = 1;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    K_Vectors kv;
    init_kv(kv, 1.0);

    std::vector<std::vector<double>> dmk;
    std::ofstream ofs_log;
    const std::string missing_dir = testing::TempDir() + "dmk_no_such_dir/";
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, missing_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk, ofs_log));
}

// ---------------------------------------------------------------------------
// out_app_flag = 2: the ionic step is part of the filename, so the reader has
// to request the matching g{istep+1} file and reject a frame whose marker
// contradicts the name.
// ---------------------------------------------------------------------------
TEST(DMKTest, ReadDMKMode2PerStepFile) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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
    ASSERT_EQ(dmk_read.size(), 1u);
    ASSERT_EQ(dmk_read[0].size(), dmk[0].size());
    for (std::size_t idx = 0; idx < dmk[0].size(); ++idx)
    {
        EXPECT_NEAR(dmk_read[0][idx].real(), dmk[0][idx].real(), 1e-8);
        EXPECT_NEAR(dmk_read[0][idx].imag(), dmk[0][idx].imag(), 1e-8);
    }

    // A different step maps to another file name, which does not exist.
    std::vector<std::vector<std::complex<double>>> dmk_other_step;
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/2, dmk_other_step, ofs_log));

    if (GlobalV::MY_RANK == 0)
    {
        remove((out_dir + "dmk1g3_nao.txt").c_str());
    }
    delete ucell;
}

// ---------------------------------------------------------------------------
// out_app_flag = 1: a fixed-name file accumulates one record per ionic step,
// so the reader has to select the requested frame instead of always taking the
// first one.
// ---------------------------------------------------------------------------
TEST(DMKTest, ReadDMKMode1FrameSelection) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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
    ASSERT_EQ(dmk_read.size(), 1u);
    for (std::size_t idx = 0; idx < dmk_first[0].size(); ++idx)
    {
        EXPECT_NEAR(dmk_read[0][idx].real(), dmk_first[0][idx].real(), 1e-8);
        EXPECT_NEAR(dmk_read[0][idx].imag(), dmk_first[0][idx].imag(), 1e-8);
    }

    dmk_read.clear();
    EXPECT_TRUE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                   /*istep=*/1, /*out_app_flag=*/1, dmk_read, ofs_log));
    ASSERT_EQ(dmk_read.size(), 1u);
    for (std::size_t idx = 0; idx < dmk_second[0].size(); ++idx)
    {
        EXPECT_NEAR(dmk_read[0][idx].real(), dmk_second[0][idx].real(), 1e-8);
        EXPECT_NEAR(dmk_read[0][idx].imag(), dmk_second[0][idx].imag(), 1e-8);
    }

    // A frame that was never appended must be reported, not silently replaced
    // by the last available record.
    std::vector<std::vector<std::complex<double>>> dmk_missing_frame;
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/5, /*out_app_flag=*/1, dmk_missing_frame, ofs_log));

    if (GlobalV::MY_RANK == 0)
    {
        remove((out_dir + "dmk1_nao.txt").c_str());
    }
    delete ucell;
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
