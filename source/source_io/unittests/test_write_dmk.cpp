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

TEST(DMKTest, ReadDMKRoundTripMultipleMatrices) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

    const int nspin = 2;
    const int nk = 2;
    const int nlocal = 2;
    Parallel_Blkcyclic pv;
    init_pv(nlocal, pv);

    K_Vectors kv;
    init_kv(kv, 0.5);
    kv.set_nkstot(nspin * nk);
    kv.set_nkstot_nospin(nk);
    kv.set_nks(nspin * nk);
    kv.kvec_c.resize(nk);
    kv.kvec_d.resize(nk);
    kv.wk.resize(nk);
    kv.isk.resize(nk);
    for (int ik = 0; ik < nk; ++ik)
    {
        kv.kvec_c[ik].x = 0.1 * ik;
        kv.kvec_d[ik].x = 0.1 * ik;
        kv.wk[ik] = 0.5;
        kv.isk[ik] = 0;
    }

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
    ASSERT_EQ(dmk_read.size(), dmk.size());
    for (int index = 0; index < nspin * nk; ++index)
    {
        ASSERT_EQ(dmk_read[index].size(), dmk[index].size());
        for (std::size_t element = 0; element < dmk[index].size(); ++element)
        {
            EXPECT_NEAR(dmk_read[index][element].real(), dmk[index][element].real(), 1e-8);
            EXPECT_NEAR(dmk_read[index][element].imag(), dmk[index][element].imag(), 1e-8);
        }
    }

    if (GlobalV::MY_RANK == 0)
    {
        remove((out_dir + "dmk1s1_nao.txt").c_str());
        remove((out_dir + "dmk2s1_nao.txt").c_str());
        remove((out_dir + "dmk1s2_nao.txt").c_str());
        remove((out_dir + "dmk2s2_nao.txt").c_str());
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

    std::vector<std::vector<double>> dmk(1, std::vector<double>(nlocal * nlocal, 7.0));
    const std::vector<std::vector<double>> original = dmk;
    std::ofstream ofs_log;
    const std::string missing_dir = testing::TempDir() + "dmk_no_such_dir/";
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, missing_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk, ofs_log));
    EXPECT_EQ(dmk, original);
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

// ---------------------------------------------------------------------------
// A failure on a later (spin, k) point must discard the matrices that were
// already read; otherwise the caller would keep a half-updated DM(k) and run
// the next SCF from a mixture of two restarts.
// ---------------------------------------------------------------------------
TEST(DMKTest, ReadDMKPartialFailureKeepsCallerData) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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
    if (GlobalV::MY_RANK == 0)
    {
        std::remove((out_dir + "dmk2_nao.txt").c_str());
    }

    const std::complex<double> sentinel(7.0, 7.0);
    std::vector<std::vector<std::complex<double>>> dmk_read(
        nspin * nk, std::vector<std::complex<double>>(pv.nrloc() * pv.ncloc(), sentinel));

    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));
    EXPECT_TRUE(all_equal_to(dmk_read, sentinel));

    if (GlobalV::MY_RANK == 0)
    {
        std::remove((out_dir + "dmk1_nao.txt").c_str());
    }
    delete ucell;
}

// ---------------------------------------------------------------------------
// A spin-down file whose header claims spin 1 must be rejected instead of being
// loaded as the down channel.
// ---------------------------------------------------------------------------
TEST(DMKTest, ReadDMKRejectsSpinMismatch) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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
    if (GlobalV::MY_RANK == 0)
    {
        EXPECT_TRUE(copy_file(out_dir + "dmk1s1_nao.txt", out_dir + "dmk1s2_nao.txt"));
    }

    const std::complex<double> sentinel(7.0, 7.0);
    std::vector<std::vector<std::complex<double>>> dmk_read(
        nspin * nk, std::vector<std::complex<double>>(pv.nrloc() * pv.ncloc(), sentinel));
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));
    EXPECT_TRUE(all_equal_to(dmk_read, sentinel));

    if (GlobalV::MY_RANK == 0)
    {
        std::remove((out_dir + "dmk1s1_nao.txt").c_str());
        std::remove((out_dir + "dmk1s2_nao.txt").c_str());
    }
    delete ucell;
}

// ---------------------------------------------------------------------------
// k-point metadata: a different k-mesh and a file holding the wrong k-point
// must both be rejected instead of being used for the requested index.
// ---------------------------------------------------------------------------
TEST(DMKTest, ReadDMKRejectsKPointMismatch) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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
    if (GlobalV::MY_RANK == 0)
    {
        EXPECT_TRUE(copy_file(out_dir + "dmk2_nao.txt", out_dir + "dmk1_nao.txt"));
    }
    dmk_read.clear();
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    if (GlobalV::MY_RANK == 0)
    {
        std::remove((out_dir + "dmk1_nao.txt").c_str());
        std::remove((out_dir + "dmk2_nao.txt").c_str());
    }
    delete ucell;
}

// ---------------------------------------------------------------------------
// k-point coordinates: a mesh of the same size but with reordered or shifted
// k-vectors must be rejected, because cal_dmr pairs the read DM(k) with the
// current KPT phases and would silently build a wrong density.
// ---------------------------------------------------------------------------
TEST(DMKTest, ReadDMKRejectsKVecMismatch) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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

    if (GlobalV::MY_RANK == 0)
    {
        std::remove((out_dir + "dmk1_nao.txt").c_str());
        std::remove((out_dir + "dmk2_nao.txt").c_str());
    }
    delete ucell;
}

// ---------------------------------------------------------------------------
// A basis size that differs from the current calculation must be rejected, not
// silently truncated or over-read.
// ---------------------------------------------------------------------------
TEST(DMKTest, ReadDMKRejectsBasisMismatch) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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

    if (GlobalV::MY_RANK == 0)
    {
        std::remove((out_dir + "dmk1_nao.txt").c_str());
    }
    delete ucell;
}

// ---------------------------------------------------------------------------
// Truncated input: a file cut inside the header and a complete file whose last
// body line was lost must both fail instead of yielding a partial matrix.
// ---------------------------------------------------------------------------
TEST(DMKTest, ReadDMKRejectsTruncatedFile) {
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
    if (GlobalV::MY_RANK == 0)
    {
        std::ofstream ofs((out_dir + "dm_nao.txt").c_str());
        ofs << "1\n";
        ofs.close();
    }
    std::vector<std::vector<double>> dmk_read;
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    // A complete record whose final line was lost: the body is short.
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

    const std::vector<double> probe = {1.0, 2.0, 3.0, 4.0};
    std::vector<std::vector<double>> dmk;
    fill_local_dmk(probe, nlocal, pv, dmk);
    const std::vector<double> efs(nspin, 0.0);
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    if (GlobalV::MY_RANK == 0)
    {
        EXPECT_TRUE(drop_last_line(out_dir + "dm_nao.txt"));
    }
    dmk_read.clear();
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    if (GlobalV::MY_RANK == 0)
    {
        std::remove((out_dir + "dm_nao.txt").c_str());
    }
    delete ucell;
}

// ---------------------------------------------------------------------------
// The flat index shared by the readers and the writers: a row-major buffer
// stores (i, j) at i * nlocal + j, a column-major one at j * nlocal + i. All
// LCAO ks_solver choices are column-major, so this mapping is what keeps the
// file row-major for every supported solver.
// ---------------------------------------------------------------------------
TEST(DMKTest, DenseBufferIndexSelectsLayout) {
    EXPECT_EQ(ModuleIO::dense_buffer_index(0, 0, 3, /*column_major=*/false), 0u);
    EXPECT_EQ(ModuleIO::dense_buffer_index(0, 1, 3, /*column_major=*/false), 1u);
    EXPECT_EQ(ModuleIO::dense_buffer_index(1, 0, 3, /*column_major=*/false), 3u);
    EXPECT_EQ(ModuleIO::dense_buffer_index(2, 1, 3, /*column_major=*/false), 7u);
    EXPECT_EQ(ModuleIO::dense_buffer_index(0, 1, 3, /*column_major=*/true), 3u);
    EXPECT_EQ(ModuleIO::dense_buffer_index(1, 0, 3, /*column_major=*/true), 1u);
    EXPECT_EQ(ModuleIO::dense_buffer_index(1, 2, 3, /*column_major=*/true), 7u);
}

// ---------------------------------------------------------------------------
// The header must state the conj-first element convention, and read_dmk must
// transparently skip that comment line when the file is read back.
// ---------------------------------------------------------------------------
TEST(DMKTest, HeaderStatesConjFirstConvention) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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

    if (GlobalV::MY_RANK == 0)
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

// ---------------------------------------------------------------------------
// With column_major = false the reader and the writer must still agree on the
// same layout, otherwise a restart would read a transposed matrix. No current
// LCAO solver takes this branch, so this only pins their mutual consistency.
// ---------------------------------------------------------------------------
TEST(DMKTest, ReadDMKRoundTripRowMajorFlag) {
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
    ASSERT_EQ(dmk_read.size(), 1u);
    ASSERT_EQ(dmk_read[0].size(), dmk[0].size());
    for (std::size_t idx = 0; idx < dmk[0].size(); ++idx)
    {
        EXPECT_NEAR(dmk_read[0][idx].real(), dmk[0][idx].real(), 1e-8);
        EXPECT_NEAR(dmk_read[0][idx].imag(), dmk[0][idx].imag(), 1e-8);
    }

    if (GlobalV::MY_RANK == 0)
    {
        std::remove((out_dir + "dmk1_nao.txt").c_str());
    }
    delete ucell;
}

// ---------------------------------------------------------------------------
// The reader must reject non-finite matrix components: strtod reports an
// overflowing literal such as "1e999" only through errno while consuming the
// full token, and both strtod and the stream extractor accept the literal
// spellings "inf"/"nan" without flagging a failure. Such a component would
// otherwise be installed as the restart density.
// ---------------------------------------------------------------------------

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

TEST(DMKTest, ReadDMKRejectsOverflowingComplexComponent) {
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

    const std::vector<double> efs(nspin, 0.0);
    K_Vectors kv;
    init_kv(kv, 1.0);

    const std::string out_dir = testing::TempDir();
    std::ofstream ofs_log;
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);

    // strtod consumes "1e999" fully and only reports the overflow via ERANGE.
    if (GlobalV::MY_RANK == 0)
    {
        EXPECT_TRUE(patch_first_body_token(out_dir + "dmk1_nao.txt", "(1e999,0)", true));
    }
    std::vector<std::vector<std::complex<double>>> dmk_read;
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    // The literal "nan" parses without any errno flag and must be rejected too.
    ModuleIO::write_dmk(dmk, kv, /*precision=*/12, efs, ucell, pv, /*column_major=*/true,
                        out_dir, /*istep=*/0, /*out_app_flag=*/0, ofs_log);
    if (GlobalV::MY_RANK == 0)
    {
        EXPECT_TRUE(patch_first_body_token(out_dir + "dmk1_nao.txt", "(nan,0)", true));
    }
    dmk_read.clear();
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    if (GlobalV::MY_RANK == 0)
    {
        std::remove((out_dir + "dmk1_nao.txt").c_str());
    }
    delete ucell;
}

TEST(DMKTest, ReadDMKRejectsNonFiniteRealComponent) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

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
    if (GlobalV::MY_RANK == 0)
    {
        EXPECT_TRUE(patch_first_body_token(out_dir + "dm_nao.txt", "inf", false));
    }
    std::vector<std::vector<double>> dmk_read;
    EXPECT_FALSE(ModuleIO::read_dmk(nspin, nk, kv, pv, /*column_major=*/true, out_dir,
                                    /*istep=*/0, /*out_app_flag=*/0, dmk_read, ofs_log));

    if (GlobalV::MY_RANK == 0)
    {
        std::remove((out_dir + "dm_nao.txt").c_str());
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
