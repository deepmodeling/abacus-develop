#include "gmock/gmock.h"
#include "gtest/gtest.h"

#define private public
#include "source_io/module_parameter/parameter.h"
#undef private

#include "source_base/global_variable.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_io/module_dm/write_dmr.h"
#include "source_io/module_hs/write_HS_R.h"
#include "source_io/module_hs/write_HS_sparse.h"
#include "source_lcao/module_hcontainer/atom_pair.h"
#include "source_lcao/module_hcontainer/hcontainer.h"

#include <complex>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace sparse_format
{
void cal_dH(const UnitCell&,
            const Parallel_Orbitals&,
            LCAO_HS_Arrays&,
            const Grid_Driver&,
            const TwoCenterBundle&,
            const LCAO_Orbitals&,
            const int&,
            const double&,
            const ModuleBase::matrix&)
{
    FAIL() << "cal_dH should not be called by writer compatibility tests.";
}

void cal_dS(const UnitCell&,
            const Parallel_Orbitals&,
            LCAO_HS_Arrays&,
            const Grid_Driver&,
            const TwoCenterBundle&,
            const LCAO_Orbitals&,
            const double&)
{
    FAIL() << "cal_dS should not be called by writer compatibility tests.";
}

void cal_TR(const UnitCell&,
            const Parallel_Orbitals&,
            LCAO_HS_Arrays&,
            const Grid_Driver&,
            const TwoCenterBundle&,
            const LCAO_Orbitals&,
            const double&)
{
    FAIL() << "cal_TR should not be called by writer compatibility tests.";
}

template <typename TK>
void cal_SR(const Parallel_Orbitals&,
            std::set<Abfs::Vector3_Order<int>>&,
            std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>&,
            std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>&,
            const Grid_Driver&,
            const double&,
            hamilt::Hamilt<TK>*)
{
    FAIL() << "cal_SR should not be called by writer compatibility tests.";
}

void destroy_dH_R_sparse(LCAO_HS_Arrays&) {}
void destroy_HS_R_sparse(LCAO_HS_Arrays&) {}
void destroy_T_R_sparse(LCAO_HS_Arrays&) {}

template void cal_SR<double>(
    const Parallel_Orbitals&,
    std::set<Abfs::Vector3_Order<int>>&,
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>&,
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>&,
    const Grid_Driver&,
    const double&,
    hamilt::Hamilt<double>*);

template void cal_SR<std::complex<double>>(
    const Parallel_Orbitals&,
    std::set<Abfs::Vector3_Order<int>>&,
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>&,
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>&,
    const Grid_Driver&,
    const double&,
    hamilt::Hamilt<std::complex<double>>*);
} // namespace sparse_format

namespace
{
std::string read_file(const std::string& filename)
{
    std::ifstream ifs(filename.c_str());
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return oss.str();
}

std::vector<int> read_binary_ints(const std::string& filename, const size_t count)
{
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    std::vector<int> values(count, 0);
    for (size_t i = 0; i < count; ++i)
    {
        ifs.read(reinterpret_cast<char*>(&values[i]), sizeof(int));
    }
    return values;
}

int count_substr(const std::string& text, const std::string& pattern)
{
    int count = 0;
    std::string::size_type pos = 0;
    while ((pos = text.find(pattern, pos)) != std::string::npos)
    {
        ++count;
        pos += pattern.size();
    }
    return count;
}

void init_unitcell(UnitCell& ucell)
{
    ucell.latName = "user_defined_lattice";
    ucell.lat0 = 10.0;
    ucell.latvec.e11 = 1.0;
    ucell.latvec.e22 = 1.0;
    ucell.latvec.e33 = 1.0;
    ucell.ntype = 1;
    ucell.nat = 1;
    ucell.atoms = new Atom[1];
    ucell.set_atom_flag = true;
    ucell.atoms[0].label = "Si";
    ucell.atoms[0].na = 1;
    ucell.atoms[0].nw = 2;
    ucell.atoms[0].taud.resize(1);
    ucell.atoms[0].taud[0] = ModuleBase::Vector3<double>(0.0, 0.25, 0.5);
}

void init_serial_orbitals(Parallel_Orbitals& pv)
{
    pv.atom_begin_row.resize(2);
    pv.atom_begin_col.resize(2);
    pv.atom_begin_row[0] = 0;
    pv.atom_begin_row[1] = 2;
    pv.atom_begin_col[0] = 0;
    pv.atom_begin_col[1] = 2;
    pv.nrow = 2;
    pv.ncol = 2;
    pv.set_serial(2, 2);
}

void fill_matrix(hamilt::HContainer<double>& matrix, Parallel_Orbitals& pv, double* values)
{
    hamilt::AtomPair<double> pair(0, 0, 0, 0, 0, &pv, values);
    matrix.insert_pair(pair);
}

bool starts_with(const std::string& text, const std::string& prefix)
{
    return text.find(prefix) == 0;
}
} // namespace

TEST(WriteHsRCompatibility, FileNameHelpersKeepCurrentContract)
{
    EXPECT_EQ(ModuleIO::hsr_gen_fname("hrs", 0, true, -1), "hrs1_nao.csr");
    EXPECT_EQ(ModuleIO::hsr_gen_fname("hrs", 1, true, 0), "hrs2_nao.csr");
    EXPECT_EQ(ModuleIO::hsr_gen_fname("hrs", 1, false, 3), "hrs2g4_nao.csr");
    EXPECT_EQ(ModuleIO::hsr_gen_fname("srs", 0, false, 0), "srs1g1_nao.csr");
    EXPECT_EQ(ModuleIO::hsr_gen_fname("srs", 0, false, -1), "srs1_nao.csr");

    EXPECT_EQ(ModuleIO::dhr_gen_fname("dhrx", 0, true, -1), "dhrxrs1_nao.csr");
    EXPECT_EQ(ModuleIO::dhr_gen_fname("dhrx", 0, false, 0), "dhrxrs1g1_nao.csr");
    EXPECT_EQ(ModuleIO::dhr_gen_fname("dsry", 1, false, 2), "dsryrs2g3_nao.csr");

    EXPECT_EQ(ModuleIO::dmr_gen_fname(1, 0, true, -1), "dmrs1_nao.csr");
    EXPECT_EQ(ModuleIO::dmr_gen_fname(1, 1, false, 2), "dmrs2g3_nao.csr");
    EXPECT_EQ(ModuleIO::dmr_gen_fname(2, 0, true, 5), "dmrs1_nao.npz");
}

TEST(WriteHsRCompatibility, HContainerCsrHeaderKeepsCurrentFormat)
{
    const std::string filename = "write_hs_r_header_h.csr";
    std::remove(filename.c_str());

    UnitCell ucell;
    init_unitcell(ucell);
    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    hamilt::HContainer<double> matrix(&pv);
    double values[4] = {1.0, 0.0, 0.5, 2.0};
    fill_matrix(matrix, pv, values);

    ModuleIO::write_hcontainer_csr(filename, &ucell, 5, &matrix, 0, 0, 1, "H");

    const std::string output = read_file(filename);
    EXPECT_THAT(output, testing::HasSubstr(" --- Ionic Step 1 ---\n"));
    EXPECT_THAT(output, testing::HasSubstr(" # print H matrix in real space H(R)\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1 # number of spin directions\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1 # spin index\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 2 # number of localized basis\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1 # number of Bravais lattice vector R\n"));
    EXPECT_THAT(output, testing::HasSubstr(" user_defined_lattice\n"));
    EXPECT_THAT(output, testing::HasSubstr(" Si\n"));
    EXPECT_THAT(output, testing::HasSubstr(" Direct\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 0 0.25 0.5\n"));
    EXPECT_THAT(output, testing::HasSubstr(" #                               CSR Format                             #\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 0 0 0 3\n"));
    EXPECT_THAT(output, testing::HasSubstr(" # CSR values\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1.00000e+00 5.00000e-01 2.00000e+00"));

    std::remove(filename.c_str());
}

TEST(WriteHsRCompatibility, HContainerCsrAppendKeepsCurrentStepSections)
{
    const std::string filename = "write_hs_r_append_s.csr";
    std::remove(filename.c_str());

    UnitCell ucell;
    init_unitcell(ucell);
    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    hamilt::HContainer<double> matrix(&pv);
    double values[4] = {1.0, 0.0, 0.0, 1.0};
    fill_matrix(matrix, pv, values);

    ModuleIO::write_hcontainer_csr(filename, &ucell, 4, &matrix, 0, 0, 1, "S");
    ModuleIO::write_hcontainer_csr(filename, &ucell, 4, &matrix, 1, 0, 1, "S");

    const std::string output = read_file(filename);
    EXPECT_EQ(count_substr(output, " --- Ionic Step "), 2);
    EXPECT_THAT(output, testing::HasSubstr(" --- Ionic Step 1 ---\n"));
    EXPECT_THAT(output, testing::HasSubstr(" --- Ionic Step 2 ---\n"));
    EXPECT_EQ(count_substr(output, " # print S matrix in real space S(R)\n"), 2);

    std::remove(filename.c_str());
}

TEST(WriteHsRCompatibility, DmrCsrHeaderKeepsCurrentFormat)
{
    std::string filename = "write_hs_r_header_dmr.csr";
    std::remove(filename.c_str());

    UnitCell ucell;
    init_unitcell(ucell);
    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    hamilt::HContainer<double> matrix(&pv);
    double values[4] = {0.25, 0.0, 0.0, 0.75};
    fill_matrix(matrix, pv, values);

    ModuleIO::write_dmr_csr(filename, &ucell, 4, &matrix, 0, 1, 2);

    const std::string output = read_file(filename);
    EXPECT_THAT(output, testing::HasSubstr(" --- Ionic Step 1 ---\n"));
    EXPECT_THAT(output, testing::HasSubstr(" # print density matrix in real space DM(R)\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 2 # number of spin directions\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 2 # spin index\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 2 # number of localized basis\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1 # number of Bravais lattice vector R\n"));

    std::remove(filename.c_str());
}

TEST(WriteHsRCompatibility, LegacySparseHeaderKeepsStepStyle)
{
    const std::string filename = "write_hs_r_legacy_s.csr";
    std::remove(filename.c_str());

    GlobalV::DRANK = 0;
    PARAM.sys.global_out_dir = "./";
    PARAM.sys.nlocal = 99;

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    std::set<Abfs::Vector3_Order<int>> all_R_coor;
    all_R_coor.insert(r_vector);
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> sparse_matrix;
    sparse_matrix[r_vector][0][1] = 0.5;
    sparse_matrix[r_vector][1][1] = 1.5;

    ModuleIO::SparseWriteOptions options;
    options.filename = filename;
    options.label = "S";
    options.threshold = 1e-10;
    options.binary = false;
    options.istep = 0;
    options.reduce = false;
    options.temp_dir = "./";
    ModuleIO::save_sparse(sparse_matrix, all_R_coor, pv, options);

    const std::string output = read_file(filename);
    EXPECT_TRUE(starts_with(output, "STEP: 0\n"));
    EXPECT_THAT(output, testing::HasSubstr("Matrix Dimension of S(R): 2\n"));
    EXPECT_THAT(output, testing::HasSubstr("Matrix number of S(R): 1\n"));
    EXPECT_THAT(output, testing::HasSubstr("0 0 0 2\n"));

    std::remove(filename.c_str());
}

TEST(WriteHsRCompatibility, LegacySparseBinaryHeaderWritesConcreteStep)
{
    const std::string filename = "write_hs_r_legacy_binary_s.csr";
    std::remove(filename.c_str());

    GlobalV::DRANK = 0;
    PARAM.sys.global_out_dir = "./";
    PARAM.sys.nlocal = 99;

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    std::set<Abfs::Vector3_Order<int>> all_R_coor;
    all_R_coor.insert(r_vector);
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> sparse_matrix;
    sparse_matrix[r_vector][0][1] = 0.5;
    sparse_matrix[r_vector][1][1] = 1.5;

    ModuleIO::SparseWriteOptions options;
    options.filename = filename;
    options.label = "S";
    options.threshold = 1e-10;
    options.binary = true;
    options.istep = 3;
    options.reduce = false;
    options.temp_dir = "./";
    ModuleIO::save_sparse(sparse_matrix, all_R_coor, pv, options);

    const std::vector<int> header_and_r = read_binary_ints(filename, 7);
    EXPECT_THAT(header_and_r, testing::ElementsAre(3, 2, 1, 0, 0, 0, 2));

    std::remove(filename.c_str());
}

TEST(WriteHsRCompatibility, HeaderStyleSamplesRemainDistinct)
{
    EXPECT_TRUE(starts_with(" --- Ionic Step 1 ---", " --- Ionic Step"));
    EXPECT_TRUE(starts_with("STEP: 0", "STEP:"));
    EXPECT_TRUE(starts_with("IONIC_STEP: 1", "IONIC_STEP:"));
}
