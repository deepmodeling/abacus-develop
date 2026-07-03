#include "gmock/gmock.h"
#include "gtest/gtest.h"

#define private public
#include "source_io/module_parameter/parameter.h"
#undef private

#include "source_base/global_variable.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_io/module_dm/write_dmr.h"
#include "source_io/module_hs/output_mat_sparse.h"
#include "source_io/module_hs/rr_sparse_writer.h"
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

#ifdef __MPI
#include <mpi.h>
#endif

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

template <typename T>
T read_binary_value(std::ifstream& ifs)
{
    T value{};
    ifs.read(reinterpret_cast<char*>(&value), sizeof(T));
    return value;
}

std::vector<std::string> read_lines(const std::string& filename)
{
    std::ifstream ifs(filename.c_str());
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(ifs, line))
    {
        lines.push_back(line);
    }
    return lines;
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

void init_sparse_output_globals(const int nspin = 1)
{
    GlobalV::DRANK = 0;
    PARAM.input.nspin = nspin;
    PARAM.input.calculation = "scf";
    PARAM.input.out_app_flag = false;
    PARAM.sys.global_out_dir = "./";
    PARAM.sys.global_matrix_dir = "./";
    PARAM.sys.nlocal = 2;
}

void remove_derivative_files(const std::string& fileflag)
{
    const std::vector<std::string> filenames = {
        "d" + fileflag + "rxs1_nao.csr",
        "d" + fileflag + "rys1_nao.csr",
        "d" + fileflag + "rzs1_nao.csr",
        "d" + fileflag + "rxs2_nao.csr",
        "d" + fileflag + "rys2_nao.csr",
        "d" + fileflag + "rzs2_nao.csr",
    };
    for (const std::string& filename: filenames)
    {
        std::remove(filename.c_str());
    }
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

TEST(WriteHsRCompatibility, LegacySparseTextCountsOnlyValuesAboveThreshold)
{
    const std::string filename = "write_hs_r_threshold_s.csr";
    std::remove(filename.c_str());

    GlobalV::DRANK = 0;
    PARAM.sys.global_out_dir = "./";

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    std::set<Abfs::Vector3_Order<int>> all_R_coor;
    all_R_coor.insert(r_vector);
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> sparse_matrix;
    sparse_matrix[r_vector][0][0] = 1.0;
    sparse_matrix[r_vector][0][1] = 1e-12;
    sparse_matrix[r_vector][1][0] = 0.0;
    sparse_matrix[r_vector][1][1] = -2.0;

    ModuleIO::SparseWriteOptions options;
    options.filename = filename;
    options.label = "S";
    options.threshold = 1e-10;
    options.binary = false;
    options.istep = 2;
    options.reduce = false;
    options.temp_dir = "./";
    ModuleIO::save_sparse(sparse_matrix, all_R_coor, pv, options);

    const std::vector<std::string> lines = read_lines(filename);
    ASSERT_GE(lines.size(), 6);
    EXPECT_EQ(lines[0], "STEP: 2");
    EXPECT_EQ(lines[1], "Matrix Dimension of S(R): 2");
    EXPECT_EQ(lines[2], "Matrix number of S(R): 1");
    EXPECT_EQ(lines[3], "0 0 0 2");

    std::istringstream value_stream(lines[4]);
    std::vector<double> values;
    double value = 0.0;
    while (value_stream >> value)
    {
        values.push_back(value);
    }
    EXPECT_THAT(values, testing::ElementsAre(1.0, -2.0));

    std::istringstream column_stream(lines[5]);
    std::vector<int> columns;
    int column = 0;
    while (column_stream >> column)
    {
        columns.push_back(column);
    }
    EXPECT_THAT(columns, testing::ElementsAre(0, 1));

    ASSERT_GE(lines.size(), 7);
    std::istringstream indptr_stream(lines[6]);
    std::vector<long long> indptr;
    long long ptr = 0;
    while (indptr_stream >> ptr)
    {
        indptr.push_back(ptr);
    }
    EXPECT_THAT(indptr, testing::ElementsAre(0, 1, 2));

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

TEST(WriteHsRCompatibility, LegacySparseBinaryCountsOnlyValuesAboveThreshold)
{
    const std::string filename = "write_hs_r_threshold_binary_s.csr";
    std::remove(filename.c_str());

    GlobalV::DRANK = 0;
    PARAM.sys.global_out_dir = "./";

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    std::set<Abfs::Vector3_Order<int>> all_R_coor;
    all_R_coor.insert(r_vector);
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> sparse_matrix;
    sparse_matrix[r_vector][0][0] = 1.0;
    sparse_matrix[r_vector][0][1] = 1e-12;
    sparse_matrix[r_vector][1][0] = 0.0;
    sparse_matrix[r_vector][1][1] = -2.0;

    ModuleIO::SparseWriteOptions options;
    options.filename = filename;
    options.label = "S";
    options.threshold = 1e-10;
    options.binary = true;
    options.istep = 4;
    options.reduce = false;
    options.temp_dir = "./";
    ModuleIO::save_sparse(sparse_matrix, all_R_coor, pv, options);

    std::ifstream ifs(filename.c_str(), std::ios::binary);
    ASSERT_TRUE(ifs.is_open());
    EXPECT_EQ(read_binary_value<int>(ifs), 4);
    EXPECT_EQ(read_binary_value<int>(ifs), 2);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 2);
    EXPECT_DOUBLE_EQ(read_binary_value<double>(ifs), 1.0);
    EXPECT_DOUBLE_EQ(read_binary_value<double>(ifs), -2.0);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_EQ(read_binary_value<long long>(ifs), 0);
    EXPECT_EQ(read_binary_value<long long>(ifs), 1);
    EXPECT_EQ(read_binary_value<long long>(ifs), 2);

    std::remove(filename.c_str());
}

TEST(WriteHsRCompatibility, SaveDHSparseTextCountsOnlyValuesAboveThreshold)
{
    remove_derivative_files("h");
    init_sparse_output_globals();

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    LCAO_HS_Arrays arrays;
    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    arrays.all_R_coor.insert(r_vector);
    arrays.dHRx_sparse[0][r_vector][0][0] = 1.0;
    arrays.dHRx_sparse[0][r_vector][0][1] = 1e-12;
    arrays.dHRx_sparse[0][r_vector][1][0] = 0.0;
    arrays.dHRx_sparse[0][r_vector][1][1] = -2.0;

    ModuleIO::save_dH_sparse(5, pv, arrays, 1e-10, false, "h", 8);

    const std::vector<std::string> lines = read_lines("dhrxs1_nao.csr");
    ASSERT_GE(lines.size(), 7);
    EXPECT_EQ(lines[0], "STEP: 5");
    EXPECT_EQ(lines[1], "Matrix Dimension of dHx(R): 2");
    EXPECT_EQ(lines[2], "Matrix number of dHx(R): 1");
    EXPECT_EQ(lines[3], "0 0 0 2");
    EXPECT_THAT(lines[4], testing::HasSubstr("1.00000000e+00"));
    EXPECT_THAT(lines[4], testing::HasSubstr("-2.00000000e+00"));

    std::istringstream column_stream(lines[5]);
    std::vector<int> columns;
    int column = 0;
    while (column_stream >> column)
    {
        columns.push_back(column);
    }
    EXPECT_THAT(columns, testing::ElementsAre(0, 1));

    std::istringstream indptr_stream(lines[6]);
    std::vector<long long> indptr;
    long long ptr = 0;
    while (indptr_stream >> ptr)
    {
        indptr.push_back(ptr);
    }
    EXPECT_THAT(indptr, testing::ElementsAre(0, 1, 2));

    const std::vector<std::string> y_lines = read_lines("dhrys1_nao.csr");
    ASSERT_GE(y_lines.size(), 4);
    EXPECT_EQ(y_lines[3], "0 0 0 0");

    remove_derivative_files("h");
}

TEST(WriteHsRCompatibility, SaveDHSparseBinaryCountsOnlyValuesAboveThreshold)
{
    remove_derivative_files("h");
    init_sparse_output_globals();

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    LCAO_HS_Arrays arrays;
    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    arrays.all_R_coor.insert(r_vector);
    arrays.dHRx_sparse[0][r_vector][0][0] = 1.0;
    arrays.dHRx_sparse[0][r_vector][0][1] = 1e-12;
    arrays.dHRx_sparse[0][r_vector][1][0] = 0.0;
    arrays.dHRx_sparse[0][r_vector][1][1] = -2.0;

    ModuleIO::save_dH_sparse(6, pv, arrays, 1e-10, true, "h", 8);

    std::ifstream ifs("dhrxs1_nao.csr", std::ios::binary);
    ASSERT_TRUE(ifs.is_open());
    EXPECT_EQ(read_binary_value<int>(ifs), 6);
    EXPECT_EQ(read_binary_value<int>(ifs), 2);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 2);
    EXPECT_DOUBLE_EQ(read_binary_value<double>(ifs), 1.0);
    EXPECT_DOUBLE_EQ(read_binary_value<double>(ifs), -2.0);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_EQ(read_binary_value<long long>(ifs), 0);
    EXPECT_EQ(read_binary_value<long long>(ifs), 1);
    EXPECT_EQ(read_binary_value<long long>(ifs), 2);

    remove_derivative_files("h");
}

TEST(WriteHsRCompatibility, SaveDSSparseSocWritesAllDirections)
{
    remove_derivative_files("s");
    init_sparse_output_globals(4);

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    LCAO_HS_Arrays arrays;
    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    arrays.all_R_coor.insert(r_vector);
    arrays.dHRx_soc_sparse[r_vector][0][0] = std::complex<double>(1.0, 0.0);
    arrays.dHRy_soc_sparse[r_vector][0][1] = std::complex<double>(2.0, -1.0);
    arrays.dHRz_soc_sparse[r_vector][1][1] = std::complex<double>(-3.0, 0.5);

    ModuleIO::save_dH_sparse(7, pv, arrays, 1e-10, false, "s", 8);

    const std::string x_output = read_file("dsrxs1_nao.csr");
    const std::string y_output = read_file("dsrys1_nao.csr");
    const std::string z_output = read_file("dsrzs1_nao.csr");
    EXPECT_THAT(x_output, testing::HasSubstr("Matrix number of dHx(R): 1\n0 0 0 1\n"));
    EXPECT_THAT(y_output, testing::HasSubstr("Matrix number of dHy(R): 1\n0 0 0 1\n"));
    EXPECT_THAT(z_output, testing::HasSubstr("Matrix number of dHz(R): 1\n0 0 0 1\n"));
    EXPECT_THAT(x_output, testing::HasSubstr("(1.00000000e+00,0.00000000e+00)"));
    EXPECT_THAT(y_output, testing::HasSubstr("(2.00000000e+00,-1.00000000e+00)"));
    EXPECT_THAT(z_output, testing::HasSubstr("(-3.00000000e+00,5.00000000e-01)"));

    remove_derivative_files("s");
}

TEST(WriteHsRCompatibility, MatSparseOutputOptionsKeepLegacyDefaults)
{
    ModuleIO::MatSparseOutputOptions options;
    EXPECT_FALSE(options.out_mat_dh);
    EXPECT_FALSE(options.out_mat_ds);
    EXPECT_FALSE(options.out_mat_t);
    EXPECT_FALSE(options.out_mat_r);
    EXPECT_EQ(options.dh_precision, 16);
    EXPECT_EQ(options.ds_precision, 16);
    EXPECT_EQ(options.t_precision, 16);
    EXPECT_EQ(options.r_precision, 16);
    EXPECT_DOUBLE_EQ(options.sparse_threshold, 1e-10);
    EXPECT_FALSE(options.binary);
}

TEST(WriteHsRCompatibility, RRSparsePayloadDetectorSkipsEmptyBlocks)
{
    int empty_counts[3] = {0, 0, 0};
    int x_only_counts[3] = {1, 0, 0};
    int z_only_counts[3] = {0, 0, 2};

    EXPECT_FALSE(ModuleIO::detail::rr_sparse_has_payload(empty_counts));
    EXPECT_TRUE(ModuleIO::detail::rr_sparse_has_payload(x_only_counts));
    EXPECT_TRUE(ModuleIO::detail::rr_sparse_has_payload(z_only_counts));
}

TEST(WriteHsRCompatibility, RRSparseTextFinalizerAllowsZeroBlocks)
{
    const std::string payload_filename = "rr_empty_payload.tmp";
    const std::string output_filename = "rr_empty.csr";
    std::remove(payload_filename.c_str());
    std::remove(output_filename.c_str());

    std::ofstream payload(payload_filename.c_str());
    payload.close();

    ModuleIO::detail::finalize_rr_sparse_file(output_filename,
                                              payload_filename,
                                              9,
                                              2,
                                              0,
                                              false,
                                              false,
                                              "WriteHsRCompatibility");

    const std::vector<std::string> lines = read_lines(output_filename);
    ASSERT_EQ(lines.size(), 3);
    EXPECT_EQ(lines[0], "STEP: 9");
    EXPECT_EQ(lines[1], "Matrix Dimension of r(R): 2");
    EXPECT_EQ(lines[2], "Matrix number of r(R): 0");

    std::remove(payload_filename.c_str());
    std::remove(output_filename.c_str());
}

TEST(WriteHsRCompatibility, RRSparseTextFinalizerKeepsSingleDirectionPayload)
{
    const std::string payload_filename = "rr_single_direction_payload.tmp";
    const std::string output_filename = "rr_single_direction.csr";
    std::remove(payload_filename.c_str());
    std::remove(output_filename.c_str());

    std::ofstream payload(payload_filename.c_str());
    payload << "1 0 -1\n";
    payload << "1\n";
    payload << " 4.00000000e+00\n";
    payload << " 0\n";
    payload << "0 1\n";
    payload << "0\n";
    payload << "0\n";
    payload.close();

    ModuleIO::detail::finalize_rr_sparse_file(output_filename,
                                              payload_filename,
                                              10,
                                              2,
                                              1,
                                              false,
                                              false,
                                              "WriteHsRCompatibility");

    const std::vector<std::string> lines = read_lines(output_filename);
    ASSERT_EQ(lines.size(), 10);
    EXPECT_EQ(lines[0], "STEP: 10");
    EXPECT_EQ(lines[1], "Matrix Dimension of r(R): 2");
    EXPECT_EQ(lines[2], "Matrix number of r(R): 1");
    EXPECT_EQ(lines[3], "1 0 -1");
    EXPECT_EQ(lines[4], "1");
    EXPECT_EQ(lines[5], " 4.00000000e+00");
    EXPECT_EQ(lines[8], "0");
    EXPECT_EQ(lines[9], "0");

    std::remove(payload_filename.c_str());
    std::remove(output_filename.c_str());
}

TEST(WriteHsRCompatibility, RRSparseBinaryFinalizerKeepsHeaderAndPayloadOrder)
{
    const std::string payload_filename = "rr_binary_payload.tmp";
    const std::string output_filename = "rr_binary.csr";
    std::remove(payload_filename.c_str());
    std::remove(output_filename.c_str());

    std::ofstream payload(payload_filename.c_str(), std::ios::binary);
    int dRx = 1;
    int dRy = 2;
    int dRz = 3;
    int x_count = 1;
    int y_count = 0;
    int z_count = 0;
    double value = 4.0;
    int column = 1;
    long long ptr0 = 0;
    long long ptr1 = 1;
    payload.write(reinterpret_cast<const char*>(&dRx), sizeof(int));
    payload.write(reinterpret_cast<const char*>(&dRy), sizeof(int));
    payload.write(reinterpret_cast<const char*>(&dRz), sizeof(int));
    payload.write(reinterpret_cast<const char*>(&x_count), sizeof(int));
    payload.write(reinterpret_cast<const char*>(&value), sizeof(double));
    payload.write(reinterpret_cast<const char*>(&column), sizeof(int));
    payload.write(reinterpret_cast<const char*>(&ptr0), sizeof(long long));
    payload.write(reinterpret_cast<const char*>(&ptr1), sizeof(long long));
    payload.write(reinterpret_cast<const char*>(&y_count), sizeof(int));
    payload.write(reinterpret_cast<const char*>(&z_count), sizeof(int));
    payload.close();

    ModuleIO::detail::finalize_rr_sparse_file(output_filename,
                                              payload_filename,
                                              11,
                                              2,
                                              1,
                                              true,
                                              false,
                                              "WriteHsRCompatibility");

    std::ifstream ifs(output_filename, std::ios::binary);
    ASSERT_TRUE(ifs.is_open());
    EXPECT_EQ(read_binary_value<int>(ifs), 11);
    EXPECT_EQ(read_binary_value<int>(ifs), 2);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_EQ(read_binary_value<int>(ifs), 2);
    EXPECT_EQ(read_binary_value<int>(ifs), 3);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_DOUBLE_EQ(read_binary_value<double>(ifs), 4.0);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_EQ(read_binary_value<long long>(ifs), 0);
    EXPECT_EQ(read_binary_value<long long>(ifs), 1);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    ifs.close();

    std::remove(payload_filename.c_str());
    std::remove(output_filename.c_str());
}

TEST(WriteHsRCompatibility, HeaderStyleSamplesRemainDistinct)
{
    EXPECT_TRUE(starts_with(" --- Ionic Step 1 ---", " --- Ionic Step"));
    EXPECT_TRUE(starts_with("STEP: 0", "STEP:"));
    EXPECT_TRUE(starts_with("IONIC_STEP: 1", "IONIC_STEP:"));
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
#endif

    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
