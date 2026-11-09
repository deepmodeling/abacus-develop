/**
 * @file test_hs_sparse_io.cpp
 * @brief Unit tests for hs_sparse_io.cpp (ModuleIO::save_sparse, both the
 *        compact layout and the csr_format text layout).
 *
 * The tests pin down the legacy CSR text/binary layout written by
 * save_sparse: header style ("Ionic_Step n"), dimension lines, threshold
 * filtering of values, and the native binary field order.
 * hs_sparse_io.cpp only reads GlobalV::DRANK, so the fixture sets it
 * directly; no PARAM state is involved.
 */
#include "csr_test_helpers.h"

#include "source_io/module_hs/hs_sparse_io.h"

#ifdef __MPI
#include <mpi.h>
#endif

TEST(HsSparseIo, SparseHeaderKeepsStepStyle)
{
    const std::string filename = "write_hs_r_legacy_s.csr";
    std::remove(filename.c_str());


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
    ModuleIO::save_sparse(sparse_matrix, all_R_coor, pv, options);

    const std::string output = read_file(filename);
    // istep = 0 is the first ionic step, printed one-based as 1
    EXPECT_TRUE(starts_with(output, "Ionic_Step 1  # ABACUS version:"));
    EXPECT_THAT(output, testing::HasSubstr("Matrix Dimension of S(R): 2\n"));
    EXPECT_THAT(output, testing::HasSubstr("Matrix number of S(R): 1\n"));
    EXPECT_THAT(output, testing::HasSubstr("0 0 0 2\n"));

    std::remove(filename.c_str());
}

TEST(HsSparseIo, SparseTextCountsOnlyValuesAboveThreshold)
{
    const std::string filename = "write_hs_r_threshold_s.csr";
    std::remove(filename.c_str());


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
    ModuleIO::save_sparse(sparse_matrix, all_R_coor, pv, options);

    const std::vector<std::string> lines = read_lines(filename);
    ASSERT_GE(lines.size(), 6);
    // istep = 2 is the third ionic step, printed one-based as 3
    EXPECT_EQ(lines[0].compare(0, 12, "Ionic_Step 3"), 0);
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

TEST(HsSparseIo, SparseBinaryHeaderWritesConcreteStep)
{
    const std::string filename = "write_hs_r_legacy_binary_s.csr";
    std::remove(filename.c_str());


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
    ModuleIO::save_sparse(sparse_matrix, all_R_coor, pv, options);

    // istep = 3 -> one-based step 4 in the binary header
    const std::vector<int> header_and_r = read_binary_ints(filename, 7);
    EXPECT_THAT(header_and_r, testing::ElementsAre(4, 2, 1, 0, 0, 0, 2));

    std::remove(filename.c_str());
}

TEST(HsSparseIo, SparseBinaryCountsOnlyValuesAboveThreshold)
{
    const std::string filename = "write_hs_r_threshold_binary_s.csr";
    std::remove(filename.c_str());


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
    ModuleIO::save_sparse(sparse_matrix, all_R_coor, pv, options);

    std::ifstream ifs(filename.c_str(), std::ios::binary);
    ASSERT_TRUE(ifs.is_open());
    // istep = 4 -> one-based step 5 in the binary header
    EXPECT_EQ(read_binary_value<int>(ifs), 5);
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

// A final-snapshot writer such as write_Vxc_R passes out_app_flag = 0 so a
// rerun truncates the previous output even when the only write happens at
// istep > 0.  This test pins that behaviour: seeding a file with stale
// content and then calling save_sparse with out_app_flag = 0 and istep = 2
// must leave a single frame, not an appended duplicate.
TEST(HsSparseIo, SparseTextOutAppFlagZeroTruncatesSeededFile)
{
    const std::string filename = testing::TempDir() + "write_hs_r_seeded_s.csr";
    // Seed a stale frame from a "previous run".
    {
        std::ofstream ofs(filename.c_str());
        ofs << "Ionic_Step 1\nstale content\n";
    }

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
    options.istep = 2;
    options.out_app_flag = 0;   // final-snapshot writer truncates
    options.reduce = false;
    ModuleIO::save_sparse(sparse_matrix, all_R_coor, pv, options);

    const std::string output = read_file(filename);
    // istep = 2 is the third ionic step, printed one-based as 3
    EXPECT_TRUE(starts_with(output, "Ionic_Step 3  # ABACUS version:"));
    EXPECT_EQ(output.find("stale content"), std::string::npos);

    std::remove(filename.c_str());
}

TEST(HsSparseIo, SparseCsrHeaderWritesCellAndBanner)
{
    const std::string filename = testing::TempDir() + "write_hs_r_csr_header_s.csr";
    std::remove(filename.c_str());

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    std::set<Abfs::Vector3_Order<int>> all_R_coor;
    all_R_coor.insert(r_vector);
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> sparse_matrix;
    sparse_matrix[r_vector][0][1] = 0.5;
    sparse_matrix[r_vector][1][1] = 1.5;

    UnitCell ucell;
    init_unitcell(ucell);

    ModuleIO::SparseWriteOptions options;
    options.filename = filename;
    options.label = "T";
    options.threshold = 1e-10;
    options.binary = false;
    options.istep = 0;
    options.reduce = false;
    options.csr_format = true;
    options.ucell = &ucell;
    ModuleIO::save_sparse(sparse_matrix, all_R_coor, pv, options);

    const std::vector<std::string> lines = read_lines(filename);
    ASSERT_GE(lines.size(), 12);
    // istep = 0 is the first ionic step, printed one-based as 1
    EXPECT_EQ(lines[0].compare(0, 12, "Ionic_Step 1"), 0);
    EXPECT_EQ(lines[1], " # print T matrix in real space T(R)");
    EXPECT_EQ(lines[2], " 1 # number of spin directions");
    EXPECT_EQ(lines[3], " 1 # spin index");
    EXPECT_EQ(lines[4], " 2 # number of localized basis");
    EXPECT_EQ(lines[5], " 1 # number of Bravais lattice vector R");
    EXPECT_EQ(lines[6], "");

    const std::string output = read_file(filename);
    EXPECT_THAT(output, testing::HasSubstr("user_defined_lattice"));
    EXPECT_THAT(output, testing::HasSubstr("Direct\n 0 0.25 0.5\n"));
    EXPECT_THAT(output, testing::HasSubstr("CSR Format"));
    EXPECT_THAT(output, testing::HasSubstr("0 0 0 2\n"));
    EXPECT_THAT(output, testing::HasSubstr(" # CSR values"));
    EXPECT_THAT(output, testing::HasSubstr(" # CSR column indices"));
    EXPECT_THAT(output, testing::HasSubstr(" # CSR row pointers"));

    std::remove(filename.c_str());
}

TEST(HsSparseIo, SparseCsrWritesEmptyRBlocks)
{
    const std::string filename = testing::TempDir() + "write_hs_r_csr_empty_s.csr";
    std::remove(filename.c_str());

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    const Abfs::Vector3_Order<int> filled_R(0, 0, 0);
    const Abfs::Vector3_Order<int> empty_R(1, 0, 0);
    std::set<Abfs::Vector3_Order<int>> all_R_coor;
    all_R_coor.insert(filled_R);
    all_R_coor.insert(empty_R);
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> sparse_matrix;
    sparse_matrix[filled_R][0][1] = 0.5;

    UnitCell ucell;
    init_unitcell(ucell);

    ModuleIO::SparseWriteOptions options;
    options.filename = filename;
    options.label = "T";
    options.threshold = 1e-10;
    options.binary = false;
    options.istep = 0;
    options.reduce = false;
    options.csr_format = true;
    options.ucell = &ucell;
    options.nspin = 2;
    options.ispin = 1;
    ModuleIO::save_sparse(sparse_matrix, all_R_coor, pv, options);

    const std::vector<std::string> lines = read_lines(filename);
    ASSERT_GE(lines.size(), 12);
    EXPECT_EQ(lines[2], " 2 # number of spin directions");
    EXPECT_EQ(lines[3], " 2 # spin index");
    EXPECT_EQ(lines[4], " 2 # number of localized basis");
    // Every R vector of the list is announced, including the empty one.
    EXPECT_EQ(lines[5], " 2 # number of Bravais lattice vector R");

    const std::string output = read_file(filename);
    EXPECT_THAT(output, testing::HasSubstr(" 0 0 0 1\n # CSR values\n"));
    // An R block without elements above the threshold is still written with
    // the same four-line record as the H(R)/S(R) writer.
    EXPECT_THAT(output, testing::HasSubstr(
        " 1 0 0 0\n # CSR values\n\n # CSR column indices\n\n # CSR row pointers\n 0 0 0\n"));

    std::remove(filename.c_str());
}

TEST(HsSparseIo, SparseCsrWritesComplexValues)
{
    const std::string filename = testing::TempDir() + "write_hs_r_csr_complex_lx.csr";
    std::remove(filename.c_str());

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    std::set<Abfs::Vector3_Order<int>> all_R_coor;
    all_R_coor.insert(r_vector);
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, std::complex<double>>>> sparse_matrix;
    sparse_matrix[r_vector][0][0] = std::complex<double>(1.5, -2.5);
    // below the threshold, must not appear in the CSR payload
    sparse_matrix[r_vector][1][0] = std::complex<double>(1e-12, 0.0);

    UnitCell ucell;
    init_unitcell(ucell);

    ModuleIO::SparseWriteOptions options;
    options.filename = filename;
    options.label = "Lx";
    options.threshold = 1e-10;
    options.binary = false;
    options.precision = 8;
    options.istep = 0;
    options.reduce = false;
    options.csr_format = true;
    options.ucell = &ucell;
    ModuleIO::save_sparse(sparse_matrix, all_R_coor, pv, options);

    const std::vector<std::string> lines = read_lines(filename);
    ASSERT_GE(lines.size(), 6);
    // istep = 0 is the first ionic step, printed one-based as 1
    EXPECT_EQ(lines[0].compare(0, 12, "Ionic_Step 1"), 0);
    EXPECT_EQ(lines[1], " # print Lx matrix in real space Lx(R)");
    EXPECT_EQ(lines[2], " 1 # number of spin directions");
    EXPECT_EQ(lines[3], " 1 # spin index");
    EXPECT_EQ(lines[4], " 2 # number of localized basis");
    EXPECT_EQ(lines[5], " 1 # number of Bravais lattice vector R");

    const std::string output = read_file(filename);
    EXPECT_THAT(output, testing::HasSubstr("0 0 0 1\n"));
    EXPECT_THAT(output, testing::HasSubstr(
        " # CSR values\n (1.50000000e+00,-2.50000000e+00)\n"));
    EXPECT_THAT(output, testing::HasSubstr(" # CSR column indices\n 0\n"));
    EXPECT_THAT(output, testing::HasSubstr(" # CSR row pointers\n 0 1 1\n"));

    std::remove(filename.c_str());
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
