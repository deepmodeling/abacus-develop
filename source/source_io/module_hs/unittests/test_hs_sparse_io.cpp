/**
 * @file test_hs_sparse_io.cpp
 * @brief Unit tests for hs_sparse_io.cpp (ModuleIO::save_sparse).
 *
 * The tests pin down the legacy CSR text/binary layout written by
 * save_sparse: header style ("STEP: n"), dimension lines, threshold
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

    GlobalV::DRANK = 0;

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
    EXPECT_TRUE(starts_with(output, "STEP: 0\n"));
    EXPECT_THAT(output, testing::HasSubstr("Matrix Dimension of S(R): 2\n"));
    EXPECT_THAT(output, testing::HasSubstr("Matrix number of S(R): 1\n"));
    EXPECT_THAT(output, testing::HasSubstr("0 0 0 2\n"));

    std::remove(filename.c_str());
}

TEST(HsSparseIo, SparseTextCountsOnlyValuesAboveThreshold)
{
    const std::string filename = "write_hs_r_threshold_s.csr";
    std::remove(filename.c_str());

    GlobalV::DRANK = 0;

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

TEST(HsSparseIo, SparseBinaryHeaderWritesConcreteStep)
{
    const std::string filename = "write_hs_r_legacy_binary_s.csr";
    std::remove(filename.c_str());

    GlobalV::DRANK = 0;

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

    const std::vector<int> header_and_r = read_binary_ints(filename, 7);
    EXPECT_THAT(header_and_r, testing::ElementsAre(3, 2, 1, 0, 0, 0, 2));

    std::remove(filename.c_str());
}

TEST(HsSparseIo, SparseBinaryCountsOnlyValuesAboveThreshold)
{
    const std::string filename = "write_hs_r_threshold_binary_s.csr";
    std::remove(filename.c_str());

    GlobalV::DRANK = 0;

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
