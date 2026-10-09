/**
 * @file test_dhs_sparse_writer.cpp
 * @brief Unit tests for dhs_sparse_writer.cpp (ModuleIO::save_dH_sparse).
 *
 * save_dH_sparse writes the derivative sparse matrices (dHx/dHy/dHz) for
 * either H or S into separate CSR files.  The tests verify that:
 * - text mode respects the sparse threshold
 * - binary mode preserves the native field order
 * - SOC spin mode (nspin == 4) writes the complex variant for all
 *   three directions.
 *
 * dhs_sparse_writer.cpp only reads GlobalV::DRANK; no PARAM state is
 * involved, so the fixture sets DRANK directly.
 */
#include "csr_test_helpers.h"

#include "source_io/module_hs/dhs_sparse_writer.h"

#ifdef __MPI
#include <mpi.h>
#endif

TEST(DhsSparseWriter, TextCountsOnlyValuesAboveThreshold)
{
    remove_derivative_files("h", 5);
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

    ModuleIO::save_dH_sparse(5, pv, arrays, 1e-10, false, "h", 8, "./", "./", "scf", false, 1, 2);

    const std::vector<std::string> lines = read_lines("dhrxs1g6_nao.csr");
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

    const std::vector<std::string> y_lines = read_lines("dhrys1g6_nao.csr");
    ASSERT_GE(y_lines.size(), 4);
    EXPECT_EQ(y_lines[3], "0 0 0 0");

    remove_derivative_files("h", 5);
}

TEST(DhsSparseWriter, BinaryCountsOnlyValuesAboveThreshold)
{
    remove_derivative_files("h", 6);
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

    ModuleIO::save_dH_sparse(6, pv, arrays, 1e-10, true, "h", 8, "./", "./", "scf", false, 1, 2);

    std::ifstream ifs("dhrxs1g7_nao.csr", std::ios::binary);
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

    remove_derivative_files("h", 6);
}

TEST(DhsSparseWriter, SocWritesAllDirections)
{
    remove_derivative_files("s", 7);
    init_sparse_output_globals();

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    LCAO_HS_Arrays arrays;
    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    arrays.all_R_coor.insert(r_vector);
    arrays.dHRx_soc_sparse[r_vector][0][0] = std::complex<double>(1.0, 0.0);
    arrays.dHRy_soc_sparse[r_vector][0][1] = std::complex<double>(2.0, -1.0);
    arrays.dHRz_soc_sparse[r_vector][1][1] = std::complex<double>(-3.0, 0.5);

    ModuleIO::save_dH_sparse(7, pv, arrays, 1e-10, false, "s", 8, "./", "./", "scf", false, 4, 2);

    const std::string x_output = read_file("dsrxs1g8_nao.csr");
    const std::string y_output = read_file("dsrys1g8_nao.csr");
    const std::string z_output = read_file("dsrzs1g8_nao.csr");
    EXPECT_THAT(x_output, testing::HasSubstr("Matrix number of dHx(R): 1\n0 0 0 1\n"));
    EXPECT_THAT(y_output, testing::HasSubstr("Matrix number of dHy(R): 1\n0 0 0 1\n"));
    EXPECT_THAT(z_output, testing::HasSubstr("Matrix number of dHz(R): 1\n0 0 0 1\n"));
    EXPECT_THAT(x_output, testing::HasSubstr("(1.00000000e+00,0.00000000e+00)"));
    EXPECT_THAT(y_output, testing::HasSubstr("(2.00000000e+00,-1.00000000e+00)"));
    EXPECT_THAT(z_output, testing::HasSubstr("(-3.00000000e+00,5.00000000e-01)"));

    remove_derivative_files("s", 7);
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
