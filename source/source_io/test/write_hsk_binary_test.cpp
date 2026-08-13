#include "source_io/module_hs/write_hs.h"

#include "gtest/gtest.h"

#include <complex>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

namespace
{

int test_rank = 0;

void initialize_distribution(Parallel_Orbitals& pv, const int dim)
{
#ifdef __MPI
    ASSERT_EQ(pv.init(dim, dim, 1, DIAG_WORLD), 0);
#else
    pv.set_serial(dim, dim);
#endif
}

template <typename T>
std::vector<T> distribute_matrix(const Parallel_Orbitals& pv, const std::vector<T>& global, const int dim)
{
    std::vector<T> local(pv.get_local_size(), T());
    for (int i = 0; i < dim; ++i)
    {
        const int ir = pv.global2local_row(i);
        if (ir < 0)
        {
            continue;
        }
        for (int j = 0; j < dim; ++j)
        {
            const int ic = pv.global2local_col(j);
            if (ic >= 0)
            {
                local[ir * pv.ncol + ic] = global[i * dim + j];
            }
        }
    }
    return local;
}

template <typename T>
std::vector<T> upper_triangle(const std::vector<T>& matrix, const int dim)
{
    std::vector<T> values;
    values.reserve(dim * (dim + 1) / 2);
    for (int i = 0; i < dim; ++i)
    {
        for (int j = i; j < dim; ++j)
        {
            values.push_back(matrix[i * dim + j]);
        }
    }
    return values;
}

template <typename T>
std::vector<T> read_record(std::ifstream& ifs, const int expected_dim)
{
    int dim = 0;
    ifs.read(reinterpret_cast<char*>(&dim), sizeof(int));
    EXPECT_TRUE(ifs.good());
    EXPECT_EQ(dim, expected_dim);

    std::vector<T> values(expected_dim * (expected_dim + 1) / 2);
    ifs.read(reinterpret_cast<char*>(values.data()), values.size() * sizeof(T));
    EXPECT_TRUE(ifs.good());
    return values;
}

void remove_test_file(const std::string& filename)
{
    if (test_rank == 0)
    {
        std::remove(filename.c_str());
    }
#ifdef __MPI
    MPI_Barrier(DIAG_WORLD);
#endif
}

} // namespace

TEST(WriteHskBinary, WritesGammaUpperTriangle)
{
    const int dim = 3;
    const std::string filename = "write_hsk_gamma.dat";
    remove_test_file(filename);

    Parallel_Orbitals pv;
    initialize_distribution(pv, dim);
    std::vector<double> global(dim * dim);
    for (int i = 0; i < dim * dim; ++i)
    {
        global[i] = i + 0.5;
    }
    const std::vector<double> local = distribute_matrix(pv, global, dim);

    ModuleIO::save_mat(0, local.data(), dim, true, 8, true, false, filename, pv, test_rank);

    if (test_rank == 0)
    {
        std::ifstream ifs(filename.c_str(), std::ios::binary | std::ios::ate);
        ASSERT_TRUE(ifs.is_open());
        EXPECT_EQ(ifs.tellg(), static_cast<std::streamoff>(sizeof(int) + 6 * sizeof(double)));
        ifs.seekg(0);
        EXPECT_EQ(read_record<double>(ifs, dim), upper_triangle(global, dim));
        EXPECT_EQ(ifs.peek(), std::ifstream::traits_type::eof());
    }

    remove_test_file(filename);
}

TEST(WriteHskBinary, WritesComplexUpperTriangleWithMpiReduction)
{
    const int dim = 4;
    const std::string filename = "write_hsk_complex.dat";
    remove_test_file(filename);

    Parallel_Orbitals pv;
    initialize_distribution(pv, dim);
    std::vector<std::complex<double>> global(dim * dim);
    for (int i = 0; i < dim * dim; ++i)
    {
        global[i] = std::complex<double>(i + 0.25, -i - 0.75);
    }
    const std::vector<std::complex<double>> local = distribute_matrix(pv, global, dim);

    ModuleIO::save_mat(0, local.data(), dim, true, 8, true, false, filename, pv, test_rank);

    if (test_rank == 0)
    {
        std::ifstream ifs(filename.c_str(), std::ios::binary | std::ios::ate);
        ASSERT_TRUE(ifs.is_open());
        const int element_count = dim * (dim + 1) / 2;
        EXPECT_EQ(ifs.tellg(), static_cast<std::streamoff>(sizeof(int) + element_count * sizeof(std::complex<double>)));
        ifs.seekg(0);
        EXPECT_EQ(read_record<std::complex<double>>(ifs, dim), upper_triangle(global, dim));
        EXPECT_EQ(ifs.peek(), std::ifstream::traits_type::eof());
    }

    remove_test_file(filename);
}

TEST(WriteHskBinary, AppendsCompleteRecordsAndCanOverwrite)
{
    const int dim = 2;
    const std::string filename = "write_hsk_append.dat";
    remove_test_file(filename);

    Parallel_Orbitals pv;
    initialize_distribution(pv, dim);
    const std::vector<double> first = {1.0, 2.0, 3.0, 4.0};
    const std::vector<double> second = {5.0, 6.0, 7.0, 8.0};
    const std::vector<double> replacement = {9.0, 10.0, 11.0, 12.0};
    const std::vector<double> first_local = distribute_matrix(pv, first, dim);
    const std::vector<double> second_local = distribute_matrix(pv, second, dim);
    const std::vector<double> replacement_local = distribute_matrix(pv, replacement, dim);

    ModuleIO::save_mat(0, first_local.data(), dim, true, 8, true, true, filename, pv, test_rank);
    ModuleIO::save_mat(1, second_local.data(), dim, true, 8, true, true, filename, pv, test_rank);

    if (test_rank == 0)
    {
        std::ifstream ifs(filename.c_str(), std::ios::binary);
        ASSERT_TRUE(ifs.is_open());
        EXPECT_EQ(read_record<double>(ifs, dim), upper_triangle(first, dim));
        EXPECT_EQ(read_record<double>(ifs, dim), upper_triangle(second, dim));
        EXPECT_EQ(ifs.peek(), std::ifstream::traits_type::eof());
    }

#ifdef __MPI
    MPI_Barrier(DIAG_WORLD);
#endif
    ModuleIO::save_mat(2, replacement_local.data(), dim, true, 8, true, false, filename, pv, test_rank);

    if (test_rank == 0)
    {
        std::ifstream ifs(filename.c_str(), std::ios::binary);
        ASSERT_TRUE(ifs.is_open());
        EXPECT_EQ(read_record<double>(ifs, dim), upper_triangle(replacement, dim));
        EXPECT_EQ(ifs.peek(), std::ifstream::traits_type::eof());
    }

    remove_test_file(filename);
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    DIAG_WORLD = MPI_COMM_WORLD;
    MPI_Comm_rank(DIAG_WORLD, &test_rank);
#else
    test_rank = 0;
#endif

    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
