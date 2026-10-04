/**
 * @file test_hs_dense_io.cpp
 * @brief Unit tests for ModuleIO::save_mat
 *
 * save_mat writes a square matrix in text or binary upper-triangle/full
 * format. When built with __MPI (the project default), save_mat uses MPI
 * collectives and MPI_Barrier. This test provides a main() that calls
 * MPI_Init first so the MPI path is safe on a single rank. Header-only
 * tests use a default-constructed/serial Parallel_2D and do not check
 * values; the binary upper-triangle tests distribute a known matrix over
 * a real Parallel_Orbitals and verify the gathered payload.
 */
#include <gtest/gtest.h>

#include "source_io/module_hs/hs_dense_io.h"
#include "source_base/parallel_2d.h"
#include "source_base/parallel_comm.h"
#include "source_basis/module_ao/parallel_orbitals.h"

#include <complex>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#ifdef __MPI
#include <mpi.h>
#endif

namespace
{
/// Rank in DIAG_WORLD; 0 without MPI.
int test_rank = 0;

/**
 * @brief Read the entire content of a text file into a string.
 */
std::string read_text_file(const std::string& filename)
{
    std::ifstream ifs(filename);
    if (!ifs.is_open())
    {
        return "";
    }
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return oss.str();
}

/**
 * @brief Read just the dimension int from a binary matrix file.
 */
int read_binary_dim(const std::string& filename)
{
    FILE* fp = std::fopen(filename.c_str(), "rb");
    if (fp == nullptr)
    {
        return -1;
    }
    int dim = 0;
    std::fread(&dim, sizeof(int), 1, fp);
    std::fclose(fp);
    return dim;
}

/**
 * @brief Create a serial Parallel_2D that owns all rows/cols on rank 0.
 */
Parallel_2D make_serial_pv(int dim)
{
    Parallel_2D pv;
    pv.set_serial(dim, dim);
    return pv;
}

/**
 * @brief Set up a Parallel_Orbitals distribution over DIAG_WORLD.
 */
void initialize_distribution(Parallel_Orbitals& pv, const int dim)
{
#ifdef __MPI
    ASSERT_EQ(pv.init(dim, dim, 1, DIAG_WORLD), 0);
#else
    pv.set_serial(dim, dim);
#endif
}

/**
 * @brief Extract the block of a global square matrix owned by this rank.
 */
template <typename T>
std::vector<T> distribute_matrix(const Parallel_Orbitals& pv,
                                 const std::vector<T>& global,
                                 const int dim)
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

/**
 * @brief Pack the upper triangle of a square matrix row by row.
 */
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

/**
 * @brief Read one binary record: dim header followed by the packed
 *        upper-triangle payload.
 */
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

/**
 * @brief Remove a test file on rank 0 and sync all ranks.
 */
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

// ---------------------------------------------------------------------------
// Text mode: verify header fields (MPI path, drank=0)
// ---------------------------------------------------------------------------

TEST(HsDenseIoSaveMat, TextHeaderContainsStepAndDimension)
{
    const std::string filename = "test_save_mat_text_header.dat";
    const int dim = 3;
    const std::vector<double> mat(dim * dim, 0.0);
    const Parallel_2D pv = make_serial_pv(dim);
    ModuleIO::save_mat<double>(0, mat.data(), dim, false, 8, false, false,
                               filename, pv, 0, "genelpa", false);
    const std::string content = read_text_file(filename);
    // MPI text path writes header lines starting with '#'
    EXPECT_NE(content.find("ionic step"), std::string::npos);
    EXPECT_NE(content.find("1"), std::string::npos); // istep+1 = 1
    EXPECT_NE(content.find("gamma only"), std::string::npos);
    EXPECT_NE(content.find("rows"), std::string::npos);
    EXPECT_NE(content.find("3"), std::string::npos); // dim
    std::remove(filename.c_str());
}

TEST(HsDenseIoSaveMat, TextHeaderStepTwo)
{
    const std::string filename = "test_save_mat_text_header_step2.dat";
    const int dim = 2;
    const std::vector<double> mat(dim * dim, 0.0);
    const Parallel_2D pv = make_serial_pv(dim);
    ModuleIO::save_mat<double>(1, mat.data(), dim, false, 8, false, false,
                               filename, pv, 0, "genelpa", false);
    const std::string content = read_text_file(filename);
    // istep=1 -> "ionic step 2"
    EXPECT_NE(content.find("ionic step 2"), std::string::npos);
    std::remove(filename.c_str());
}

TEST(HsDenseIoSaveMat, TextHeaderGammaOnlyDouble)
{
    const std::string filename = "test_save_mat_text_gamma.dat";
    const int dim = 2;
    const std::vector<double> mat(dim * dim, 0.0);
    const Parallel_2D pv = make_serial_pv(dim);
    ModuleIO::save_mat<double>(0, mat.data(), dim, false, 8, false, false,
                               filename, pv, 0, "genelpa", false);
    const std::string content = read_text_file(filename);
    // gamma_only is true for double (std::is_same<T, double>::value)
    EXPECT_NE(content.find("gamma only 1"), std::string::npos);
    std::remove(filename.c_str());
}

TEST(HsDenseIoSaveMat, TextHeaderGammaOnlyComplex)
{
    const std::string filename = "test_save_mat_text_gamma_complex.dat";
    const int dim = 2;
    const std::vector<std::complex<double>> mat(dim * dim);
    const Parallel_2D pv = make_serial_pv(dim);
    ModuleIO::save_mat<std::complex<double>>(0, mat.data(), dim, false, 8,
                               false, false, filename, pv, 0, "genelpa", false);
    const std::string content = read_text_file(filename);
    // gamma_only is false for complex<double>
    EXPECT_NE(content.find("gamma only 0"), std::string::npos);
    std::remove(filename.c_str());
}

// ---------------------------------------------------------------------------
// Binary mode: verify dim header (MPI path, drank=0)
// ---------------------------------------------------------------------------

TEST(HsDenseIoSaveMat, BinaryDimHeaderDouble)
{
    const std::string filename = "test_save_mat_bin_dim_double.dat";
    const int dim = 4;
    const std::vector<double> mat(dim * dim, 0.0);
    const Parallel_2D pv = make_serial_pv(dim);
    ModuleIO::save_mat<double>(0, mat.data(), dim, true, 8, false, false,
                               filename, pv, 0, "genelpa", false);
    const int dim_out = read_binary_dim(filename);
    EXPECT_EQ(dim_out, dim);
    std::remove(filename.c_str());
}

TEST(HsDenseIoSaveMat, BinaryDimHeaderComplex)
{
    const std::string filename = "test_save_mat_bin_dim_complex.dat";
    const int dim = 3;
    const std::vector<std::complex<double>> mat(dim * dim);
    const Parallel_2D pv = make_serial_pv(dim);
    ModuleIO::save_mat<std::complex<double>>(0, mat.data(), dim, true, 8,
                               false, false, filename, pv, 0, "genelpa", false);
    const int dim_out = read_binary_dim(filename);
    EXPECT_EQ(dim_out, dim);
    std::remove(filename.c_str());
}

// ---------------------------------------------------------------------------
// Binary upper triangle: distribute a known matrix over a Parallel_Orbitals
// and verify the rank-0 file contains the complete gathered payload.
// ---------------------------------------------------------------------------

TEST(HsDenseIoSaveMat, BinaryGammaUpperTriangle)
{
    const int dim = 3;
    const std::string filename = "save_mat_gamma_upper.dat";
    remove_test_file(filename);

    Parallel_Orbitals pv;
    initialize_distribution(pv, dim);
    std::vector<double> global(dim * dim);
    for (int i = 0; i < dim * dim; ++i)
    {
        global[i] = i + 0.5;
    }
    const std::vector<double> local = distribute_matrix(pv, global, dim);

    ModuleIO::save_mat(0, local.data(), dim, true, 8, true, false, filename,
                       pv, test_rank, "cg");

    if (test_rank == 0)
    {
        std::ifstream ifs(filename.c_str(), std::ios::binary | std::ios::ate);
        ASSERT_TRUE(ifs.is_open());
        EXPECT_EQ(ifs.tellg(),
                  static_cast<std::streamoff>(sizeof(int) + 6 * sizeof(double)));
        ifs.seekg(0);
        EXPECT_EQ(read_record<double>(ifs, dim), upper_triangle(global, dim));
        EXPECT_EQ(ifs.peek(), std::ifstream::traits_type::eof());
    }

    remove_test_file(filename);
}

TEST(HsDenseIoSaveMat, BinaryComplexUpperTriangleMpiReduction)
{
    const int dim = 4;
    const std::string filename = "save_mat_complex_upper.dat";
    remove_test_file(filename);

    Parallel_Orbitals pv;
    initialize_distribution(pv, dim);
    std::vector<std::complex<double>> global(dim * dim);
    for (int i = 0; i < dim * dim; ++i)
    {
        global[i] = std::complex<double>(i + 0.25, -i - 0.75);
    }
    const std::vector<std::complex<double>> local
        = distribute_matrix(pv, global, dim);

    ModuleIO::save_mat(0, local.data(), dim, true, 8, true, false, filename,
                       pv, test_rank, "cg");

    if (test_rank == 0)
    {
        std::ifstream ifs(filename.c_str(), std::ios::binary | std::ios::ate);
        ASSERT_TRUE(ifs.is_open());
        const int element_count = dim * (dim + 1) / 2;
        EXPECT_EQ(ifs.tellg(),
                  static_cast<std::streamoff>(sizeof(int)
                                              + element_count * sizeof(std::complex<double>)));
        ifs.seekg(0);
        EXPECT_EQ((read_record<std::complex<double>>(ifs, dim)),
                  upper_triangle(global, dim));
        EXPECT_EQ(ifs.peek(), std::ifstream::traits_type::eof());
    }

    remove_test_file(filename);
}

TEST(HsDenseIoSaveMat, BinaryAppendCompleteRecordsAndOverwrite)
{
    const int dim = 2;
    const std::string filename = "save_mat_append.dat";
    remove_test_file(filename);

    Parallel_Orbitals pv;
    initialize_distribution(pv, dim);
    const std::vector<double> first = {1.0, 2.0, 3.0, 4.0};
    const std::vector<double> second = {5.0, 6.0, 7.0, 8.0};
    const std::vector<double> replacement = {9.0, 10.0, 11.0, 12.0};
    const std::vector<double> first_local = distribute_matrix(pv, first, dim);
    const std::vector<double> second_local = distribute_matrix(pv, second, dim);
    const std::vector<double> replacement_local
        = distribute_matrix(pv, replacement, dim);

    ModuleIO::save_mat(0, first_local.data(), dim, true, 8, true, true,
                       filename, pv, test_rank, "cg");
    ModuleIO::save_mat(1, second_local.data(), dim, true, 8, true, true,
                       filename, pv, test_rank, "cg");

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
    ModuleIO::save_mat(2, replacement_local.data(), dim, true, 8, true, false,
                       filename, pv, test_rank, "cg");

    if (test_rank == 0)
    {
        std::ifstream ifs(filename.c_str(), std::ios::binary);
        ASSERT_TRUE(ifs.is_open());
        EXPECT_EQ(read_record<double>(ifs, dim), upper_triangle(replacement, dim));
        EXPECT_EQ(ifs.peek(), std::ifstream::traits_type::eof());
    }

    remove_test_file(filename);
}

// ---------------------------------------------------------------------------
// Main: initialize MPI before running tests (required by save_mat's MPI path)
// ---------------------------------------------------------------------------

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
