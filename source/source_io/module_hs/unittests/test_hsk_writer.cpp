/**
 * @file test_hsk_writer.cpp
 * @brief Unit tests for ModuleIO::write_hsk
 *
 * write_hsk loops over the k points owned by this rank, calls
 * Hamilt::updateHk followed by Hamilt::matrix, and writes one H(k) file
 * per k point and one S(k) file per spin-up k point (the spin-down
 * overlap is skipped because both spin channels share the same S).
 *
 * A lightweight FakeHamilt supplies controllable local H/S blocks; the
 * tests verify the call order (updateHk once per k point), the k-index
 * mapping (ik -> ik2iktot global index appears in the filename), the
 * spin-skip rule, and that the binary payload matches the distributed
 * input after MPI reduction.
 */
#include "source_io/module_hs/hsk_writer.h"

#include "source_base/matrix_block.h"
#include "source_base/module_out/filename.h"
#include "source_base/parallel_comm.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_hamilt/hamilt.h"

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
/// Rank in DIAG_WORLD; 0 without MPI.
int test_rank = 0;

/// Output options shared by the tests (binary, no append, step 0).
const int k_out_type = 2;
const bool k_out_app = false;
const int k_istep = 0;

/**
 * @brief Test double of hamilt::Hamilt: replays per-k-point local H/S
 *        buffers set by the test and counts updateHk calls.
 */
template <typename T>
class FakeHamilt : public hamilt::Hamilt<T>
{
  public:
    int dim = 0;
    int update_count = 0;
    std::vector<std::vector<T>> hk_local;
    std::vector<std::vector<T>> sk_local;

    void set_k_matrix(const int ik, std::vector<T> hk, std::vector<T> sk)
    {
        if (static_cast<int>(hk_local.size()) <= ik)
        {
            hk_local.resize(ik + 1);
            sk_local.resize(ik + 1);
        }
        hk_local[ik] = std::move(hk);
        sk_local[ik] = std::move(sk);
    }

    void updateHk(const int ik) override
    {
        ++update_count;
        current_ik_ = ik;
    }

    void matrix(ModuleBase::MatrixBlock<T>& hk,
                ModuleBase::MatrixBlock<T>& sk) override
    {
        EXPECT_GE(current_ik_, 0);
        hk = ModuleBase::MatrixBlock<T> {nullptr, 0, 0, nullptr};
        hk.p = hk_local[current_ik_].data();
        hk.row = static_cast<size_t>(dim);
        hk.col = static_cast<size_t>(dim);
        sk = ModuleBase::MatrixBlock<T> {nullptr, 0, 0, nullptr};
        sk.p = sk_local[current_ik_].data();
        sk.row = static_cast<size_t>(dim);
        sk.col = static_cast<size_t>(dim);
    }

  private:
    int current_ik_ = -1;
};

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
 * @brief Whether a regular file can be opened.
 */
bool file_exists(const std::string& filename)
{
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    return ifs.is_open();
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

/**
 * @brief Build the same filename write_hsk asks filename_output for.
 */
std::string expected_filename(const std::string& property,
                              const int ik,
                              const std::vector<int>& ik2iktot,
                              const int nspin,
                              const int nkstot,
                              const bool gamma_only)
{
    return ModuleIO::filename_output("", property, "nao", ik, ik2iktot,
                                     nspin, nkstot, k_out_type, k_out_app,
                                     gamma_only, k_istep);
}

/**
 * @brief Open the running-log stream on a per-rank file.
 */
std::ofstream open_running_log()
{
    const std::string filename
        = "write_hsk_running_" + std::to_string(test_rank) + ".log";
    return std::ofstream(filename.c_str());
}

} // namespace

// ---------------------------------------------------------------------------
// nspin = 1, multik, real matrices: both H(k) and S(k) are written for
// every k point; payload equals the upper triangle of the distributed input.
// ---------------------------------------------------------------------------

TEST(WriteHsk, WritesHkAndSkForEveryKPoint)
{
    const int dim = 3;
    const int nspin = 1;
    const int nks = 2;
    const int nkstot = 2;
    const std::vector<int> ik2iktot = {0, 1};
    const std::vector<int> isk = {0, 0};
    const bool gamma_only = false;

    Parallel_Orbitals pv;
    initialize_distribution(pv, dim);

    FakeHamilt<double> fake;
    fake.dim = dim;
    std::vector<std::vector<double>> hk_global(nks);
    std::vector<std::vector<double>> sk_global(nks);
    for (int ik = 0; ik < nks; ++ik)
    {
        hk_global[ik].resize(dim * dim);
        sk_global[ik].resize(dim * dim);
        for (int i = 0; i < dim * dim; ++i)
        {
            hk_global[ik][i] = (ik + 1) * 100.0 + i + 0.5;
            sk_global[ik][i] = -((ik + 1) * 100.0 + i) - 0.25;
        }
        fake.set_k_matrix(ik,
                          distribute_matrix(pv, hk_global[ik], dim),
                          distribute_matrix(pv, sk_global[ik], dim));
    }

    std::vector<std::string> filenames;
    for (int ik = 0; ik < nks; ++ik)
    {
        filenames.push_back(expected_filename("hk", ik, ik2iktot, nspin,
                                              nkstot, gamma_only));
        filenames.push_back(expected_filename("sk", ik, ik2iktot, nspin,
                                              nkstot, gamma_only));
    }
    for (size_t i = 0; i < filenames.size(); ++i)
    {
        remove_test_file(filenames[i]);
    }

    std::ofstream ofs = open_running_log();
    ModuleIO::write_hsk("", nspin, nks, nkstot, ik2iktot, isk, &fake, pv,
                        gamma_only, k_out_app, k_istep, k_out_type, 8, dim,
                        "cg", test_rank, ofs);

    EXPECT_EQ(fake.update_count, nks);
    if (test_rank == 0)
    {
        for (int ik = 0; ik < nks; ++ik)
        {
            const std::string h_fn = expected_filename("hk", ik, ik2iktot,
                                                       nspin, nkstot, gamma_only);
            std::ifstream hfs(h_fn.c_str(), std::ios::binary);
            ASSERT_TRUE(hfs.is_open()) << h_fn;
            EXPECT_EQ(read_record<double>(hfs, dim),
                      upper_triangle(hk_global[ik], dim));
            EXPECT_EQ(hfs.peek(), std::ifstream::traits_type::eof());

            const std::string s_fn = expected_filename("sk", ik, ik2iktot,
                                                       nspin, nkstot, gamma_only);
            std::ifstream sfs(s_fn.c_str(), std::ios::binary);
            ASSERT_TRUE(sfs.is_open()) << s_fn;
            EXPECT_EQ(read_record<double>(sfs, dim),
                      upper_triangle(sk_global[ik], dim));
            EXPECT_EQ(sfs.peek(), std::ifstream::traits_type::eof());
        }
    }

    for (size_t i = 0; i < filenames.size(); ++i)
    {
        remove_test_file(filenames[i]);
    }
    ofs.close();
    remove_test_file("write_hsk_running_" + std::to_string(test_rank) + ".log");
}

// ---------------------------------------------------------------------------
// nspin = 2: H(k) is written for both spin channels, but S(k) is skipped
// for the spin-down k point (isk == 1).
// ---------------------------------------------------------------------------

TEST(WriteHsk, SkipsOverlapForSpinDownChannel)
{
    const int dim = 2;
    const int nspin = 2;
    const int nks = 2;
    const int nkstot = 2;
    const std::vector<int> ik2iktot = {0, 1};
    const std::vector<int> isk = {0, 1};
    const bool gamma_only = false;

    Parallel_Orbitals pv;
    initialize_distribution(pv, dim);

    FakeHamilt<double> fake;
    fake.dim = dim;
    for (int ik = 0; ik < nks; ++ik)
    {
        std::vector<double> hk(dim * dim, static_cast<double>(ik + 1));
        std::vector<double> sk(dim * dim, -static_cast<double>(ik + 1));
        fake.set_k_matrix(ik, distribute_matrix(pv, hk, dim),
                          distribute_matrix(pv, sk, dim));
    }

    const std::string hk_up = expected_filename("hk", 0, ik2iktot, nspin,
                                                nkstot, gamma_only);
    const std::string hk_dn = expected_filename("hk", 1, ik2iktot, nspin,
                                                nkstot, gamma_only);
    const std::string sk_up = expected_filename("sk", 0, ik2iktot, nspin,
                                                nkstot, gamma_only);
    const std::string sk_dn = expected_filename("sk", 1, ik2iktot, nspin,
                                                nkstot, gamma_only);
    std::vector<std::string> cleanup = {hk_up, hk_dn, sk_up, sk_dn};
    for (size_t i = 0; i < cleanup.size(); ++i)
    {
        remove_test_file(cleanup[i]);
    }

    std::ofstream ofs = open_running_log();
    ModuleIO::write_hsk("", nspin, nks, nkstot, ik2iktot, isk, &fake, pv,
                        gamma_only, k_out_app, k_istep, k_out_type, 8, dim,
                        "cg", test_rank, ofs);

    EXPECT_EQ(fake.update_count, nks);
    if (test_rank == 0)
    {
        EXPECT_TRUE(file_exists(hk_up));
        EXPECT_TRUE(file_exists(hk_dn));
        EXPECT_TRUE(file_exists(sk_up));
        EXPECT_FALSE(file_exists(sk_dn));
    }

    for (size_t i = 0; i < cleanup.size(); ++i)
    {
        remove_test_file(cleanup[i]);
    }
    ofs.close();
    remove_test_file("write_hsk_running_" + std::to_string(test_rank) + ".log");
}

// ---------------------------------------------------------------------------
// Complex multik case: write_hsk<std::complex<double>> must instantiate and
// write the packed upper triangle of complex double pairs.
// ---------------------------------------------------------------------------

TEST(WriteHsk, WritesComplexKMatrix)
{
    const int dim = 2;
    const int nspin = 1;
    const int nks = 1;
    const int nkstot = 1;
    const std::vector<int> ik2iktot = {0};
    const std::vector<int> isk = {0};
    const bool gamma_only = false;

    Parallel_Orbitals pv;
    initialize_distribution(pv, dim);

    std::vector<std::complex<double>> hk_global(dim * dim);
    std::vector<std::complex<double>> sk_global(dim * dim);
    for (int i = 0; i < dim * dim; ++i)
    {
        hk_global[i] = std::complex<double>(i + 0.25, -i - 0.75);
        sk_global[i] = std::complex<double>(-i - 0.5, i + 0.125);
    }

    FakeHamilt<std::complex<double>> fake;
    fake.dim = dim;
    fake.set_k_matrix(0, distribute_matrix(pv, hk_global, dim),
                      distribute_matrix(pv, sk_global, dim));

    const std::string h_fn = expected_filename("hk", 0, ik2iktot, nspin,
                                               nkstot, gamma_only);
    const std::string s_fn = expected_filename("sk", 0, ik2iktot, nspin,
                                               nkstot, gamma_only);
    remove_test_file(h_fn);
    remove_test_file(s_fn);

    std::ofstream ofs = open_running_log();
    ModuleIO::write_hsk<std::complex<double>>("", nspin, nks, nkstot,
                                              ik2iktot, isk, &fake, pv,
                                              gamma_only, k_out_app, k_istep,
                                              k_out_type, 8, dim, "cg",
                                              test_rank, ofs);

    EXPECT_EQ(fake.update_count, nks);
    if (test_rank == 0)
    {
        std::ifstream hfs(h_fn.c_str(), std::ios::binary | std::ios::ate);
        ASSERT_TRUE(hfs.is_open());
        const int element_count = dim * (dim + 1) / 2;
        EXPECT_EQ(hfs.tellg(),
                  static_cast<std::streamoff>(sizeof(int)
                                              + element_count * sizeof(std::complex<double>)));
        hfs.seekg(0);
        EXPECT_EQ((read_record<std::complex<double>>(hfs, dim)),
                  upper_triangle(hk_global, dim));

        std::ifstream sfs(s_fn.c_str(), std::ios::binary);
        ASSERT_TRUE(sfs.is_open());
        EXPECT_EQ((read_record<std::complex<double>>(sfs, dim)),
                  upper_triangle(sk_global, dim));
    }

    remove_test_file(h_fn);
    remove_test_file(s_fn);
    ofs.close();
    remove_test_file("write_hsk_running_" + std::to_string(test_rank) + ".log");
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
