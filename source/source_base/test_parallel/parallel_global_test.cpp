#ifdef __MPI
#include "source_base/parallel_global.h"

#include "mpi.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <complex>
#include <cstring>
#include <string>
#include <unistd.h>

#include "source_base/global_variable.h"

/************************************************
 *  unit test of functions in parallel_global.cpp
 ***********************************************/

/**
 * The tested functions are:
 *   i. Parallel_Global::split_diag_world(), which is
 *   used in David diagonalization in pw basis
 *   calculation.
 *   ii. Parallel_Global::split_grid_world()
 *   iii. Parallel_Global::MyProd(std::complex<double> *in,std::complex<double> *inout,int *len,MPI_Datatype *dptr);
 *   iv. Parallel_Global::init_pools();
 *   v. Parallel_Global::divide_pools(void);
 */

class MPIContext
{
  public:
    MPIContext()
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &_size);
    }

    int GetRank() const
    {
        return _rank;
    }
    int GetSize() const
    {
        return _size;
    }

    int drank;
    int dsize;
    int dcolor;

    int grank;
    int gsize;

    int kpar;
    int nproc_in_pool;
    int my_pool;
    int rank_in_pool;

    int nstogroup;
    int MY_BNDGROUP;
    int rank_in_stogroup;
    int nproc_in_stogroup;

  private:
    int _rank;
    int _size;
};

// --- Normal Test ---
class ParaGlobal : public ::testing::Test
{
  protected:
    MPIContext mpi;
    int nproc;
    int my_rank;
    void SetUp() override
    {
        nproc = mpi.GetSize();
        my_rank = mpi.GetRank();
    }
};


TEST_F(ParaGlobal, SplitGrid)
{
    // NPROC is set to 4 in parallel_global_test.sh
    if (nproc == 4)
    {
        Parallel_Global::split_grid_world(2, nproc, my_rank, mpi.grank, mpi.gsize);
        EXPECT_EQ(mpi.gsize, 2);
        if (my_rank == 0)
            EXPECT_EQ(mpi.grank, 0);
        if (my_rank == 1)
            EXPECT_EQ(mpi.grank, 1);
        if (my_rank == 2)
            EXPECT_EQ(mpi.grank, 0);
        if (my_rank == 3)
            EXPECT_EQ(mpi.grank, 1);
        Parallel_Global::split_grid_world(4, nproc, my_rank, mpi.grank, mpi.gsize);
        EXPECT_EQ(mpi.gsize, 1);
        if (my_rank == 0)
            EXPECT_EQ(mpi.grank, 0);
        if (my_rank == 1)
            EXPECT_EQ(mpi.grank, 0);
        if (my_rank == 2)
            EXPECT_EQ(mpi.grank, 0);
        if (my_rank == 3)
            EXPECT_EQ(mpi.grank, 0);
    }
    else
    {
        Parallel_Global::split_grid_world(nproc, nproc, my_rank, mpi.grank, mpi.gsize);
        EXPECT_EQ(mpi.gsize, 1);
        EXPECT_EQ(mpi.grank, 0);
    }
    // std::cout<<my_rank<<" "<<nproc<<" ";
    // std::cout<<mpi.grank<<" "<<mpi.gsize<<std::endl;
}

TEST_F(ParaGlobal, SplitDiag)
{
    // NPROC is set to 4 in parallel_global_test.sh
    if (nproc == 4)
    {
        Parallel_Global::split_diag_world(2, nproc, my_rank, mpi.drank, mpi.dsize, mpi.dcolor);
        EXPECT_EQ(mpi.dsize, 2);
        if (my_rank == 0)
            EXPECT_EQ(mpi.drank, 0);
        if (my_rank == 1)
            EXPECT_EQ(mpi.drank, 0);
        if (my_rank == 2)
            EXPECT_EQ(mpi.drank, 1);
        if (my_rank == 3)
            EXPECT_EQ(mpi.drank, 1);
        Parallel_Global::split_diag_world(4, nproc, my_rank, mpi.drank, mpi.dsize, mpi.dcolor);
        EXPECT_EQ(mpi.dsize, 4);
        if (my_rank == 0)
            EXPECT_EQ(mpi.drank, 0);
        if (my_rank == 1)
            EXPECT_EQ(mpi.drank, 1);
        if (my_rank == 2)
            EXPECT_EQ(mpi.drank, 2);
        if (my_rank == 3)
            EXPECT_EQ(mpi.drank, 3);
    }
    else
    {
        Parallel_Global::split_diag_world(nproc, nproc, my_rank, mpi.drank, mpi.dsize, mpi.dcolor);
        EXPECT_EQ(mpi.dsize, nproc);
    }
    // std::cout<<my_rank<<" "<<nproc<<" ";
    // std::cout<<mpi.drank<<" "<<mpi.dsize<<std::endl;
}

TEST_F(ParaGlobal, MyProd)
{
    std::complex<double> in[2] = {std::complex<double>(1.0, 2.0), std::complex<double>(-1, -2)};
    std::complex<double> inout[2] = {std::complex<double>(2.0, 1.0), std::complex<double>(-2, -1)};

    int len = 2;
    MPI_Datatype dptr = MPI_DOUBLE_COMPLEX;
    Parallel_Global::myProd(in, inout, &len, &dptr);
    EXPECT_EQ(inout[0], std::complex<double>(3.0, 3.0));
    EXPECT_EQ(inout[1], std::complex<double>(-3.0, -3.0));
}



TEST_F(ParaGlobal, DivideMPIPools)
{
    this->nproc = 12;
    mpi.kpar = 3;
    this->my_rank = 5;
    Parallel_Global::divide_mpi_groups(this->nproc,
                                       mpi.kpar,
                                       this->my_rank,
                                       mpi.nproc_in_pool,
                                       mpi.my_pool,
                                       mpi.rank_in_pool);
    EXPECT_EQ(mpi.nproc_in_pool, 4);
    EXPECT_EQ(mpi.my_pool, 1);
    EXPECT_EQ(mpi.rank_in_pool, 1);
}


class FakeMPIContext
{
  public:
    FakeMPIContext()
    {
        _rank = 0;
        _size = 1;
    }

    int GetRank() const
    {
        return _rank;
    }
    int GetSize() const
    {
        return _size;
    }

    int drank;
    int dsize;
    int dcolor;

    int grank;
    int gsize;

    int kpar;
    int nproc_in_pool;
    int my_pool;
    int rank_in_pool;

    int nstogroup;
    int MY_BNDGROUP;
    int rank_in_stogroup;
    int nproc_in_stogroup;

  private:
    int _rank;
    int _size;
};

// --- DeathTest: Single thread ---
// Since these precondition checks cause the processes to die, we call such tests death tests.
// convention of naming the test suite: *DeathTest
// Death tests should be run in a single-threaded context.
// Such DeathTest will be run before all other tests.
class ParaGlobalDeathTest : public ::testing::Test
{
  protected:
    FakeMPIContext mpi;
    int nproc;
    int my_rank;
    int real_rank;

    // DeathTest SetUp:
    // Init variable, single thread
    void SetUp() override
    {
        int is_init = 0;
        MPI_Initialized(&is_init);
        if (is_init) {
             MPI_Comm_rank(MPI_COMM_WORLD, &real_rank);
        } else {
             real_rank = 0;
        }

        if (real_rank != 0) return;

        nproc = mpi.GetSize();
        my_rank = mpi.GetRank();

        // init log file needed by WARNING_QUIT
        GlobalV::ofs_warning.open("warning.log");


    }

    // clean log file
    void TearDown() override
    {
        if (real_rank != 0) return;

        GlobalV::ofs_warning.close();
        remove("warning.log");
    }
};

TEST_F(ParaGlobalDeathTest, InitPools)
{
    if (real_rank != 0) return;
    nproc = 12;
    mpi.kpar = 3;
    mpi.nstogroup = 3;
    my_rank = 5;
    EXPECT_EXIT(
    // This gtest Macro expect that a given `statement` causes the program to exit, with an
    // integer exit status that satisfies `predicate`(Here ::testing::ExitedWithCode(1)),
    // and emitting error output that matches `matcher`(Here "Error").
        {
            // redirect stdout to stderr to capture WARNING_QUIT output
            dup2(STDERR_FILENO, STDOUT_FILENO);
            Parallel_Global::init_pools(nproc,
                                my_rank,
                                mpi.nstogroup,
                                mpi.kpar,
                                mpi.nproc_in_stogroup,
                                mpi.rank_in_stogroup,
                                mpi.MY_BNDGROUP,
                                mpi.nproc_in_pool,
                                mpi.rank_in_pool,
                                mpi.my_pool);
        },
        ::testing::ExitedWithCode(1),
        "Error");
}

TEST_F(ParaGlobalDeathTest, DivideMPIPoolsNgEqZero)
{
    if (real_rank != 0) return;
    // test for num_groups == 0,
    // Num_group Equals 0
    // WARNING_QUIT
    this->nproc = 12;
    mpi.kpar = 0;
    EXPECT_EXIT(
        {
            // redirect stdout to stderr to capture WARNING_QUIT output
            dup2(STDERR_FILENO, STDOUT_FILENO);
            Parallel_Global::divide_mpi_groups(this->nproc,
                                       mpi.kpar,
                                       this->my_rank,
                                       mpi.nproc_in_pool,
                                       mpi.my_pool,
                                       mpi.rank_in_pool);
        },
        ::testing::ExitedWithCode(1),
        "Number of groups must be greater than 0."
    );
}

TEST_F(ParaGlobalDeathTest, DivideMPIPoolsNgGtProc)
{
    if (real_rank != 0) return;
    // test for procs < num_groups
    // Num_group GreaterThan Processors
    // WARNING_QUIT
    this->nproc = 12;
    mpi.kpar = 24;
    this->my_rank = 5;
    EXPECT_EXIT(
        {
            // redirect stdout to stderr to capture WARNING_QUIT output
            dup2(STDERR_FILENO, STDOUT_FILENO);
            Parallel_Global::divide_mpi_groups(this->nproc,
                                        mpi.kpar,
                                        this->my_rank,
                                        mpi.nproc_in_pool,
                                        mpi.my_pool,
                                        mpi.rank_in_pool);
        },
        testing::ExitedWithCode(1),
        "Error: Number of processes.*must be greater than the number of groups"
    );
}

int main(int argc, char** argv)
{
    bool is_death_test_child = false;
    for (int i = 0; i < argc; ++i) {
        if (std::string(argv[i]).find("gtest_internal_run_death_test") != std::string::npos) {
            is_death_test_child = true;
            break;
        }
    }

    if (!is_death_test_child)
    {
        MPI_Init(&argc, &argv);
    }

    testing::InitGoogleTest(&argc, argv);
    testing::FLAGS_gtest_death_test_style = "threadsafe";
    int result = RUN_ALL_TESTS();

    if (!is_death_test_child) {
        MPI_Finalize();
    }
    return result;
}
#endif // __MPI
