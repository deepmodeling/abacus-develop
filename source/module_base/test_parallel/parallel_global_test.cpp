#ifdef __MPI
#include "module_base/parallel_global.h"

#include "mpi.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <complex>
#include <cstring>
#include <string>

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

TEST_F(ParaGlobal, InitPools)
{
    nproc = 12;
    mpi.kpar = 3;
    mpi.nstogroup = 3;
    my_rank = 5;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(Parallel_Global::init_pools(nproc,
                                my_rank,
                                mpi.nstogroup,
                                mpi.kpar,
                                mpi.nproc_in_stogroup,
                                mpi.rank_in_stogroup,
                                mpi.MY_BNDGROUP,
                                mpi.nproc_in_pool,
                                mpi.rank_in_pool,
                                mpi.my_pool), ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Error:"));
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

int main(int argc, char** argv)
{

    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
#endif
