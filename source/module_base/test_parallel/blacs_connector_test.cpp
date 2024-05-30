#ifdef __MPI

#include "../blacs_connector.h"
#include <mpi.h>
#include "gtest/gtest.h"

/************************************************
 *  unit test of functions in blacs_connector.h
 ***********************************************/

/**
 * - Tested Function
 *   - Cblacs_gridinit
 *     - Initializes a grid of processors with a given number of rows and columns.
 *        The function creates a cartesian topology of all the processors initialized
 *        by the BLS library. In this topology, each processor is identified by its
 *        coordinates (row, col) in the grid.
 */

class BLACSTest: public testing::Test
{
protected:
    void SetUp();

    int rank = 0;
    int nprocs = 0;
    char layout = 'R';

    // number of rows and columns in the process grid
    int nprow = 0;
    int npcol = 0;

    // process coordinate
    int iprow = -1;
    int ipcol = -1;
};

void BLACSTest::SetUp()
{
    Cblacs_pinfo(&rank, &nprocs);
}


TEST_F(BLACSTest, WorldGrid)
{
    // generate a grid of size 1 x nproc
    nprow = 1;
    npcol = nprocs;

    int ictxt_row = Csys2blacs_handle(MPI_COMM_WORLD);
    Cblacs_gridinit(&ictxt_row, &layout, nprow, npcol);
    Cblacs_gridinfo(ictxt_row, &nprow, &npcol, &iprow, &ipcol);

    EXPECT_EQ(iprow, 0);
    EXPECT_EQ(ipcol, rank);

    // generate a grid of size nproc x 1
    nprow = nprocs;
    npcol = 1;

    int ictxt_col = Csys2blacs_handle(MPI_COMM_WORLD);
    Cblacs_gridinit(&ictxt_col, &layout, nprow, npcol);
    Cblacs_gridinfo(ictxt_col, &nprow, &npcol, &iprow, &ipcol);

    EXPECT_EQ(iprow, rank);
    EXPECT_EQ(ipcol, 0);


    // two BLACS grids should have difference context index
    EXPECT_NE(ictxt_row, ictxt_col);
}

TEST_F(BLACSTest, SplitGrid)
{
    // this test create BLACS grids based on a disjoint communicator

    const int n_blacs = 2;
    int rank_sub = -1;
    int nprocs_sub = 0;

    // sub communicators are divided based on odd / even ranks
    MPI_Comm comm_sub;
    MPI_Comm_split(MPI_COMM_WORLD, rank % n_blacs, rank, &comm_sub);
    MPI_Comm_rank(comm_sub, &rank_sub);
    MPI_Comm_size(comm_sub, &nprocs_sub);

    int ctxt_sub = Csys2blacs_handle(comm_sub);

    nprow = 1, npcol = nprocs_sub; // row-like grids
    Cblacs_gridinit(&ctxt_sub, &layout, nprow, npcol);
    Cblacs_gridinfo(ctxt_sub, &nprow, &npcol, &iprow, &ipcol);

    // verifies that the BLACS grid is created based on comm_sub instead of MPI_COMM_WORLD
    EXPECT_EQ(iprow, 0);
    EXPECT_EQ(ipcol, rank_sub);
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
