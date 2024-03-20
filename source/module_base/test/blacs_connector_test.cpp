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

TEST(blacs_connector, Cblacs_gridinit)
{
    int icontxt;
    char layout[] = "ROW";
    int nprow = 1;
    int npcol = 1;

    int myid, nprocs;
    Cblacs_pinfo(&myid, &nprocs);
    Cblacs_get(-1, 0, &icontxt);

    // Call the Cblacs_gridinit() function
    Cblacs_gridinit(&icontxt, layout, nprow, npcol);

    // Check if the grid context handle is created successfully
    EXPECT_EQ(icontxt, 0);
}

int main(int argc, char** argv)
{
    int myrank;
    int mysize;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mysize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    testing::InitGoogleTest(&argc, argv);

    int result = 0;
    result = RUN_ALL_TESTS();
    MPI_Finalize();
    return 0;
}
