#include "tddft_test.h"

#include <gtest/gtest.h>
#include <mpi.h>

#include "module_base/blacs_connector.h"
#include "module_base/scalapack_connector.h"
#include "module_basis/module_ao/parallel_orbitals.h"

/************************************************
 *  unit test of module_tddft
 ***********************************************/

int myprow, nprow, ictxt, mypcol, npcol;

void MPIInit()
{
    int myrank;
    int mysize;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &mysize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Set up BLACS context and grid
    int nprocs;
    nprow = 1;
    npcol = 1;
    Cblacs_pinfo(&myrank, &mysize);
    Cblacs_get(-1, 0, &ictxt);
    char order[] = "Row";
    Cblacs_gridinit(&ictxt, order, nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myprow, &mypcol);
}

/************************************************
 *  unit test of module_tddft
 ***********************************************/
int main(int argc, char** argv)
{
    MPIInit();

    int result = 0;
    testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();

    Cblacs_exit(ictxt);

    // MPI_Finalize();
    return 0;
}
