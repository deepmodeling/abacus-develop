//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-10
//==========================================================
#include "driver.h"
#include "src_parallel/parallel_global.h"
#include <ctime>

int main(int argc, char **argv)
{
    Parallel_Global::read_mpi_parameters(argc,argv);
	// main program for doing electronic structure calculations
	//----------------------------------------------------------
	Driver DD;
	DD.init();

#ifdef __MPI
	MPI_Comm_free(&POOL_WORLD);
	MPI_Comm_free(&STO_WORLD);
	MPI_Comm_free(&PARAPW_WORLD);
	MPI_Comm_free(&GRID_WORLD);
	MPI_Comm_free(&DIAG_WORLD);
    MPI_Finalize();
#endif

    return 0;
}