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
	Parallel_Global::finalize_mpi();
#endif

    return 0;
}