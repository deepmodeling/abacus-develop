//==========================================================
// AUTHOR : liuyu
// DATE : 2021-07-15
//==========================================================
#include "driver_classic.h"
#include "../src_parallel/parallel_global.h"
#include "../module_base/timer.h"
#include <ctime>

int main(int argc, char **argv)
{
    srand(time(0));
    Parallel_Global::read_mpi_parameters(argc,argv);

    //----------------------------------------------------------
	// main program for doing CMD calculations
	//----------------------------------------------------------
    Driver_classic::init();

#ifdef __MPI
    MPI_Finalize();
#endif

    return 0;
}