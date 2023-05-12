//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-10
//==========================================================
#include "driver.h"
#include "module_base/parallel_global.h"
#include <ctime>
#include <iostream>

int main(int argc, char **argv)
{
    if (argc > 1 && std::string(argv[1]) == "--version")
    {
#ifdef VERSION
        const char* version = VERSION;
#else
        const char* version = "unknown";
#endif
        std::cout << "ABACUS version " << version << std::endl;
        std::exit(0);
    }

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