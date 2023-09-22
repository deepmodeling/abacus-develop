//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-10
//==========================================================

#include <ctime>
#include "driver.h"
#include "module_base/parallel_global.h"
#include "module_io/parse_args.h"

int main(int argc, char** argv)
{
    ModuleIO::parse_args(argc, argv);

    Parallel_Global::read_mpi_parameters(argc, argv);
    // main program for doing electronic structure calculations
    //----------------------------------------------------------
    Driver DD;
    DD.init();

#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif

    return 0;
}