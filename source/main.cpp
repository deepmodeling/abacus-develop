//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-10
//==========================================================

#include "driver.h"
#include "module_base/parallel_global.h"
#include "module_io/parse_args.h"

int main(int argc, char** argv)
{
    /*
    read the arguement in the command-line
    with "abacus -v", the program exit and returns version info
    with no arguments, the program continues.
    */
    ModuleIO::parse_args(argc, argv);

    /*
    read the mpi parameters in the command-line
    */
    Parallel_Global::read_mpi_parameters(argc, argv);

    /*
    main program for doing electronic structure calculations
    */
    Driver DD;
    DD.init();

    /*
    After running mpi version of abacus, release the mpi resources.
    */
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif

    return 0;
}