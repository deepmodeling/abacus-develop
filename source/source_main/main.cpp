//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-10
//==========================================================

#include <unistd.h>
#include "source_main/driver.h"
#include "source_base/parallel_global.h"
#include "source_io/parse_args.h"
#include "source_io/module_parameter/parameter.h"
#ifdef _OPENMP
#include <fftw3.h>
#include <omp.h>
#endif

int main(int argc, char** argv)
{
    /*
    read the arguement in the command-line,
    with "abacus -v", the program exit and returns version info,
    with no arguments, the program continues.
    */
    ModuleIO::parse_args(argc, argv);

    /*
    read the mpi parameters in the command-line,
    initialize the mpi environment.
    */
    int nproc = 1;
    int my_rank = 0;
    int nthread_per_proc = 1;
    Parallel_Global::read_pal_param(argc, argv, nproc, nthread_per_proc, my_rank);
#ifdef _OPENMP
    // ref: https://www.fftw.org/fftw3_doc/Usage-of-Multi_002dthreaded-FFTW.html
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
    PARAM.set_pal_param(my_rank, nproc, nthread_per_proc);

    /*
    main program for doing electronic structure calculations.
    */
    Driver DD;
    DD.init();

    /*
    Skip MPI_Finalize to avoid OpenMPI 4.0.3 hwloc segfault.
    OS will reclaim all resources on _exit.
    For non-MPI builds, clean up FFTW threads normally.
    */
#ifdef __MPI
    _exit(0);
#endif
#ifdef _OPENMP
    fftw_cleanup_threads();
#endif

    return 0;
}
