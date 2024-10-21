//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-10
//==========================================================

#include "driver.h"
#include "fftw3.h"
#include "module_base/parallel_global.h"
#include "module_io/parse_args.h"
#include "module_parameter/parameter.h"
#ifdef _OPENMP
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
    After running mpi version of abacus, release the mpi resources.
    */
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
#ifdef _OPENMP
    fftw_cleanup_threads();
#endif

    return 0;
}