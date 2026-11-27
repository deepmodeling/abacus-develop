//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-10
//==========================================================

#include "fftw3.h"
#include "source_base/parallel_global.h"
#include "source_io/module_parameter/parameter.h"
#include "source_io/parse_args.h"
#include "source_io/parse_command_line.h"
#include "source_main/driver.h"
#include "source_main/version.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include <ctime>
#include <exception>
#include <iostream>
#include <new>
#include <stdexcept>

namespace
{
void print_welcome_banner()
{
#ifdef VERSION
    const char* version = VERSION;
#else
    const char* version = "unknown";
#endif
#ifdef COMMIT_INFO
#include "commit.h"
    const char* commit = COMMIT;
#else
    const char* commit = "unknown";
#endif
    std::cout << std::endl
              << "                     ▄████▄ █████▄ ▄████▄ ▄█████ ██  ██ ▄█████ " << std::endl
              << "                     ██▄▄██ ██▄▄██ ██▄▄██ ██     ██  ██ ▀▀▀▄▄▄ " << std::endl
              << "                     ██  ██ ██▄▄█▀ ██  ██ ▀█████ ▀████▀ █████▀ " << std::endl
              << std::endl;
    std::cout << "                                                                                     " << std::endl
              << "                                " << version << std::endl
              << std::endl
              << "               Atomic-orbital Based Ab-initio Computation at UStc                    " << std::endl
              << std::endl
              << "                     Website: http://abacus.ustc.edu.cn/                             " << std::endl
              << "               Documentation: https://abacus.deepmodeling.com/                       " << std::endl
              << "                  Repository: https://github.com/abacusmodeling/abacus-develop       " << std::endl
              << "                              https://github.com/deepmodeling/abacus-develop         " << std::endl
              << "                      Commit: " << commit << std::endl
              << std::endl;
    time_t time_now = time(nullptr);
    std::cout << " " << ctime(&time_now);
}
} // namespace

int main(int argc, char** argv)
{
    // Finish all CLI validation before acquiring MPI or FFTW resources.
    // Informational flags keep their existing early-exit precedence.
    ModuleIO::CommandLineArgs cli;
    try
    {
        ModuleIO::parse_args(argc, argv);
        cli = ModuleIO::parse_command_line(argc, argv);
    }
    catch (const std::exception& error)
    {
        std::cerr << "Command line error: " << error.what() << std::endl;
        ModuleIO::print_help(argv[0]);
        return 1;
    }

    /*
    read the mpi parameters in the command-line,
    initialize the mpi environment.
    */
    int nproc = 1;
    int my_rank = 0;
    int nthread_per_proc = 1;
    Parallel_Global::read_pal_param(argc, argv, nproc, nthread_per_proc, my_rank);
    if (my_rank == 0)
    {
        print_welcome_banner();
    }
#ifdef _OPENMP
    // ref: https://www.fftw.org/fftw3_doc/Usage-of-Multi_002dthreaded-FFTW.html
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
    PARAM.set_pal_param(my_rank, nproc, nthread_per_proc);

    if (my_rank == 0 && !cli.vars.empty())
    {
        std::cerr << "Warning: -p/--parameter variables are parsed but are not applied to INPUT in PART I."
                  << std::endl;
    }

    /*
    main program for doing electronic structure calculations.
    */
    Driver DD;
    // Catch allocation failures from the large LCAO dense matrices (issue
    // #8342) so the run exits with a clear message instead of an uncaught
    // exception or an MPI abort.  Driver stays in this scope so its destructor
    // still runs after finalize_mpi, preserving the original cleanup order.
    try
    {
        DD.init(cli);
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << " ABACUS: out of memory: " << e.what() << std::endl;
#ifdef __MPI
        // The allocation failure may be rank-local (e.g. a root-only gathered
        // buffer), so peers can be blocked inside a collective.  finalize_mpi
        // runs collective cleanup and would hang; abort the whole job instead.
        Parallel_Global::abort_mpi(1);
#endif
        return 1;
    }
    catch (const std::length_error& e)
    {
        std::cerr << " ABACUS: allocation length error (likely a too-large dense"
                     " matrix): " << e.what() << std::endl;
#ifdef __MPI
        Parallel_Global::abort_mpi(1);
#endif
        return 1;
    }
    catch (const std::exception& e)
    {
        std::cerr << " ABACUS: fatal error: " << e.what() << std::endl;
#ifdef __MPI
        Parallel_Global::abort_mpi(1);
#endif
        return 1;
    }

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
