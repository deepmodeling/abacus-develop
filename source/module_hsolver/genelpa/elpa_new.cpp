#include "elpa_new.h"

#include "elpa_solver.h"
#include "my_math.hpp"
#include "utils.h"
#include <cfloat>
#include <complex>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <mpi.h>
#include <regex>
#include <sstream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

std::map<int, elpa_t> NEW_ELPA_HANDLE_POOL;

ELPA_Solver::ELPA_Solver(const bool isReal,
                         const MPI_Comm comm,
                         const int nev,
                         const int narows,
                         const int nacols,
                         const int* desc,
                         const bool reuse_handle_0)
{
    this->isReal = isReal;
    this->comm = comm;
    this->nev = nev;
    this->narows = narows;
    this->nacols = nacols;
    for (int i = 0; i < 9; ++i)
        this->desc[i] = desc[i];
    cblacs_ctxt = desc[1];
    nFull = desc[2];
    nblk = desc[4];
    lda = desc[8];
    // cout<<"parameters are passed\n";
    MPI_Comm_rank(comm, &myid);
    Cblacs_gridinfo(cblacs_ctxt, &nprows, &npcols, &myprow, &mypcol);
    // cout<<"blacs grid is inited\n";
    allocate_work();
    // cout<<"work array is inited\n";
    if (isReal)
        kernel_id = read_real_kernel();
    else
        kernel_id = read_complex_kernel();
    // cout<<"kernel id is inited as "<<kernel_id<<"\n";
    int error;

    static int total_handle = 0;

#ifdef _OPENMP
    int num_threads = omp_get_max_threads();
#else
    int num_threads = 1;
#endif

    elpa_init(20210430);

    handle_id = ++total_handle;

    //delete the old elpa_handle and reuse the handle_id=0
    if(reuse_handle_0 && total_handle>0)
    {
        NEW_ELPA_HANDLE_POOL.erase(0);
        handle_id = 0;
    }

    elpa_t handle;

    handle = elpa_allocate(&error);
    NEW_ELPA_HANDLE_POOL[handle_id] = handle;

#ifdef _OPENMP
    elpa_set_integer(NEW_ELPA_HANDLE_POOL[handle_id], "omp_threads", num_threads, &error);
#endif
    elpa_set_integer(NEW_ELPA_HANDLE_POOL[handle_id], "na", nFull, &error);
    elpa_set_integer(NEW_ELPA_HANDLE_POOL[handle_id], "nev", nev, &error);
    elpa_set_integer(NEW_ELPA_HANDLE_POOL[handle_id], "local_nrows", narows, &error);
    elpa_set_integer(NEW_ELPA_HANDLE_POOL[handle_id], "local_ncols", nacols, &error);
    elpa_set_integer(NEW_ELPA_HANDLE_POOL[handle_id], "nblk", nblk, &error);
    elpa_set_integer(NEW_ELPA_HANDLE_POOL[handle_id], "mpi_comm_parent", MPI_Comm_c2f(comm), &error);
    elpa_set_integer(NEW_ELPA_HANDLE_POOL[handle_id], "process_row", myprow, &error);
    elpa_set_integer(NEW_ELPA_HANDLE_POOL[handle_id], "process_col", mypcol, &error);

    error = elpa_setup(NEW_ELPA_HANDLE_POOL[handle_id]);
    // cout<<"elpa handle is setup\n";
    elpa_set_integer(NEW_ELPA_HANDLE_POOL[handle_id], "solver", ELPA_SOLVER_2STAGE, &error);
    this->setQR(0);
    this->setKernel(isReal, kernel_id);
    // cout<<"elpa kernel is setup\n";
    this->setLoglevel(0);
    // cout<<"log level is setup\n";
}

ELPA_Solver::ELPA_Solver(const bool isReal,
                         const MPI_Comm comm,
                         const int nev,
                         const int narows,
                         const int nacols,
                         const int* desc,
                         const int* otherParameter)
{
    this->isReal = isReal;
    this->comm = comm;
    this->nev = nev;
    this->narows = narows;
    this->nacols = nacols;
    for (int i = 0; i < 9; ++i)
        this->desc[i] = desc[i];

    kernel_id = otherParameter[0];
    useQR = otherParameter[1];
    loglevel = otherParameter[2];

    cblacs_ctxt = desc[1];
    nFull = desc[2];
    nblk = desc[4];
    lda = desc[8];
    MPI_Comm_rank(comm, &myid);
    Cblacs_gridinfo(cblacs_ctxt, &nprows, &npcols, &myprow, &mypcol);
    allocate_work();

    int error;
    static std::map<int, elpa_t> NEW_ELPA_HANDLE_POOL;
    static int total_handle;

#ifdef _OPENMP
    int num_threads = omp_get_max_threads();
#else
    int num_threads = 1;
#endif

    elpa_init(20210430);

    handle_id = ++total_handle;
    elpa_t handle;
    handle = elpa_allocate(&error);
    NEW_ELPA_HANDLE_POOL[handle_id] = handle;

#ifdef _OPENMP
    elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "omp_threads", num_threads, &error);
#endif
    elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "na", nFull, &error);
    elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "nev", nev, &error);
    elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "local_nrows", narows, &error);
    elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "local_ncols", nacols, &error);
    elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "nblk", nblk, &error);
    elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "mpi_comm_parent", MPI_Comm_c2f(comm), &error);
    elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "process_row", myprow, &error);
    elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "process_col", mypcol, &error);
    elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "blacs_context", cblacs_ctxt, &error);
    elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "solver", ELPA_SOLVER_2STAGE, &error);
    elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "debug", wantDebug, &error);
    elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "qr", useQR, &error);
    this->setQR(useQR);
    this->setKernel(isReal, kernel_id);
    this->setLoglevel(loglevel);
}

void ELPA_Solver::setLoglevel(int loglevel)
{
    int error;
    this->loglevel = loglevel;
    static bool isLogfileInited = false;

    if (loglevel >= 2)
    {
        wantDebug = 1;
        elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "verbose", 1, &error);
        elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "debug", wantDebug, &error);
        if (!isLogfileInited)
        {
            std::stringstream logfilename;
            logfilename.str("");
            logfilename << "GenELPA_" << myid << ".log";
            logfile.open(logfilename.str());
            logfile << "logfile inited\n";
            isLogfileInited = true;
        }
    }
    else
    {
        wantDebug = 0;
    }
}

void ELPA_Solver::setKernel(bool isReal, int kernel)
{
    this->kernel_id = kernel;
    int error;
    if (isReal)
        elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "real_kernel", kernel, &error);
    else
        elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "complex_kernel", kernel, &error);
}

void ELPA_Solver::setQR(int useQR)
{
    this->useQR = useQR;
    int error;
    elpa_set(NEW_ELPA_HANDLE_POOL[handle_id], "qr", useQR, &error);
}

void ELPA_Solver::exit()
{
    // delete[] dwork;
    // delete[] zwork;
    if (loglevel > 2)
        logfile.close();
    int error;
    elpa_deallocate(NEW_ELPA_HANDLE_POOL[handle_id], &error);
    elpa_uninit(&error);
}

int ELPA_Solver::read_cpuflag()
{
    int cpuflag = 0;

    std::ifstream f_cpuinfo("/proc/cpuinfo");
    std::string cpuinfo_line;
    std::regex cpuflag_ex("flags.*");
    std::regex cpuflag_avx512(".*avx512.*");
    std::regex cpuflag_avx2(".*avx2.*");
    std::regex cpuflag_avx(".*avx.*");
    std::regex cpuflag_sse(".*sse.*");
    while (getline(f_cpuinfo, cpuinfo_line))
    {
        if (std::regex_match(cpuinfo_line, cpuflag_ex))
        {
            // cout<<cpuinfo_line<<endl;
            if (std::regex_match(cpuinfo_line, cpuflag_avx512))
            {
                cpuflag = 4;
            }
            else if (std::regex_match(cpuinfo_line, cpuflag_avx2))
            {
                cpuflag = 3;
            }
            else if (std::regex_match(cpuinfo_line, cpuflag_avx))
            {
                cpuflag = 2;
            }
            else if (std::regex_match(cpuinfo_line, cpuflag_sse))
            {
                cpuflag = 1;
            }
            break;
        }
    }
    f_cpuinfo.close();
    return cpuflag;
}

int ELPA_Solver::read_real_kernel()
{
    int kernel_id;

    if (const char* env = getenv("ELPA_DEFAULT_real_kernel"))
    {
        if (strcmp(env, "ELPA_2STAGE_REAL_GENERIC_SIMPLE") == 0)
            kernel_id = ELPA_2STAGE_REAL_GENERIC_SIMPLE;
        else if (strcmp(env, "ELPA_2STAGE_REAL_BGP") == 0)
            kernel_id = ELPA_2STAGE_REAL_BGP;
        else if (strcmp(env, "ELPA_2STAGE_REAL_BGQ") == 0)
            kernel_id = ELPA_2STAGE_REAL_BGQ;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SSE_ASSEMBLY") == 0)
            kernel_id = ELPA_2STAGE_REAL_SSE_ASSEMBLY;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SSE_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_REAL_SSE_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SSE_BLOCK4") == 0)
            kernel_id = ELPA_2STAGE_REAL_SSE_BLOCK4;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SSE_BLOCK6") == 0)
            kernel_id = ELPA_2STAGE_REAL_SSE_BLOCK6;
        else if (strcmp(env, "ELPA_2STAGE_REAL_AVX_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_REAL_AVX_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_REAL_AVX_BLOCK4") == 0)
            kernel_id = ELPA_2STAGE_REAL_AVX_BLOCK4;
        else if (strcmp(env, "ELPA_2STAGE_REAL_AVX_BLOCK6") == 0)
            kernel_id = ELPA_2STAGE_REAL_AVX_BLOCK6;
        else if (strcmp(env, "ELPA_2STAGE_REAL_AVX2_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_REAL_AVX2_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_REAL_AVX2_BLOCK4") == 0)
            kernel_id = ELPA_2STAGE_REAL_AVX2_BLOCK4;
        else if (strcmp(env, "ELPA_2STAGE_REAL_AVX2_BLOCK6") == 0)
            kernel_id = ELPA_2STAGE_REAL_AVX2_BLOCK6;
        else if (strcmp(env, "ELPA_2STAGE_REAL_AVX512_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_REAL_AVX512_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_REAL_AVX512_BLOCK4") == 0)
            kernel_id = ELPA_2STAGE_REAL_AVX512_BLOCK4;
        else if (strcmp(env, "ELPA_2STAGE_REAL_AVX512_BLOCK6") == 0)
            kernel_id = ELPA_2STAGE_REAL_AVX512_BLOCK6;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SPARC64_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_REAL_SPARC64_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SPARC64_BLOCK4") == 0)
            kernel_id = ELPA_2STAGE_REAL_SPARC64_BLOCK4;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SPARC64_BLOCK6") == 0)
            kernel_id = ELPA_2STAGE_REAL_SPARC64_BLOCK6;
        else if (strcmp(env, "ELPA_2STAGE_REAL_NEON_ARCH64_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_REAL_NEON_ARCH64_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_REAL_NEON_ARCH64_BLOCK4") == 0)
            kernel_id = ELPA_2STAGE_REAL_NEON_ARCH64_BLOCK4;
        else if (strcmp(env, "ELPA_2STAGE_REAL_NEON_ARCH64_BLOCK6") == 0)
            kernel_id = ELPA_2STAGE_REAL_NEON_ARCH64_BLOCK6;
        else if (strcmp(env, "ELPA_2STAGE_REAL_VSX_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_REAL_VSX_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_REAL_VSX_BLOCK4") == 0)
            kernel_id = ELPA_2STAGE_REAL_VSX_BLOCK4;
        else if (strcmp(env, "ELPA_2STAGE_REAL_VSX_BLOCK6") == 0)
            kernel_id = ELPA_2STAGE_REAL_VSX_BLOCK6;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SVE128_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_REAL_SVE128_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SVE128_BLOCK4") == 0)
            kernel_id = ELPA_2STAGE_REAL_SVE128_BLOCK4;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SVE128_BLOCK6") == 0)
            kernel_id = ELPA_2STAGE_REAL_SVE128_BLOCK6;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SVE256_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_REAL_SVE256_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SVE256_BLOCK4") == 0)
            kernel_id = ELPA_2STAGE_REAL_SVE256_BLOCK4;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SVE256_BLOCK6") == 0)
            kernel_id = ELPA_2STAGE_REAL_SVE256_BLOCK6;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SVE512_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_REAL_SVE512_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SVE512_BLOCK4") == 0)
            kernel_id = ELPA_2STAGE_REAL_SVE512_BLOCK4;
        else if (strcmp(env, "ELPA_2STAGE_REAL_SVE512_BLOCK6") == 0)
            kernel_id = ELPA_2STAGE_REAL_SVE512_BLOCK6;
        else if (strcmp(env, "ELPA_2STAGE_REAL_GENERIC_SIMPLE_BLOCK4") == 0)
            kernel_id = ELPA_2STAGE_REAL_GENERIC_SIMPLE_BLOCK4;
        else if (strcmp(env, "ELPA_2STAGE_REAL_GENERIC_SIMPLE_BLOCK6") == 0)
            kernel_id = ELPA_2STAGE_REAL_GENERIC_SIMPLE_BLOCK6;
        else
            kernel_id = ELPA_2STAGE_REAL_GENERIC;
    }
    else
    {
        int cpuflag = read_cpuflag();
        switch (cpuflag)
        {
        case 4:
            kernel_id = ELPA_2STAGE_REAL_AVX512_BLOCK4;
            break;
        case 3:
            kernel_id = ELPA_2STAGE_REAL_AVX2_BLOCK2;
            break;
        case 2:
            kernel_id = ELPA_2STAGE_REAL_AVX_BLOCK2;
            break;
        case 1:
            kernel_id = ELPA_2STAGE_REAL_SSE_BLOCK2;
            break;
        default:
            kernel_id = ELPA_2STAGE_REAL_GENERIC_SIMPLE_BLOCK6;
            break;
        }
    }
    return kernel_id;
}

int ELPA_Solver::read_complex_kernel()
{
    int kernel_id;
    if (const char* env = getenv("ELPA_DEFAULT_complex_kernel"))
    {
        if (strcmp(env, "ELPA_2STAGE_COMPLEX_GENERIC_SIMPLE") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_GENERIC_SIMPLE;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_BGP") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_BGP;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_BGQ") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_BGQ;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_SSE_ASSEMBLY") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_SSE_ASSEMBLY;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_SSE_BLOCK1") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_SSE_BLOCK1;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_SSE_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_SSE_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_AVX_BLOCK1") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_AVX_BLOCK1;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_AVX_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_AVX_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_AVX2_BLOCK1") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_AVX2_BLOCK1;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_AVX2_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_AVX2_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_AVX512_BLOCK1") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_AVX512_BLOCK1;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_AVX512_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_AVX512_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_SVE128_BLOCK1") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_SVE128_BLOCK1;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_SVE128_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_SVE128_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_SVE256_BLOCK1") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_SVE256_BLOCK1;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_SVE256_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_SVE256_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_SVE512_BLOCK1") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_SVE512_BLOCK1;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_SVE512_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_SVE512_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_NEON_ARCH64_BLOCK1") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_NEON_ARCH64_BLOCK1;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_NEON_ARCH64_BLOCK2") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_NEON_ARCH64_BLOCK2;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_NVIDIA_GPU") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_NVIDIA_GPU;
        else if (strcmp(env, "ELPA_2STAGE_COMPLEX_AMD_GPU") == 0)
            kernel_id = ELPA_2STAGE_COMPLEX_AMD_GPU;
        else
            kernel_id = ELPA_2STAGE_COMPLEX_GENERIC;
    }
    else
    {
        int cpuflag = read_cpuflag();
        switch (cpuflag)
        {
        case 4:
            kernel_id = ELPA_2STAGE_COMPLEX_AVX512_BLOCK2;
            break;
        case 3:
            kernel_id = ELPA_2STAGE_COMPLEX_AVX2_BLOCK2;
            break;
        case 2:
            kernel_id = ELPA_2STAGE_COMPLEX_AVX_BLOCK2;
            break;
        case 1:
            kernel_id = ELPA_2STAGE_COMPLEX_SSE_BLOCK2;
            break;
        default:
            kernel_id = ELPA_2STAGE_COMPLEX_GENERIC_SIMPLE;
            break;
        }
    }
    return kernel_id;
}

int ELPA_Solver::allocate_work()
{
    unsigned long nloc = static_cast<unsigned long>(narows) * nacols; // local size
    unsigned long maxloc; // maximum local size
    MPI_Allreduce(&nloc, &maxloc, 1, MPI_UNSIGNED_LONG, MPI_MAX, comm);
    maxloc = nloc;

    if (isReal)
        dwork.resize(maxloc);
    else
        zwork.resize(maxloc);
    return 0;
}

void ELPA_Solver::timer(int myid, const char function[], const char step[], double& t0)
{
    double t1;
    if (t0 < 0) // t0 < 0 means this is the init call before the function
    {
        t0 = MPI_Wtime();
        t0 = (double)clock()/CLOCKS_PER_SEC;
        logfile << "DEBUG: Process " << myid << " Call " << function << std::endl;
    }
    else
    {
        t1 = MPI_Wtime();
        t1 = (double)clock()/CLOCKS_PER_SEC;
        logfile << "DEBUG: Process " << myid << " Step " << step << " " << function << " time: " << t1 - t0 << " s"
                << std::endl;
    }
}

void ELPA_Solver::outputParameters()
{
    logfile << "myid " << myid << ": comm id(in FORTRAN):" << MPI_Comm_c2f(comm) << std::endl;
    logfile << "myid " << myid << ": nprows: " << nprows << " npcols: " << npcols << std::endl;
    logfile << "myid " << myid << ": myprow: " << myprow << " mypcol: " << mypcol << std::endl;
    logfile << "myid " << myid << ": nFull: " << nFull << " nev: " << nev << std::endl;
    logfile << "myid " << myid << ": narows: " << narows << " nacols: " << nacols << std::endl;
    logfile << "myid " << myid << ": blacs parameters setting" << std::endl;
    logfile << "myid " << myid << ": blacs ctxt:" << cblacs_ctxt << std::endl;
    logfile << "myid " << myid << ": desc: ";
    for (int i = 0; i < 9; ++i)
        logfile << desc[i] << " ";
    logfile << std::endl;
    logfile << "myid " << myid << ": nblk: " << nblk << " lda: " << lda << std::endl;
    logfile << "myid " << myid << ": useQR: " << useQR << " kernel:" << kernel_id << std::endl;
    ;
    logfile << "myid " << myid << ": wantDebug: " << wantDebug << " loglevel: " << loglevel << std::endl;
}
