#pragma once

#include "hamilt_bse_solver.h"
namespace BSE
{
template <typename T>
void solve_tda(const int my_rank,
               const std::vector<T>& A,
               const Parallel_2D& pA,
               std::vector<double>& ev,
               std::vector<T>& v)
{
    ModuleBase::TITLE("HamiltBSE", "elpa_solve_tda1");
    ModuleBase::timer::start("HamiltBSE", "elpa_solve_tda");

    assert(pA.get_global_row_size() == pA.get_global_col_size());
    const int nA = pA.get_global_row_size();

    elpa_t elpaInstance;
    int status;

    if (elpa_init(20210430) != ELPA_OK)
    {
        fprintf(stderr, "Error: ELPA API version not supported");
        exit(1);
    }

    elpaInstance = elpa_allocate(&status);
    if (status != ELPA_OK)
    {
        std::cout << "Could not allocate elpa instance" << std::endl;
        exit(1);
    }
    elpa_set(elpaInstance, "na", nA, &status);
    elpa_set(elpaInstance, "nev", static_cast<int>(ev.size()), &status);
    elpa_set(elpaInstance, "local_nrows", pA.get_row_size(), &status);
    elpa_set(elpaInstance, "local_ncols", pA.get_col_size(), &status);
    elpa_set(elpaInstance, "nblk", pA.get_block_size(), &status);
    elpa_set(elpaInstance, "mpi_comm_parent", MPI_Comm_c2f(MPI_COMM_WORLD), &status);
    elpa_set(elpaInstance, "process_row", pA.get_coord_row(), &status);
    elpa_set(elpaInstance, "process_col", pA.get_coord_col(), &status);
    elpa_set(elpaInstance, "solver", ELPA_SOLVER_2STAGE, &status);
#ifdef _OPENMP
    int num_threads = omp_get_max_threads();
#else
    int num_threads = 1;
#endif
    elpa_set(elpaInstance, "omp_threads", num_threads, &status);
    status = elpa_setup(elpaInstance);
    if (status != ELPA_OK)
    {
        fprintf(stderr, "Could not set up the ELPA object");
    }
    ModuleBase::TITLE("HamiltBSE", "elpa_solve_tda2");
    elpa_eigenvectors(elpaInstance, const_cast<T*>(A.data()), ev.data(), v.data(), &status);
    ModuleBase::TITLE("HamiltBSE", "elpa_solve_tda3");
    elpa_deallocate(elpaInstance, &status);
    elpa_uninit(&status);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "elpa_solve_tda");
    ModuleBase::timer::end("HamiltBSE", "elpa_solve_tda");
}
}