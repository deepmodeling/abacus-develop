#include "diago_elpa_native.h"

#include "module_base/blacs_connector.h"
#include "module_base/global_variable.h"
#include "module_parameter/parameter.h"
#include "module_base/lapack_connector.h"
#include "module_base/scalapack_connector.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_hsolver/genelpa/elpa_new.h"
#include "omp.h"

namespace hsolver
{
#ifdef __MPI
template <typename T>
MPI_Comm DiagoElpaNative<T>::setmpicomm()
{
    if (this->elpa_num_thread == -1)
    {
        return MPI_COMM_WORLD;
    }
    else
    {
        int _num = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &_num);
        if (elpa_num_thread > _num || elpa_num_thread <= 0)
        {
            return MPI_COMM_WORLD;
        }
        else
        {
            lastmpinum++;
            int* _ranks = new int[elpa_num_thread];
            for (int i = 0; i < elpa_num_thread; i++)
            {
                _ranks[i] = (lastmpinum + i) % _num;
            }
            MPI_Group _tempgroup, _oldgroup;
            MPI_Comm_group(MPI_COMM_WORLD, &_oldgroup);
            MPI_Group_incl(_oldgroup, elpa_num_thread, _ranks, &_tempgroup);
            MPI_Comm _new_comm;
            MPI_Comm_create(MPI_COMM_WORLD, _tempgroup, &_new_comm);
            delete[] _ranks;
            return _new_comm;
        }
    }
}
#endif

#ifdef __MPI
template <typename T>
void DiagoElpaNative<T>::diag_pool(hamilt::MatrixBlock<T>& h_mat,
                                   hamilt::MatrixBlock<T>& s_mat,
                                   psi::Psi<T>& psi,
                                   Real* eigenvalue_in,
                                   MPI_Comm& comm)
{

    ModuleBase::timer::tick("DiagoElpaNative", "elpa_solve");

    int nev = PARAM.inp.nbands;
    int narows = h_mat.row;
    int nacols = h_mat.col;

    int cblacs_ctxt = h_mat.desc[1];
    int nFull = h_mat.desc[2];
    int nblk = h_mat.desc[4];
    // cout<<"parameters are passed\n";
    int nprows, npcols, myprow, mypcol;

    Cblacs_gridinfo(cblacs_ctxt, &nprows, &npcols, &myprow, &mypcol);
    std::vector<Real> eigen(PARAM.globalv.nlocal, 0.0);
    std::vector<T> eigenvectors(narows * nacols);

    if (elpa_init(20210430) != ELPA_OK)
    {
        fprintf(stderr, "Error: ELPA API version not supported");
    }

    // elpa_init(20210430);
    int success;
    elpa_t handle = elpa_allocate(&success);
    int num_threads = omp_get_max_threads();
    elpa_set(handle, "omp_threads", num_threads, &success);
    elpa_set(handle, "na", (int)nFull, &success);
    elpa_set(handle, "nev", (int)nev, &success);
    elpa_set(handle, "local_nrows", (int)narows, &success);
    elpa_set(handle, "local_ncols", (int)nacols, &success);
    elpa_set(handle, "nblk", (int)nblk, &success);
    elpa_set(handle, "mpi_comm_parent", (int)(MPI_Comm_c2f(comm)), &success);
    elpa_set(handle, "process_row", (int)myprow, &success);
    elpa_set(handle, "process_col", (int)mypcol, &success);
    elpa_set(handle, "blacs_context", (int)cblacs_ctxt, &success);
    elpa_setup(handle);
    elpa_set(handle, "solver", ELPA_SOLVER_1STAGE, &success);

/*  ELPA_WITH_NVIDIA_GPU_VERSION is a symbol defined in elpa/elpa_configured_options.h
    For example:
    cat elpa/elpa_configured_options.h
    #define ELPA_WITH_NVIDIA_GPU_VERSION 1
    #define ELPA_WITH_AMD_GPU_VERSION 0
    #define ELPA_WITH_SYCL_GPU_VERSION 0
 */
#if ELPA_WITH_NVIDIA_GPU_VERSION
    if (PARAM.inp.device == "gpu")
    {
        elpa_set(handle, "nvidia-gpu", 1, &success);
        elpa_set(handle, "real_kernel", ELPA_2STAGE_REAL_NVIDIA_GPU, &success);
        elpa_setup_gpu(handle);
    }
#endif

    elpa_generalized_eigenvectors(handle,
                                  h_mat.p,
                                  s_mat.p,
                                  eigen.data(),
                                  eigenvectors.data(),
                                  this->DecomposedState,
                                  &success);
    elpa_deallocate(handle, &success);
    elpa_uninit(&success);

    ModuleBase::timer::tick("DiagoElpaNative", "elpa_solve");
    if (std::is_same<T, double>::value)
    {
        // for gamma only, the decomposed s_mat will be reused
        this->DecomposedState = 1;
    }
    else
    {
        // for k pointer, the decomposed s_mat can not be reused
        this->DecomposedState = 0;
    }

    const int inc = 1;
    BlasConnector::copy(PARAM.inp.nbands, eigen.data(), inc, eigenvalue_in, inc);
    const int size = psi.get_nbands() * psi.get_nbasis();
    BlasConnector::copy(size, eigenvectors.data(), inc, psi.get_pointer(), inc);
}
#endif

template <typename T>
void DiagoElpaNative<T>::diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoElpaNative", "diag");
#ifdef __MPI
    hamilt::MatrixBlock<T> h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);
    MPI_Comm COMM_DIAG = setmpicomm(); // set mpi_comm needed
    diag_pool(h_mat, s_mat, psi, eigenvalue_in, COMM_DIAG);
#else
    ModuleBase::WARNING_QUIT("DiagoElpaNative", "DiagoElpaNative only can be used with macro __MPI");
#endif
}

template class DiagoElpaNative<double>;
template class DiagoElpaNative<std::complex<double>>;

} // namespace hsolver
