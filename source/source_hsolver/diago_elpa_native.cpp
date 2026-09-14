#include "diago_elpa_native.h"

#include "source_base/global_function.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/blacs_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_hsolver/module_genelpa/elpa_new.h"
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
void DiagoElpaNative<T>::diag_pool(ModuleBase::MatrixBlock<T>& h_mat,
                                   ModuleBase::MatrixBlock<T>& s_mat,
                                   psi::Psi<T>& psi,
                                   Real* eigenvalue_in,
                                   MPI_Comm& comm)
{

    ModuleBase::timer::start("DiagoElpaNative", "elpa_solve");

    int nev = this->nbands;
    int narows = h_mat.row;
    int nacols = h_mat.col;

    int cblacs_ctxt = h_mat.desc[1];
    int nFull = h_mat.desc[2];
    int nblk = h_mat.desc[4];
    // cout<<"parameters are passed\n";
    int nprows, npcols, myprow, mypcol;

    Cblacs_gridinfo(cblacs_ctxt, &nprows, &npcols, &myprow, &mypcol);
    std::vector<Real> eigen(this->nlocal, 0.0);
    std::vector<T> eigenvectors(narows * nacols);

    // The complex LCAO matrices follow the LAPACK UPLO='U' convention: only
    // their upper triangles are guaranteed to contain the Hermitian matrix.
    // ELPA's native generalized solver consumes both triangles, so complete
    // private Hermitian copies before the solve.
    std::vector<T> h_work;
    std::vector<T> s_work;
    T* h_elpa = h_mat.p;
    T* s_elpa = s_mat.p;
    int decomposed_state = this->DecomposedState;
    if (!std::is_same<T, double>::value)
    {
        h_work.resize(narows * nacols);
        s_work.resize(narows * nacols);
        const int one = 1;
        ScalapackConnector::tranc(nFull,
                                  nFull,
                                  T(1.0),
                                  h_mat.p,
                                  one,
                                  one,
                                  h_mat.desc,
                                  T(0.0),
                                  h_work.data(),
                                  one,
                                  one,
                                  h_mat.desc);
        ScalapackConnector::tranc(nFull,
                                  nFull,
                                  T(1.0),
                                  s_mat.p,
                                  one,
                                  one,
                                  s_mat.desc,
                                  T(0.0),
                                  s_work.data(),
                                  one,
                                  one,
                                  s_mat.desc);
        const auto local_to_global = [](const int local_index,
                                        const int block_size,
                                        const int process_coordinate,
                                        const int source_coordinate,
                                        const int process_count) {
            if (source_coordinate < 0)
            {
                return local_index;
            }
            const int process_offset
                = (process_coordinate - source_coordinate + process_count) % process_count;
            return ((local_index / block_size) * process_count + process_offset) * block_size
                   + local_index % block_size;
        };
        for (int local_col = 0; local_col < nacols; ++local_col)
        {
            const int global_col
                = local_to_global(local_col, h_mat.desc[5], mypcol, h_mat.desc[7], npcols);
            for (int local_row = 0; local_row < narows; ++local_row)
            {
                const int global_row
                    = local_to_global(local_row, h_mat.desc[4], myprow, h_mat.desc[6], nprows);
                if (global_row <= global_col)
                {
                    const int local_index = local_row + local_col * narows;
                    h_work[local_index] = h_mat.p[local_index];
                    s_work[local_index] = s_mat.p[local_index];
                }
            }
        }
        h_elpa = h_work.data();
        s_elpa = s_work.data();
        decomposed_state = 0;
    }

    if (elpa_init(20210430) != ELPA_OK)
    {
        fprintf(stderr, "Error: ELPA API version not supported");
    }

    // elpa_init(20210430);
    int success = 0;
    elpa_t handle = elpa_allocate(&success);
#ifdef _OPENMP
    int num_threads = omp_get_max_threads();
#else
    int num_threads = 1;
#endif
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
    if (this->use_gpu)
    {
        elpa_set(handle, "nvidia-gpu", 1, &success);
        elpa_set(handle, "real_kernel", ELPA_2STAGE_REAL_NVIDIA_GPU, &success);
        elpa_setup_gpu(handle);
    }
#endif

    elpa_generalized_eigenvectors(handle,
                                  h_elpa,
                                  s_elpa,
                                  eigen.data(),
                                  eigenvectors.data(),
                                  decomposed_state,
                                  &success);
    elpa_deallocate(handle, &success);
    elpa_uninit(&success);

    ModuleBase::timer::end("DiagoElpaNative", "elpa_solve");
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
    BlasConnector::copy(this->nbands, eigen.data(), inc, eigenvalue_in, inc);
    const int size = psi.get_nbands() * psi.get_nbasis();
    BlasConnector::copy(size, eigenvectors.data(), inc, psi.get_pointer(), inc);
}
#endif

template <typename T>
void DiagoElpaNative<T>::diag(ModuleBase::MatrixBlock<T>& h_mat,
                              ModuleBase::MatrixBlock<T>& s_mat,
                              psi::Psi<T>& psi,
                              Real* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoElpaNative", "diag");
#ifdef __MPI
    MPI_Comm COMM_DIAG = setmpicomm(); // set mpi_comm needed
    diag_pool(h_mat, s_mat, psi, eigenvalue_in, COMM_DIAG);
#else
    ModuleBase::WARNING_QUIT("DiagoElpaNative", "DiagoElpaNative only can be used with macro __MPI");
#endif
}

template class DiagoElpaNative<double>;
template class DiagoElpaNative<std::complex<double>>;

} // namespace hsolver
