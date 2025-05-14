#include "module_hsolver/diag_hs_para.h"

#include "module_base/scalapack_connector.h"
#include "module_base/parallel_2d.h"
#include "module_hsolver/diago_pxxxgvx.h"

#ifdef __ELPA
#include "module_hsolver/genelpa/elpa_solver.h"
#endif

#include <iostream>

namespace hsolver
{

#ifdef __ELPA
void elpa_diag(MPI_Comm comm,
               const int nband,
               std::complex<double>* h_local,
               std::complex<double>* s_local,
               double* ekb,
               std::complex<double>* wfc_2d,
               Parallel_2D& para2d_local)
{
    int DecomposedState = 0;
    ELPA_Solver es(false, comm, nband, para2d_local.get_row_size(), para2d_local.get_col_size(), para2d_local.desc);
    es.generalized_eigenvector(h_local, s_local, DecomposedState, ekb, wfc_2d);
    es.exit();
}

void elpa_diag(MPI_Comm comm,
               const int nband,
               double* h_local,
               double* s_local,
               double* ekb,
               double* wfc_2d,
               Parallel_2D& para2d_local)
{
    int DecomposedState = 0;
    ELPA_Solver es(true, comm, nband, para2d_local.get_row_size(), para2d_local.get_col_size(), para2d_local.desc);
    es.generalized_eigenvector(h_local, s_local, DecomposedState, ekb, wfc_2d);
    es.exit();
}

void elpa_diag(MPI_Comm comm,
               const int nband,
               std::complex<float>* h_local,
               std::complex<float>* s_local,
               float* ekb,
               std::complex<float>* wfc_2d,
               Parallel_2D& para2d_local)
{
    std::cout << "Error: ELPA do not support single precision. " << std::endl;
    exit(1);
}

void elpa_diag(MPI_Comm comm,
               const int nband,
               float* h_local,
               float* s_local,
               float* ekb,
               float* wfc_2d,
               Parallel_2D& para2d_local)
{
    std::cout << "Error: ELPA do not support single precision. " << std::endl;
    exit(1);
}

#endif

#ifdef __MPI

template <typename T>
void diago_hs_para(T* h,
                   T* s,
                   const int lda,
                   const int nband,
                   typename GetTypeReal<T>::type* const ekb,
                   T* const wfc,
                   const MPI_Comm& comm,
                   const int diag_subspace,
                   const int block_size)
{
    int myrank = 0;
    MPI_Comm_rank(comm, &myrank);
    Parallel_2D para2d_global;
    Parallel_2D para2d_local;
    para2d_global.init(lda, lda, lda, comm);

    int max_nb = block_size;
    if (block_size == 0)
    {
        if (nband > 500)
        {
            max_nb = 32;
        }
        else
        {
            max_nb = 16;
        }
    }
    else if (block_size < 0)
    {
        std::cout << "Error: block_size in diago_subspace should be a positive integer. " << std::endl;
        exit(1);
    }

    // for genelpa, if the block size is too large that some cores have no data, then it will cause error.
    if (diag_subspace == 1)
    {
        if (max_nb * (std::max(para2d_global.dim0, para2d_global.dim1) - 1) >= lda)
        {
            max_nb = lda / std::max(para2d_global.dim0, para2d_global.dim1);
        }
    }

    para2d_local.init(lda, lda, max_nb, comm);
    std::vector<T> h_local(para2d_local.get_col_size() * para2d_local.get_row_size());
    std::vector<T> s_local(para2d_local.get_col_size() * para2d_local.get_row_size());
    std::vector<T> wfc_2d(para2d_local.get_col_size() * para2d_local.get_row_size());

    // distribute h and s to 2D
    Cpxgemr2d(lda, lda, h, 1, 1, para2d_global.desc, h_local.data(), 1, 1, para2d_local.desc, para2d_local.blacs_ctxt);
    Cpxgemr2d(lda, lda, s, 1, 1, para2d_global.desc, s_local.data(), 1, 1, para2d_local.desc, para2d_local.blacs_ctxt);

    if (diag_subspace == 1)
    {
#ifdef __ELPA
        elpa_diag(comm, nband, h_local.data(), s_local.data(), ekb, wfc_2d.data(), para2d_local);
#else
        std::cout << "ERROR: try to use ELPA to solve the generalized eigenvalue problem, but ELPA is not compiled. "
                  << std::endl;
        exit(1);
#endif
    }
    else if (diag_subspace == 2)
    {
        hsolver::pxxxgvx_diag(para2d_local.desc,
                              para2d_local.get_row_size(),
                              para2d_local.get_col_size(),
                              nband,
                              h_local.data(),
                              s_local.data(),
                              ekb,
                              wfc_2d.data());
    }
    else
    {
        std::cout << "Error: parallel diagonalization method is not supported. " << "diag_subspace = " << diag_subspace
                  << std::endl;
        exit(1);
    }

    // gather wfc
    Cpxgemr2d(lda, lda, wfc_2d.data(), 1, 1, para2d_local.desc, wfc, 1, 1, para2d_global.desc, para2d_local.blacs_ctxt);

    // free the context
    Cblacs_gridexit(para2d_local.blacs_ctxt);
    Cblacs_gridexit(para2d_global.blacs_ctxt);
}

// template instantiation
template void diago_hs_para<double>(double* h,
                                    double* s,
                                    const int lda,
                                    const int nband,
                                    typename GetTypeReal<double>::type* const ekb,
                                    double* const wfc,
                                    const MPI_Comm& comm,
                                    const int diag_subspace,
                                    const int block_size);

template void diago_hs_para<std::complex<double>>(std::complex<double>* h,
                                                  std::complex<double>* s,
                                                  const int lda,
                                                  const int nband,
                                                  typename GetTypeReal<std::complex<double>>::type* const ekb,
                                                  std::complex<double>* const wfc,
                                                  const MPI_Comm& comm,
                                                  const int diag_subspace,
                                                  const int block_size);

template void diago_hs_para<float>(float* h,
                                   float* s,
                                   const int lda,
                                   const int nband,
                                   typename GetTypeReal<float>::type* const ekb,
                                   float* const wfc,
                                   const MPI_Comm& comm,
                                   const int diag_subspace,
                                   const int block_size);

template void diago_hs_para<std::complex<float>>(std::complex<float>* h,
                                                 std::complex<float>* s,
                                                 const int lda,
                                                 const int nband,
                                                 typename GetTypeReal<std::complex<float>>::type* const ekb,
                                                 std::complex<float>* const wfc,
                                                 const MPI_Comm& comm,
                                                 const int diag_subspace,
                                                 const int block_size);

#endif

} // namespace hsolver