#ifndef GATHER_MAT_H
#define GATHER_MAT_H

#include "source_base/module_external/scalapack_connector.h" // Cpxgemr2d
#include "source_hamilt/matrixblock.h"

namespace module_rt
{
//------------------------ MPI gathering and distributing functions ------------------------//
// This struct is used for collecting matrices from all processes to root process
template <typename T>
struct Matrix_g
{
    std::shared_ptr<T> p;
    size_t row;
    size_t col;
    std::shared_ptr<int> desc;
};

#ifdef __MPI
// Collect matrices from all processes to root process
template <typename T>
void gatherMatrix(const int myid, const int root_proc, const hamilt::MatrixBlock<T>& mat_l, Matrix_g<T>& mat_g)
{
    const int* desca = mat_l.desc; // Obtain the descriptor of the local matrix
    int ctxt = desca[1];           // BLACS context
    int nrows = desca[2];          // Global matrix row number
    int ncols = desca[3];          // Global matrix column number

    if (myid == root_proc)
    {
        mat_g.p.reset(new T[nrows * ncols]); // No need to delete[] since it is a shared_ptr
    }
    else
    {
        mat_g.p.reset(new T[nrows * ncols]); // Placeholder for non-root processes
    }

    // Set the descriptor of the global matrix
    mat_g.desc.reset(new int[9]{1, ctxt, nrows, ncols, nrows, ncols, 0, 0, nrows});
    mat_g.row = nrows;
    mat_g.col = ncols;

    // Call the Cpxgemr2d function in ScaLAPACK to collect the matrix data
    Cpxgemr2d(nrows, ncols, mat_l.p, 1, 1, const_cast<int*>(desca), mat_g.p.get(), 1, 1, mat_g.desc.get(), ctxt);
}

template <typename T>
void gatherPsi(const int myid,
               const int root_proc,
               T* psi_l,
               const Parallel_Orbitals& para_orb,
               module_rt::Matrix_g<T>& psi_g)
{
    const int* desc_psi = para_orb.desc_wfc; // Obtain the descriptor from Parallel_Orbitals
    int ctxt = desc_psi[1];                  // BLACS context
    int nrows = desc_psi[2];                 // Global matrix row number
    int ncols = desc_psi[3];                 // Global matrix column number

    if (myid == root_proc)
    {
        psi_g.p.reset(new T[nrows * ncols]); // No need to delete[] since it is a shared_ptr
    }
    else
    {
        psi_g.p.reset(new T[nrows * ncols]); // Placeholder for non-root processes
    }

    // Set the descriptor of the global psi
    psi_g.desc.reset(new int[9]{1, ctxt, nrows, ncols, nrows, ncols, 0, 0, nrows});
    psi_g.row = nrows;
    psi_g.col = ncols;

    // Call the Cpxgemr2d function in ScaLAPACK to collect the matrix data
    Cpxgemr2d(nrows, ncols, psi_l, 1, 1, const_cast<int*>(desc_psi), psi_g.p.get(), 1, 1, psi_g.desc.get(), ctxt);
}

template <typename T>
void distributePsi(const Parallel_Orbitals& para_orb, T* psi_l, const module_rt::Matrix_g<T>& psi_g)
{
    const int* desc_psi = para_orb.desc_wfc; // Obtain the descriptor from Parallel_Orbitals
    int ctxt = desc_psi[1];                  // BLACS context
    int nrows = desc_psi[2];                 // Global matrix row number
    int ncols = desc_psi[3];                 // Global matrix column number

    // Call the Cpxgemr2d function in ScaLAPACK to distribute the matrix data
    Cpxgemr2d(nrows, ncols, psi_g.p.get(), 1, 1, psi_g.desc.get(), psi_l, 1, 1, const_cast<int*>(desc_psi), ctxt);
}
//------------------------ MPI gathering and distributing functions ------------------------//

#endif // __MPI
} // namespace module_rt
#endif // GATHER_MAT_H