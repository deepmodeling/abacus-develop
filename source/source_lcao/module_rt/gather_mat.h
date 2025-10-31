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
//------------------------ MPI gathering and distributing functions ------------------------//

} // namespace module_rt
#endif // GATHER_MAT_H