#ifndef DISTBCDMATRIX_H
#define DISTBCDMATRIX_H

#include <mpi.h>

#include "module_hsolver/module_pexsi/dist_matrix_transformer.h"
// a Block Cyclic Data Distribution matrix
// http://www.netlib.org/utk/papers/factor/node3.html
// local matrix elements is stored in column major
// used for pexsi
namespace pexsi
{
class DistBCDMatrix
{

  public:
    // DistBCDMatrix(MPI_Comm comm, MPI_Group group, int nprow, int npcol, int size, int nblk, int nrow, int ncol);
    // DistBCDMatrix(MPI_Comm comm, MPI_Group group, int nprow, int npcol, int size, int nblk, int nrow, int ncol, char
    // LAYOUT);

    // DistBCDMatrix(MPI_Comm comm, MPI_Group group, int blacs_ctxt, int size, int nblk, int nrow, int ncol);
    DistBCDMatrix(MPI_Comm comm, MPI_Group group, int blacs_ctxt, int size, int nblk, int nrow, int ncol, char LAYOUT);
    ~DistBCDMatrix();

    int globalRow(const int localRow);
    int globalCol(const int localCol);
    int localRow(const int globalRow, int& myprow);
    int localCol(const int globalCol, int& mypcol);
    int pnum(const int prow, const int pcol);
    //~DistBCDMatrix();

    const MPI_Comm get_comm() const
    {
        return comm;
    };
    const MPI_Group get_group() const
    {
        return group;
    };
    const int get_nrow() const
    {
        return nrow;
    };
    const int get_ncol() const
    {
        return ncol;
    };
    const char get_LAYOUT() const
    {
        return LAYOUT;
    };

  private:
    // MPI communicator
    MPI_Comm comm;
    MPI_Group group;

    // blacs context
    int blacs_ctxt;

    // row and column of process grid
    int nprows;
    int npcols;

    // total number of processes
    int nprocs;

    // Matrix size
    int size;

    // block size
    int nblk;

    // row and c0lumn of Local matrix part
    int nrow;
    int ncol;

    // protected:

    // private:

    // current process row and column
    int myprow;
    int mypcol;

    // current process id
    int myproc;

    int* prowpcol2pnum;
    // the local data layout
    // 'R' or 'r' for row-major, which is used in C/C++
    // 'C' or 'c' for column-major, which is used in Fortran
    char LAYOUT;
};
} // namespace pexsi
#endif // DISTBCDMATRIX_H