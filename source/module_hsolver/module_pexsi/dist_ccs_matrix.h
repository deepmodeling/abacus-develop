#ifndef DISTCCSMATRIX_H
#define DISTCCSMATRIX_H

#include <mpi.h>
// Distributed Compressed Column Storage Matrix format
// used for PEXSI
namespace pexsi
{
class DistCCSMatrix
{

  public:
    DistCCSMatrix();
    DistCCSMatrix(MPI_Comm comm);
    DistCCSMatrix(int size, int nnzLocal);
    DistCCSMatrix(MPI_Comm comm, int size, int nnzLocal);
    DistCCSMatrix(MPI_Comm comm, int size, int nnzLocal, double* valLocal, int* index);

    int globalCol(int localCol);
    int localCol(int globalCol, int& mypcol);
    void setnnz(int nnzLocal);
    ~DistCCSMatrix();

  private:
    // MPI communicator
    MPI_Comm comm;
    MPI_Group group;

    // total number of processes and the processes with data in
    int nprocs;
    int nproc_data;
    MPI_Group group_data;
    MPI_Comm comm_data;

    // Matrix size
    int size;

    // Number of non-zero values in the matrix
    int nnz;

    // Number of non-zero values in the matrix of the local process
    int nnzLocal;

    // number of columns in current process
    int numColLocal;

    // the first column index in current process
    int firstCol;

    // Array stores the indices to the nonzero row indices in rowptrLocal and nzvalLocal
    int* colptrLocal;
    int* rowindLocal;
};
} // namespace pexsi
#endif // DISTCCSMATRIX_H
