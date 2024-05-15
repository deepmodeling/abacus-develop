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

    const MPI_Comm get_comm() const
    {
        return comm;
    };
    const MPI_Group get_group() const
    {
        return group;
    };
    const MPI_Group get_group_data() const
    {
        return group_data;
    };
    const int get_size() const
    {
        return size;
    };
    const int get_nnz() const
    {
        return nnz;
    };
    const int get_nnzlocal() const
    {
        return nnzLocal;
    };
    const int get_numcol_local() const
    {
        return numColLocal;
    };
    int* get_colptr_local() const
    {
        return colptrLocal;
    };
    int* get_rowind_local() const
    {
        return rowindLocal;
    };

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

    // friend class DistMatrixTransformer;
};
} // namespace pexsi
#endif // DISTCCSMATRIX_H
