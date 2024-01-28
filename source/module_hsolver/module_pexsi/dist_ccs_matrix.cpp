#ifdef __PEXSI
#include "dist_ccs_matrix.h"

#include <mpi.h>

namespace pexsi
{
DistCCSMatrix::DistCCSMatrix(void)
{
    this->comm = MPI_COMM_WORLD;
    this->size = 0;
    this->nnz = 0;
    this->nnzLocal = 0;
    this->numColLocal = 0;
    this->colptrLocal = NULL;
    this->rowindLocal = NULL;
}

DistCCSMatrix::DistCCSMatrix(MPI_Comm comm_in)
{
    this->comm = comm_in;
    this->size = 0;
    this->nnz = 0;
    this->nnzLocal = 0;
    this->numColLocal = 0;
    this->colptrLocal = NULL;
    this->rowindLocal = NULL;
}

DistCCSMatrix::DistCCSMatrix(int size_in, int nnzLocal_in)
{
    this->comm = MPI_COMM_WORLD;
    this->size = size_in;
    this->nnzLocal = nnzLocal_in;
    MPI_Request req;
    MPI_Iallreduce(&nnzLocal, &this->nnz, 1, MPI_INT, MPI_SUM, this->comm, &req);
    this->numColLocal = 0;
    this->colptrLocal = new int[size];
    this->rowindLocal = new int[nnzLocal];

    MPI_Status req_status;
    MPI_Wait(&req, &req_status);
}

DistCCSMatrix::DistCCSMatrix(MPI_Comm comm_in, int nproc_data_in, int size_in)
{
    this->comm = comm_in;
    this->nproc_data = nproc_data_in;
    int nproc_data_range[3] = {0, this->nproc_data - 1, 1};
    // create processes group with data: this->group_data and associated communicator
    MPI_Comm_group(this->comm, &this->group);
    MPI_Group_range_incl(this->group, 1, &nproc_data_range, &this->group_data);
    this->comm_data = MPI_COMM_NULL;
    MPI_Comm_create(this->comm, this->group_data, &this->comm_data);
    this->size = size_in;
    this->nnz = 0;
    this->nnzLocal = 0;
    int myproc;
    if (comm != MPI_COMM_NULL)
    {
        MPI_Comm_size(comm, &nprocs);
        MPI_Comm_rank(comm, &myproc);
        if (myproc < nproc_data - 1)
        {
            this->numColLocal = size / nproc_data;
            this->firstCol = size / nproc_data * myproc;
            this->colptrLocal = new int[this->numColLocal + 1];
            this->rowindLocal = NULL;
        }
        else if (myproc == nproc_data - 1)
        {
            this->numColLocal = size - myproc * (size / nproc_data);
            this->firstCol = size / nproc_data * myproc;
            this->colptrLocal = new int[this->numColLocal + 1];
            this->rowindLocal = NULL;
        }
        else
        {
            this->numColLocal = 0;
            this->firstCol = size - 1;
            this->colptrLocal = new int[this->numColLocal + 1];
            this->rowindLocal = NULL;
        }
    }
}

int DistCCSMatrix::globalCol(int localCol)
{
    return this->firstCol + localCol;
}

// NOTE: the process id is 0-based
int DistCCSMatrix::localCol(int globalCol, int& mypcol)
{
    mypcol = int(globalCol / int(this->size / this->nproc_data));
    if (mypcol >= this->nproc_data)
        mypcol = this->nproc_data - 1;

    return mypcol > 0 ? globalCol - (this->size / this->nproc_data) * mypcol : globalCol;
}

void DistCCSMatrix::setnnz(int nnzLocal_in)
{
    if (this->comm_data != MPI_COMM_NULL)
    {
        MPI_Allreduce(&nnzLocal_in, &this->nnz, 1, MPI_INT, MPI_SUM, this->comm_data);
        this->nnzLocal = nnzLocal_in;
        this->rowindLocal = new int[nnzLocal];
        this->colptrLocal[this->numColLocal] = nnzLocal_in + 1;
    }
}

DistCCSMatrix::~DistCCSMatrix()
{
    delete[] colptrLocal;
    delete[] rowindLocal;
}
} // namespace pexsi
#endif