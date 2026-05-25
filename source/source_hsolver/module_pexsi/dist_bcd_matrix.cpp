#ifdef __PEXSI
#include "dist_bcd_matrix.h"
#include <cctype>
#include <mpi.h>
extern "C"
{
    void Cblacs_gridinfo(int icontxt, int* nprow, int* npcol, int* myprow, int* mypcol);
    int Cblacs_pnum(int blacs_ctxt, int prow, int pcol);
};

namespace pexsi
{
DistBCDMatrix::DistBCDMatrix(MPI_Comm comm,
                             MPI_Group group,
                             int blacs_ctxt,
                             int size,
                             int nblk,
                             int nrow,
                             int ncol,
                             char layout)
{
    this->comm = comm;
    this->group = group;
    this->blacs_ctxt = blacs_ctxt;
    this->size = size;
    this->nblk = nblk;
    this->nrow = nrow;
    this->ncol = ncol;
    if (layout == 'r' || layout == 'c')
    {
        this->layout = layout;
    }
    else
    {
        throw("The layout must be 'r' or 'c'");
    }

    if (comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank(comm, &this->myproc);
        Cblacs_gridinfo(blacs_ctxt, &this->nprows, &this->npcols, &this->myprow, &this->mypcol);
    }
    else
    {
        this->myproc = -1;
        this->myprow = -1;
        this->mypcol = -1;
        this->nprocs = 0;
        this->nprows = 0;
        this->npcols = 0;
        this->prowpcol2pnum = nullptr;
        return;
    }

    // Synchronize matrix parameters inside the BCD communicator.  Older code used
    // MPI_COMM_WORLD here, which prevents multiple independent k-point pools from
    // constructing BCD descriptors concurrently.
    int myid_in_comm = 0;
    MPI_Comm_rank(comm, &myid_in_comm);
    if (myid_in_comm == 0)
    {
        MPI_Comm_size(comm, &this->nprocs);
        int PARA_BCAST[4] = {this->nblk, this->nprocs, this->nprows, this->npcols};
        MPI_Bcast(&PARA_BCAST[0], 4, MPI_INT, 0, comm);
    }
    else
    {
        int PARA_BCAST[4];
        MPI_Bcast(&PARA_BCAST[0], 4, MPI_INT, 0, comm);
        this->nblk = PARA_BCAST[0];
        this->nprocs = PARA_BCAST[1];
        this->nprows = PARA_BCAST[2];
        this->npcols = PARA_BCAST[3];
    }
    this->prowpcol2pnum = new int[this->nprocs];
    if (myid_in_comm == 0)
    {
        for (int i = 0; i < this->nprows; ++i)
        {
            for (int j = 0; j < this->npcols; ++j)
            {
                this->prowpcol2pnum[i * this->npcols + j] = Cblacs_pnum(this->blacs_ctxt, i, j);
            }
        }
    }
    MPI_Bcast(this->prowpcol2pnum, this->nprocs, MPI_INT, 0, comm);
}

DistBCDMatrix::~DistBCDMatrix()
{
    delete[] prowpcol2pnum;
}

int DistBCDMatrix::globalRow(const int localRow)
{
    return (localRow / nblk * nprows + myprow) * nblk + localRow % nblk;
}

int DistBCDMatrix::globalCol(const int localCol)
{
    return (localCol / nblk * npcols + mypcol) * nblk + localCol % nblk;
}

int DistBCDMatrix::localRow(const int globalRow, int& myprow)
{
    myprow = int((globalRow % (nblk * nprows)) / nblk);
    return int(globalRow / (nblk * nprows)) * nblk + globalRow % nblk;
}

int DistBCDMatrix::localCol(const int globalCol, int& mypcol)
{
    mypcol = int((globalCol % (nblk * npcols)) / nblk);
    return int(globalCol / (nblk * npcols)) * nblk + globalCol % nblk;
}

int DistBCDMatrix::pnum(const int prow, const int pcol)
{
    return this->prowpcol2pnum[prow * this->npcols + pcol];
}
} // namespace pexsi
#endif
