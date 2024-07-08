#ifndef DIAG_COMM
#define DIAG_COMM

#ifdef __MPI
#include "mpi.h"
#endif

namespace hsolver
{

struct diag_comm_info
{

    const int rank;
    const int nproc;

#ifndef __MPI
    diag_comm_info(const int rank_in, const int nproc_in) : rank(rank_in), nproc(nproc_in)
    {
    }
#else
    const MPI_Comm comm;
    diag_comm_info(const MPI_Comm& comm_in, const int rank_in, const int nproc_in)
        : comm(comm_in), rank(rank_in), nproc(nproc_in)
    {
    }
#endif
};

} // namespace hsolver

#endif