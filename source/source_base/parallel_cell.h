#ifndef PARALLEL_CELL_H
#define PARALLEL_CELL_H

#ifdef __MPI
#include <mpi.h>
#endif

namespace ModuleBase
{
class CommunicationDomain
{
public:
    CommunicationDomain();
#ifdef __MPI
    explicit CommunicationDomain(MPI_Comm communicator);
    MPI_Comm communicator() const;
#endif
    int rank() const;
    int size() const;

private:
#ifdef __MPI
    MPI_Comm communicator_ = MPI_COMM_NULL;
#endif
    int rank_ = 0;
    int size_ = 1;
};

CommunicationDomain world_communication_domain();
} // namespace ModuleBase

#endif
