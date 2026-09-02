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
    CommunicationDomain() = default;
#ifdef __MPI
    void initialize(MPI_Comm communicator);
    MPI_Comm communicator() const;
#endif
    int rank() const;

private:
#ifdef __MPI
    MPI_Comm communicator_ = MPI_COMM_NULL;
#endif
    int rank_ = 0;
};

CommunicationDomain world_communication_domain();
} // namespace ModuleBase

#endif
