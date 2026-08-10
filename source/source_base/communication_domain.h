#ifndef COMMUNICATION_DOMAIN_H
#define COMMUNICATION_DOMAIN_H

#ifdef __MPI
#include <mpi.h>
#endif

namespace ModuleBase
{
/**
 * @brief A non-owning MPI communication domain supplied to a subsystem.
 *
 * The creator remains responsible for the lifetime of the MPI communicator.
 */
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
