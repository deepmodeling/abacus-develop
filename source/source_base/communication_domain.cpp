#include "source_base/communication_domain.h"

namespace ModuleBase
{
CommunicationDomain::CommunicationDomain()
{
}

#ifdef __MPI
CommunicationDomain::CommunicationDomain(MPI_Comm communicator) : communicator_(communicator)
{
    if (communicator_ != MPI_COMM_NULL)
    {
        MPI_Comm_rank(communicator_, &rank_);
        MPI_Comm_size(communicator_, &size_);
    }
}

MPI_Comm CommunicationDomain::communicator() const
{
    return communicator_;
}
#endif

int CommunicationDomain::rank() const
{
    return rank_;
}

int CommunicationDomain::size() const
{
    return size_;
}

CommunicationDomain world_communication_domain()
{
#ifdef __MPI
    return CommunicationDomain(MPI_COMM_WORLD);
#else
    return CommunicationDomain();
#endif
}
} // namespace ModuleBase
