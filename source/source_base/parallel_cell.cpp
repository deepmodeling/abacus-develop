#include "source_base/parallel_cell.h"

namespace ModuleBase
{
#ifdef __MPI
void CommunicationDomain::initialize(MPI_Comm communicator)
{
    communicator_ = communicator;
    rank_ = 0;
    if (communicator_ != MPI_COMM_NULL)
    {
        MPI_Comm_rank(communicator_, &rank_);
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

CommunicationDomain world_communication_domain()
{
    CommunicationDomain communication_domain;
#ifdef __MPI
    communication_domain.initialize(MPI_COMM_WORLD);
#endif
    return communication_domain;
}
} // namespace ModuleBase
