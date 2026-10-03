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

int CommunicationDomain::size() const
{
    int result = 1;
#ifdef __MPI
    if (communicator_ != MPI_COMM_NULL)
    {
        MPI_Comm_size(communicator_, &result);
    }
#endif
    return result;
}

int CommunicationDomain::max(int value) const
{
#ifdef __MPI
    if (communicator_ != MPI_COMM_NULL)
    {
        MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, communicator_);
    }
#endif
    return value;
}

double CommunicationDomain::max(double value) const
{
#ifdef __MPI
    if (communicator_ != MPI_COMM_NULL)
    {
        MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_MAX, communicator_);
    }
#endif
    return value;
}

CommunicationDomain world_comm_domain()
{
    CommunicationDomain comm_domain;
#ifdef __MPI
    comm_domain.initialize(MPI_COMM_WORLD);
#endif
    return comm_domain;
}
} // namespace ModuleBase
