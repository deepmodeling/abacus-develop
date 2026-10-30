#include "source_base/parallel_cell.h"

#include <stdexcept>
#include <vector>

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

int CommunicationDomain::sum(int value) const
{
#ifdef __MPI
    if (communicator_ != MPI_COMM_NULL)
    {
        MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_SUM, communicator_);
    }
#endif
    return value;
}

void CommunicationDomain::sum(int* values, int count) const
{
#ifdef __MPI
    if (communicator_ != MPI_COMM_NULL)
    {
        MPI_Allreduce(MPI_IN_PLACE, values, count, MPI_INT, MPI_SUM, communicator_);
    }
#endif
}

void CommunicationDomain::sum(double* values, int count) const
{
#ifdef __MPI
    if (communicator_ != MPI_COMM_NULL)
    {
        MPI_Allreduce(MPI_IN_PLACE, values, count, MPI_DOUBLE, MPI_SUM, communicator_);
    }
#endif
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

std::string CommunicationDomain::gather_text(const std::string& value, int maximum_count) const
{
    if (maximum_count <= 0)
    {
        throw std::invalid_argument("CommunicationDomain text batch limit must be positive.");
    }
    const std::size_t count_limit = static_cast<std::size_t>(maximum_count);
    const int local_overflow = value.size() > count_limit;
    if (max(local_overflow) != 0)
    {
        throw std::length_error("CommunicationDomain text exceeds the batch limit.");
    }
#ifdef __MPI
    if (communicator_ != MPI_COMM_NULL)
    {
        const int local_size = static_cast<int>(value.size());
        const int process_count = size();
        std::vector<int> lengths(process_count, 0);
        std::vector<int> offsets(process_count, 0);
        MPI_Gather(&local_size, 1, MPI_INT, lengths.data(), 1, MPI_INT, 0, communicator_);
        std::size_t total_size = 0;
        int total_overflow = 0;
        if (rank_ == 0)
        {
            for (int rank = 0; rank < process_count; ++rank)
            {
                if (total_size + lengths[rank] > count_limit)
                {
                    total_overflow = 1;
                    break;
                }
                offsets[rank] = static_cast<int>(total_size);
                total_size += lengths[rank];
            }
        }
        MPI_Bcast(&total_overflow, 1, MPI_INT, 0, communicator_);
        if (total_overflow != 0)
        {
            throw std::length_error("CommunicationDomain gathered text exceeds the batch limit.");
        }
        std::string gathered;
        if (rank_ == 0)
        {
            gathered.resize(total_size);
        }
        char* buffer = nullptr;
        if (!gathered.empty())
        {
            buffer = &gathered[0];
        }
        MPI_Gatherv(value.data(), local_size, MPI_CHAR, buffer, lengths.data(), offsets.data(), MPI_CHAR, 0, communicator_);
        return gathered;
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
