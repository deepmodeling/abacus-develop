#include "para_world.h"

namespace Parallel
{

ParaWorld::ParaWorld(const std::string& tag) : tag_(tag), rank_(0), size_(1)
{
#ifdef __MPI
    if (!tag.empty())
    {
        comm_ = MPI_COMM_SELF;
    }
    else
    {
        comm_ = MPI_COMM_NULL;
    }
#endif
}

#ifdef __MPI
ParaWorld::ParaWorld(const std::string& tag, const MPI_Comm& comm) : tag_(tag), comm_(comm)
{
    if (comm == MPI_COMM_NULL)
    {
        rank_ = -1;
        size_ = 0;
        return;
    }
    MPI_Comm_rank(comm, &rank_);
    MPI_Comm_size(comm, &size_);
}
#endif

bool ParaWorld::valid() const
{
#ifdef __MPI
    return comm_ != MPI_COMM_NULL;
#else
    return !tag_.empty();
#endif
}

} // namespace Parallel
