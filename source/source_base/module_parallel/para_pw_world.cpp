#include "para_pw_world.h"

#include <algorithm>
#include <cassert>

namespace Parallel
{

ParaPwWorld::ParaPwWorld(int npw)
    : ParaWorld("pw"), npw_(npw), npwtot_(npw), npw_per_(1, npw)
{
}

#ifdef __MPI
ParaPwWorld::ParaPwWorld(const MPI_Comm& comm, const std::vector<int>& npw_per)
    : ParaWorld("pw", comm), npw_per_(npw_per)
{
    npw_ = npw_per_[rank()];
    npwtot_ = 0;
    for (int n : npw_per_)
    {
        npwtot_ += n;
    }
}
#endif

int ParaPwWorld::npw_per(int p) const
{
    assert(p >= 0 && p < static_cast<int>(npw_per_.size()));
    return npw_per_[p];
}

} // namespace Parallel
