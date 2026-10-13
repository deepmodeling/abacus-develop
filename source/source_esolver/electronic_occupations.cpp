#include "electronic_occupations.h"
#include "source_base/matrix.h"
#include "source_base/parallel_reduce.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace ModuleESolver
{
double checked_occupation_electron_count(const ModuleBase::matrix& weights,
                                         double requested_nelect,
                                         int npool,
                                         int ranks_per_pool)
{
    int invalid = !std::isfinite(requested_nelect) || requested_nelect <= 0.0
                  || npool <= 0 || ranks_per_pool <= 0;
    double occupied = 0.0;
    double compensation = 0.0;
    for (int ik = 0; ik < weights.nr; ++ik)
    {
        for (int ib = 0; ib < weights.nc; ++ib)
        {
            const double weight = weights(ik, ib);
            if (!std::isfinite(weight))
            {
                invalid = 1;
                continue;
            }
            // Compensated sum of the local k/spin and band shard. Generalized
            // smearing weights may be negative; only their total is constrained.
            const double corrected = weight - compensation;
            const double next = occupied + corrected;
            compensation = (next - occupied) - corrected;
            occupied = next;
        }
    }
    if (!std::isfinite(occupied)) { invalid = 1; }
    // Synchronize failures before a rank can enter a different collective path.
    Parallel_Reduce::reduce_all(invalid);
    if (invalid)
    {
        throw std::runtime_error("constant-potential electron count requires finite occupations "
                                 "and a valid target/occupation partition");
    }

    // Same distribution as the native band-energy sum: remove the replicated
    // copies inside each pool before combining KPAR and BNDPAR contributions.
    Parallel_Reduce::reduce_double_allpool(npool, ranks_per_pool, occupied);
    const double tolerance = 1.0e-8;
    int mismatch = !std::isfinite(occupied) || occupied <= 0.0
                   || std::abs(occupied - requested_nelect) > tolerance;
    Parallel_Reduce::reduce_all(mismatch);
    if (mismatch)
    {
        std::ostringstream message;
        message << std::setprecision(17)
                << "constant-potential occupation electron count mismatch: requested=" << requested_nelect
                << ", occupied=" << occupied << ", absolute tolerance=" << tolerance;
        throw std::runtime_error(message.str());
    }
    return occupied;
}
}
