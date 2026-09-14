#ifndef MD_STATISTICS_H
#define MD_STATISTICS_H

#include "source_base/matrix.h"

namespace MD_func
{

struct MDKineticState
{
    double kinetic = 0.0;
    double temperature = 0.0;
};

struct MDStressState
{
    ModuleBase::matrix stress;
    ModuleBase::matrix temperature_tensor;
};

} // namespace MD_func

#endif // MD_STATISTICS_H
