#ifndef ELECSTATE_TOOLS_H
#define ELECSTATE_TOOLS_H
#include "elecstate.h"

namespace elecstate
{
    void calEBand(const ModuleBase::matrix& ekb,const ModuleBase::matrix& wg,fenergy& f_en);
}

#endif