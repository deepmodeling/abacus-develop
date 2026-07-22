#ifndef ABACUS_SOURCE_CONTEXT_ORCHESTRATION_CONTEXT_H
#define ABACUS_SOURCE_CONTEXT_ORCHESTRATION_CONTEXT_H

#ifndef ABACUS_CAN_READ_SIMULATION_CONTEXT
#error "orchestration_context.h is restricted to L0/L1 orchestration targets"
#endif

#include "source_context/simulation_context.h"

namespace ModuleContext
{

const SimulationContext& current_simulation_context();

} // namespace ModuleContext

#endif
