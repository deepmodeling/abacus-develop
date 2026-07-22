#ifndef ABACUS_SOURCE_CONTEXT_SIMULATION_CONTEXT_BUILDER_H
#define ABACUS_SOURCE_CONTEXT_SIMULATION_CONTEXT_BUILDER_H

#ifndef ABACUS_CAN_BIND_SIMULATION_CONTEXT
#error "simulation_context_builder.h is restricted to the assembly layer"
#endif

#include "source_context/simulation_context.h"

struct Input_para;
struct System_para;

namespace ModuleContext
{

class SimulationContextBuilder
{
  public:
    SimulationContextBuilder(const Input_para& input, const System_para& system);

    void capture_runtime(const System_para& system);
    SimulationContext finalize(const System_para& system);
    bool runtime_captured() const { return runtime_captured_; }
    bool finalized() const { return finalized_; }

  private:
    SimulationContext context_;
    bool runtime_captured_ = false;
    bool finalized_ = false;
};

} // namespace ModuleContext

#endif
