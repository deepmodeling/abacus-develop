#ifndef ABACUS_SOURCE_CONTEXT_SIMULATION_CONTEXT_BINDING_H
#define ABACUS_SOURCE_CONTEXT_SIMULATION_CONTEXT_BINDING_H

#ifndef ABACUS_CAN_BIND_SIMULATION_CONTEXT
#error "simulation_context_binding.h is restricted to the assembly layer"
#endif

namespace ModuleContext
{

struct SimulationContext;

class ScopedSimulationContextBinding
{
  public:
    explicit ScopedSimulationContextBinding(const SimulationContext& context);
    ~ScopedSimulationContextBinding();

    ScopedSimulationContextBinding(const ScopedSimulationContextBinding&) = delete;
    ScopedSimulationContextBinding& operator=(const ScopedSimulationContextBinding&) = delete;

  private:
    const SimulationContext* context_;
};

} // namespace ModuleContext

#endif
