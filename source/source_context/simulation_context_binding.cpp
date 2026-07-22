#include "source_context/orchestration_context.h"
#include "source_context/simulation_context_binding.h"

#include <atomic>
#include <stdexcept>

namespace
{

std::atomic<const ModuleContext::SimulationContext*> active_context(nullptr);

} // namespace

namespace ModuleContext
{

ScopedSimulationContextBinding::ScopedSimulationContextBinding(const SimulationContext& context) : context_(&context)
{
    const SimulationContext* expected = nullptr;
    if (!active_context.compare_exchange_strong(expected,
                                                context_,
                                                std::memory_order_release,
                                                std::memory_order_relaxed))
    {
        throw std::logic_error("a SimulationContext is already bound");
    }
}

ScopedSimulationContextBinding::~ScopedSimulationContextBinding()
{
    const SimulationContext* expected = context_;
    active_context.compare_exchange_strong(expected,
                                           nullptr,
                                           std::memory_order_release,
                                           std::memory_order_relaxed);
}

const SimulationContext& current_simulation_context()
{
    const SimulationContext* context = active_context.load(std::memory_order_acquire);
    if (context == nullptr)
    {
        throw std::logic_error("SimulationContext has not been bound");
    }
    return *context;
}

} // namespace ModuleContext
