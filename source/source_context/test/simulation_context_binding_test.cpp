#include "source_context/orchestration_context.h"
#include "source_context/simulation_context_binding.h"
#include "source_hamilt/module_xc/exx_info.h"
#include "source_io/module_restart/restart.h"

#include <gtest/gtest.h>

#include <future>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace GlobalC
{
Exx_Info exx_info;
Restart restart;
} // namespace GlobalC

TEST(SimulationContextBinding, RejectsReadOutsideBindingLifetime)
{
    EXPECT_THROW(ModuleContext::current_simulation_context(), std::logic_error);

    ModuleContext::SimulationContext context;
    {
        ModuleContext::ScopedSimulationContextBinding binding(context);
        EXPECT_EQ(&context, &ModuleContext::current_simulation_context());
    }

    EXPECT_THROW(ModuleContext::current_simulation_context(), std::logic_error);
}

TEST(SimulationContextBinding, RejectsDuplicateBinding)
{
    ModuleContext::SimulationContext first;
    ModuleContext::SimulationContext second;
    ModuleContext::ScopedSimulationContextBinding binding(first);

    EXPECT_THROW(ModuleContext::ScopedSimulationContextBinding duplicate(second), std::logic_error);
    EXPECT_EQ(&first, &ModuleContext::current_simulation_context());
}

TEST(SimulationContextBinding, SupportsSequentialRuns)
{
    ModuleContext::SimulationContext first;
    ModuleContext::SimulationContext second;
    {
        ModuleContext::ScopedSimulationContextBinding binding(first);
        EXPECT_EQ(&first, &ModuleContext::current_simulation_context());
    }
    {
        ModuleContext::ScopedSimulationContextBinding binding(second);
        EXPECT_EQ(&second, &ModuleContext::current_simulation_context());
    }
}

TEST(SimulationContextBinding, SharesTheDriverOwnedObjectWithWorkerThreads)
{
    ModuleContext::SimulationContext context;
    ModuleContext::ScopedSimulationContextBinding binding(context);

    std::future<const ModuleContext::SimulationContext*> worker
        = std::async(std::launch::async, []() { return &ModuleContext::current_simulation_context(); });
    EXPECT_EQ(&context, worker.get());
}

#ifdef _OPENMP
TEST(SimulationContextBinding, IsVisibleToOpenMpWorkers)
{
    ModuleContext::SimulationContext context;
    ModuleContext::ScopedSimulationContextBinding binding(context);
    int matching_workers = 0;
    int worker_count = 0;
#pragma omp parallel reduction(+ : matching_workers, worker_count)
    {
        ++worker_count;
        if (&ModuleContext::current_simulation_context() == &context)
        {
            ++matching_workers;
        }
    }
    EXPECT_EQ(worker_count, matching_workers);
}
#endif
