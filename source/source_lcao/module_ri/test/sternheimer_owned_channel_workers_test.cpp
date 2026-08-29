#include "source_lcao/module_ri/sternheimer_channel_resources.h"

#include <cassert>

int main()
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::available,
                                                     1000000,
                                                     0,
                                                     1,
                                                     "test"};

    const auto sparse_plan = ModuleRI::plan_sternheimer_owned_channel_workers(36, 6, 30, 1, 0, memory);
    const auto bandwidth_plan = ModuleRI::plan_sternheimer_owned_channel_workers(163, 13, 30, 1, 0, memory);

    assert(sparse_plan.automatic_workers == 6);
    assert(sparse_plan.effective_workers == 1);
    assert(bandwidth_plan.automatic_workers == 13);
    assert(bandwidth_plan.effective_workers == 13);
    return 0;
}
