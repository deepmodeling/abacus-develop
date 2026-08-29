#include "source_lcao/module_ri/sternheimer_channel_resources.h"

#include <cassert>

int main()
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::available,
                                                     1000000,
                                                     0,
                                                     1,
                                                     "test"};

    const auto plan = ModuleRI::plan_sternheimer_owned_channel_workers(163, 9, 30, 1, 0, memory);

    assert(plan.automatic_workers == 9);
    assert(plan.effective_workers == 1);
    return 0;
}
