#include "source_base/timer.h"

#include <gtest/gtest.h>
#include <stdexcept>

namespace
{
int timed_return(const bool early)
{
    ModuleBase::timer::Scope timer("TimerScopeTest", "timed_return");
    if (early)
    {
        return 1;
    }
    return 2;
}

void timed_exception()
{
    ModuleBase::timer::Scope timer("TimerScopeTest", "timed_exception");
    throw std::runtime_error("test exception");
}
} // namespace

TEST(TimerScope, EarlyReturnClosesAndAllowsReuse)
{
    ModuleBase::timer::enable();
    EXPECT_EQ(timed_return(true), 1);
    EXPECT_EQ(timed_return(false), 2);
    const auto& timer = ModuleBase::timer::timer_pool.at("TimerScopeTest").at("timed_return");
    EXPECT_TRUE(timer.start_flag);
    EXPECT_EQ(timer.calls, 2U);
}

TEST(TimerScope, ExceptionClosesAndAllowsReuse)
{
    ModuleBase::timer::enable();
    EXPECT_THROW(timed_exception(), std::runtime_error);
    EXPECT_THROW(timed_exception(), std::runtime_error);
    const auto& timer = ModuleBase::timer::timer_pool.at("TimerScopeTest").at("timed_exception");
    EXPECT_TRUE(timer.start_flag);
    EXPECT_EQ(timer.calls, 2U);
}
