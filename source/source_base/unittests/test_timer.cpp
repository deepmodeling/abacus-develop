#include "source_base/timer.h"

#include <gtest/gtest.h>
#include <stdexcept>

namespace
{
int timed_return(const bool early)
{
    ModuleBase::timer::Scope scope("TimerScopeTest", "timed_return");
    if (early)
    {
        return 1;
    }
    return 2;
}

void timed_exception()
{
    ModuleBase::timer::Scope scope("TimerScopeTest", "timed_exception");
    throw std::runtime_error("test exception");
}
} // namespace

TEST(TimerScope, EarlyReturnClosesAndAllowsReuse)
{
    EXPECT_EQ(timed_return(true), 1);
    const ModuleBase::timer::Timer_One& timer = ModuleBase::timer::timer_pool.at("TimerScopeTest").at("timed_return");
    EXPECT_TRUE(timer.start_flag);
    EXPECT_EQ(timer.calls, 1U);
    EXPECT_EQ(timed_return(false), 2);
    EXPECT_TRUE(timer.start_flag);
    EXPECT_EQ(timer.calls, 2U);
    ModuleBase::timer::timer_pool.erase("TimerScopeTest");
}

TEST(TimerScope, ExceptionClosesAndAllowsReuse)
{
    EXPECT_THROW(timed_exception(), std::runtime_error);
    const ModuleBase::timer::Timer_One& timer = ModuleBase::timer::timer_pool.at("TimerScopeTest").at("timed_exception");
    EXPECT_TRUE(timer.start_flag);
    EXPECT_EQ(timer.calls, 1U);
    EXPECT_THROW(timed_exception(), std::runtime_error);
    EXPECT_TRUE(timer.start_flag);
    EXPECT_EQ(timer.calls, 2U);
    ModuleBase::timer::timer_pool.erase("TimerScopeTest");
}
