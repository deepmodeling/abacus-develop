#ifndef ABACUS_LR_NACV_TIMER_H
#define ABACUS_LR_NACV_TIMER_H
#include "source_base/timer.h"
namespace LR
{
/** @brief Time a complete NACV operation, including exceptional exits. */
class NacvTimerScope
{
  public:
    /**
     * @param[in] category Owning class or module label, valid for this scope.
     * @param[in] name Enclosing function name, valid for this scope.
     */
    NacvTimerScope(const char *category, const char *name) : category_(category), name_(name)
    {
        ModuleBase::timer::start(category_, name_);
    }
    ~NacvTimerScope() { ModuleBase::timer::end(category_, name_); }
    NacvTimerScope(const NacvTimerScope &) = delete;
    NacvTimerScope &operator=(const NacvTimerScope &) = delete;

  private:
    const char *const category_;
    const char *const name_;
};
} // namespace LR
#endif
