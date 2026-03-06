#include "gint_precision.h"

#include <atomic>

namespace ModuleGint
{

namespace
{
std::atomic<GintRealPrecision> current_precision{GintRealPrecision::fp64};
}

GintExecConfig current_exec_config()
{
    return GintExecConfig{current_precision.load(std::memory_order_relaxed)};
}

ScopedExecConfig::ScopedExecConfig(const GintExecConfig& cfg)
    : previous_cfg_(current_exec_config())
{
    current_precision.store(cfg.cpu_internal_real, std::memory_order_relaxed);
}

ScopedExecConfig::~ScopedExecConfig()
{
    current_precision.store(previous_cfg_.cpu_internal_real, std::memory_order_relaxed);
}

} // namespace ModuleGint
