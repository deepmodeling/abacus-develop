#pragma once

namespace ModuleGint
{

enum class GintRealPrecision
{
    fp64,
    fp32
};

struct GintExecConfig
{
    GintRealPrecision cpu_internal_real = GintRealPrecision::fp64;
};

GintExecConfig current_exec_config();

class ScopedExecConfig
{
  public:
    explicit ScopedExecConfig(const GintExecConfig& cfg);
    ~ScopedExecConfig();

  private:
    GintExecConfig previous_cfg_;
};

} // namespace ModuleGint
