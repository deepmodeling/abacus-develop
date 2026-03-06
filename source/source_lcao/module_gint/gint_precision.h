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

} // namespace ModuleGint
