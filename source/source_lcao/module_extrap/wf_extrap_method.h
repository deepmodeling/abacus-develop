#ifndef SOURCE_LCAO_MODULE_EXTRAP_WF_EXTRAP_METHOD_H
#define SOURCE_LCAO_MODULE_EXTRAP_WF_EXTRAP_METHOD_H

#include <string>

namespace ModuleExtrap
{

enum class WfcExtrapMethod
{
    None,
    UsePrevWf
};

enum class WfcExtrapStatus
{
    Success,
    Disabled,
    EmptyHistory,
    DimensionMismatch,
    InvalidInput,
    Unsupported,
    OrthogonalizationFailed
};

inline WfcExtrapMethod wfc_extrap_method_from_string(const std::string& method)
{
    if (method == "use_prev_wf")
    {
        return WfcExtrapMethod::UsePrevWf;
    }
    return WfcExtrapMethod::None;
}

inline const char* to_string(const WfcExtrapMethod method) noexcept
{
    switch (method)
    {
    case WfcExtrapMethod::UsePrevWf:
        return "use_prev_wf";
    case WfcExtrapMethod::None:
    default:
        return "none";
    }
}

inline const char* to_string(const WfcExtrapStatus status) noexcept
{
    switch (status)
    {
    case WfcExtrapStatus::Success:
        return "success";
    case WfcExtrapStatus::Disabled:
        return "disabled";
    case WfcExtrapStatus::EmptyHistory:
        return "empty_history";
    case WfcExtrapStatus::DimensionMismatch:
        return "dimension_mismatch";
    case WfcExtrapStatus::InvalidInput:
        return "invalid_input";
    case WfcExtrapStatus::Unsupported:
        return "unsupported";
    case WfcExtrapStatus::OrthogonalizationFailed:
        return "orthogonalization_failed";
    default:
        return "unknown";
    }
}

} // namespace ModuleExtrap

#endif // SOURCE_LCAO_MODULE_EXTRAP_WF_EXTRAP_METHOD_H
