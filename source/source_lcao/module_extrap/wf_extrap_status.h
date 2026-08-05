#ifndef SOURCE_LCAO_MODULE_EXTRAP_WF_EXTRAP_STATUS_H
#define SOURCE_LCAO_MODULE_EXTRAP_WF_EXTRAP_STATUS_H

namespace ModuleExtrap
{

enum class WfcExtrapStatus
{
    Success,
    EmptyHistory,
    DimensionMismatch,
    InvalidInput,
    Unsupported,
    OrthogonalizationFailed
};

inline const char* to_string(const WfcExtrapStatus status) noexcept
{
    switch (status)
    {
    case WfcExtrapStatus::Success:
        return "success";
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

#endif // SOURCE_LCAO_MODULE_EXTRAP_WF_EXTRAP_STATUS_H
