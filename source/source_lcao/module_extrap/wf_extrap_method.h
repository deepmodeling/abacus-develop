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

} // namespace ModuleExtrap

#endif // SOURCE_LCAO_MODULE_EXTRAP_WF_EXTRAP_METHOD_H
