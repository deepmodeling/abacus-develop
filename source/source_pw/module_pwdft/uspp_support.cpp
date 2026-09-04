#include "uspp_support.h"

#include "source_base/tool_quit.h"

#include <sstream>
#include <vector>

namespace pw
{

void validate_uspp_support(const bool use_uspp,
                           const std::string& basis_type,
                           const std::string& esolver_type,
                           const int nspin,
                           const int xc_func_type,
                           const bool berry_phase,
                           const bool towannier90,
                           const bool cal_cond)
{
    if (!use_uspp)
    {
        return;
    }

    std::vector<std::string> violations;
    if (basis_type != "pw")
    {
        violations.push_back("basis_type=" + basis_type + " (only pw is supported)");
    }
    if (esolver_type != "ksdft")
    {
        violations.push_back("esolver_type=" + esolver_type + " (only ksdft is supported)");
    }
    if (nspin != 1 && nspin != 2)
    {
        violations.push_back("nspin=" + std::to_string(nspin) + " (only 1 and 2 are supported)");
    }
    if (xc_func_type != 1 && xc_func_type != 2)
    {
        violations.push_back("XC functional type=" + std::to_string(xc_func_type) + " (only LDA and GGA are supported)");
    }
    if (berry_phase)
    {
        violations.push_back("berry_phase=true is not supported");
    }
    if (towannier90)
    {
        violations.push_back("towannier90=true is not supported");
    }
    if (cal_cond)
    {
        violations.push_back("cal_cond=true is not supported");
    }

    if (violations.empty())
    {
        return;
    }

    std::ostringstream message;
    message << "Unsupported USPP configuration: ";
    for (std::size_t index = 0; index < violations.size(); ++index)
    {
        if (index != 0)
        {
            message << "; ";
        }
        message << violations[index];
    }
    message << ".";
    ModuleBase::WARNING_QUIT("pw::validate_uspp_support", message.str());
}

} // namespace pw
