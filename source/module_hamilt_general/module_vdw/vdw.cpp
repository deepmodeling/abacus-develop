#include <algorithm>
#include <cassert>

#include "vdw.h"
#include "vdwd2.h"
#include "vdwd3.h"
#include "module_base/tool_quit.h"

std::string parse_xcname(const std::string &xc_input,
                         const std::vector<std::string> &xc_psp)
{
    if (xc_input != "default")
    {
        return xc_input;
    }

    if (xc_psp.size() <= 0)
    {
        ModuleBase::WARNING_QUIT("ModuleHamiltGeneral::ModuleVDW::parse_xcname", 
        "XC name automatic inference failed: no pseudopotential files are found");
    }
    std::vector<std::string> xc_psp_uniq = xc_psp;
    std::sort(xc_psp_uniq.begin(), xc_psp_uniq.end());
    auto last = std::unique(xc_psp_uniq.begin(), xc_psp_uniq.end());
    xc_psp_uniq.erase(last, xc_psp_uniq.end());
    
    if (xc_psp_uniq.size() > 1)
    {
        ModuleBase::WARNING_QUIT("ModuleHamiltGeneral::ModuleVDW::parse_xcname", 
        "XC name automatic inference failed: inconsistency in XC names is found"
        " in the pseudopotential files");
    }
    const std::string xc = xc_psp_uniq[0];
    std::cout << " ***WARNING*** ModuleHamiltGeneral::ModuleVDW::parse_xcname: "
              << "XC name is automatically inferred from pseudopotential as `" 
              << xc << "`" << std::endl;
    return xc;
}

namespace vdw
{

std::unique_ptr<Vdw> make_vdw(const UnitCell &ucell, 
                              const Input_para &input,
                              std::ofstream* plog)
{
    // if (ucell.nat < 2 && input.vdw_method != "none")
    // {
    //     ModuleBase::WARNING("VDW", "Only one atom in this system, and will not do the calculation of VDW");
    //     return nullptr;
    // }
    if (input.vdw_method == "d2")
    {
        std::unique_ptr<Vdwd2> vdw_ptr = make_unique<Vdwd2>(ucell);
        vdw_ptr->parameter().initial_parameters(input, plog);
        vdw_ptr->parameter().initset(ucell);
        return vdw_ptr;
    }
    else if (input.vdw_method == "d3_0" || input.vdw_method == "d3_bj")
    {
        std::vector<std::string> xc_psp(ucell.ntype);
        for (int it = 0; it < ucell.ntype; it++)
        {
            xc_psp[it] = ucell.atoms[it].ncpp.xc_func;
        }
        std::unique_ptr<Vdwd3> vdw_ptr = make_unique<Vdwd3>(ucell);
        vdw_ptr->parameter().initial_parameters(parse_xcname(input.dft_functional, xc_psp), input, plog);
        return vdw_ptr;
    }
    else if (input.vdw_method != "none")
    {
        ModuleBase::WARNING_QUIT("ModuleHamiltGeneral::ModuleVDW::make_vdw", 
        "Unrecognized Van der Waals correction method: " + input.vdw_method);
        return nullptr;
    }
    return nullptr; // "none" method
}

} // namespace vdw