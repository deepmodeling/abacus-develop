
#include "vdw.h"
#include "vdwd2.h"
#include "vdwd3.h"

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
        std::unique_ptr<Vdwd3> vdw_ptr = make_unique<Vdwd3>(ucell);
        vdw_ptr->parameter().initial_parameters(input, plog);
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