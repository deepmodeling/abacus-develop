
#include "vdw.h"
#include "vdwd2.h"
#include "vdwd3.h"

namespace vdw
{

std::unique_ptr<Vdw> make_vdw(const UnitCell &ucell, const Input &input)
{
    if (ucell.nat < 2 && input.vdw_method != "none")
    {
        ModuleBase::WARNING("VDW", "Only one atom in this system, and will not do the calculation of VDW");
        return nullptr;
    }
    else if (input.vdw_method == "d2")
    {
        std::unique_ptr<Vdwd2> vdw_ptr = std::make_unique<Vdwd2>(ucell);
        vdw_ptr->parameter().initial_parameters(input);
        vdw_ptr->parameter().initset(ucell);
        return vdw_ptr;
    }
    else if (input.vdw_method == "d3_0" || input.vdw_method == "d3_bj")
    {
        std::unique_ptr<Vdwd3> vdw_ptr = std::make_unique<Vdwd3>(ucell);
        vdw_ptr->parameter().initial_parameters(input);
        return vdw_ptr;
    }
    else
    {
        return nullptr;
    }
}

} // namespace vdw