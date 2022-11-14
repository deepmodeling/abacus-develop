
#include "vdw.h"
#include "vdwd2.h"
#include "vdwd3.h"

namespace vdw
{

std::unique_ptr<Vdw> make_vdw(const UnitCell &ucell, const Input &input)
{
    if (INPUT.vdw_method == "d2")
    {
        std::unique_ptr<Vdwd2> vdw_ptr = make_unique<Vdwd2>(ucell);
        vdw_ptr->parameter().initial_parameters(input);
        vdw_ptr->parameter().initset(ucell);
        return vdw_ptr;
    }
    else if (INPUT.vdw_method == "d3_0" || INPUT.vdw_method == "d3_bj")
    {
        std::unique_ptr<Vdwd3> vdw_ptr = make_unique<Vdwd3>(ucell);
        vdw_ptr->parameter().initial_parameters(input);
        return vdw_ptr;
    }
    else
    {
        return nullptr;
    }
}

} // namespace vdw