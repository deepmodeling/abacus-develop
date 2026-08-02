#ifndef XC_FUNCTIONAL_GGA_NONCOL_SF_BUILTIN_H
#define XC_FUNCTIONAL_GGA_NONCOL_SF_BUILTIN_H

#include "source_base/matrix.h"
#include "source_estate/module_charge/charge.h"

#include <tuple>
#include <vector>

namespace ModulePW
{
class PW_Basis;
}
struct UnitCell;

namespace ModuleXC
{
namespace NCGGA_SF_Builtin
{

std::tuple<double, double, ModuleBase::matrix> v_xc_ncgga_sf_builtin(
    const int& nrxx, const double& omega, const double tpiba, const Charge* const chr);

void gradcorr_ncgga_sf_builtin(const Charge* const chr, ModulePW::PW_Basis* rhopw,
                                const UnitCell* ucell, std::vector<double>& stress_gga);

} // namespace NCGGA_SF_Builtin
} // namespace ModuleXC

#endif
