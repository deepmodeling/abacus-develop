#ifndef XC_FUNCTIONAL_GGA_NONCOL_SF_BUILTIN_H
#define XC_FUNCTIONAL_GGA_NONCOL_SF_BUILTIN_H

#include "source_base/matrix.h"

#include <tuple>
#include <vector>

class Charge;
namespace ModulePW
{
class PW_Basis;
}
struct UnitCell;

namespace ModuleXC
{
namespace NCGGA_SF_Builtin
{

// gga_grad: 2 = projected divergence of h, 3 = full Scalmani-Frisch divergence
std::tuple<double, double, ModuleBase::matrix> v_xc_ncgga_sf_builtin(
    const int& nrxx, const double& omega, const double tpiba, const Charge* const chr,
    const int gga_grad);

void gradcorr_ncgga_sf_builtin(const Charge* const chr, ModulePW::PW_Basis* rhopw,
                                const UnitCell* ucell, std::vector<double>& stress_gga);

} // namespace NCGGA_SF_Builtin
} // namespace ModuleXC

#endif
