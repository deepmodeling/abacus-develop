#ifndef XC_FUNCTIONAL_NCGGA_SF_H
#define XC_FUNCTIONAL_NCGGA_SF_H

#include "source_base/matrix.h"

#include <tuple>
#include <vector>

class Charge;
namespace ModulePW
{
class PW_Basis;
}

namespace ModuleXC
{
namespace NCGGA_SF_Builtin
{

// Exact discrete reverse of the regularized projected LCA graph (gga_grad=2).
std::tuple<double, double, ModuleBase::matrix> v_xc_ncgga_sf_builtin(const int& nrxx,
                                                                     const double& omega,
                                                                     const double tpiba,
                                                                     const Charge* const chr);

// Gradient-metric stress of the exact gga_grad=2 projected-LCA graph.  The
// returned lower triangle is the unnormalised real-grid sum; Stress_Func
// applies the existing pool reduction and 1/nxyz normalisation.
void gradcorr_ncgga_lca_builtin(const Charge* const chr,
                                ModulePW::PW_Basis* rhopw,
                                const double tpiba,
                                std::vector<double>& stress_gga);

} // namespace NCGGA_SF_Builtin
} // namespace ModuleXC

#endif
