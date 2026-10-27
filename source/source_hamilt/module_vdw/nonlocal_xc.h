#ifndef NONLOCAL_XC_H
#define NONLOCAL_XC_H

#include <tuple>

class Charge;
namespace ModuleBase
{
class matrix;
}

namespace vdw
{

/**
 * @brief Add the rVV10 nonlocal correlation, evaluated by libvdwxc, to the results of XC_Functional::v_xc.
 *
 * The energy and the potential are those of the total density of chr, valence plus core. The
 * potential is added to every spin channel that carries density, and vtxc counts the valence
 * density only, as in v_xc. b and c are the rVV10 parameters; all results are in Ry. The complete
 * grid is collected on the root of the k-point pool, where libvdwxc runs serially.
 */
void add_rvv10(const double b,
               const double c,
               const Charge& chr,
               std::tuple<double, double, ModuleBase::matrix>& etxc_vtxc_v);

} // namespace vdw

#endif
