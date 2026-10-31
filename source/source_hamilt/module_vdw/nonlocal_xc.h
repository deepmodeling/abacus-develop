#ifndef NONLOCAL_XC_H
#define NONLOCAL_XC_H

class Charge;
struct XCFunctionalParameters;
namespace ModuleBase
{
class matrix;
}

namespace vdw
{

// rVV10 nonlocal correlation of vdw_method rvv10, evaluated by libvdwxc. Both functions do nothing
// unless parameters.rvv10 is set. They use the total density of chr (valence plus core) and are
// collective over the pool communicator of chr.rhopw. The libvdwxc object, its FFT plans and its
// buffers are kept between calls.

/// Add the rVV10 energy, potential and vtxc (Ry) to the pool-reduced results of an XC evaluation.
void add_rvv10(const XCFunctionalParameters& parameters,
               const Charge& chr,
               double& etxc,
               double& vtxc,
               ModuleBase::matrix& v);

/// Add the rVV10 stress (Ry/Bohr^3) that the XC diagonal and potential terms do not contain.
void add_rvv10_stress(const XCFunctionalParameters& parameters, const Charge& chr, ModuleBase::matrix& sigma);

} // namespace vdw

#endif
