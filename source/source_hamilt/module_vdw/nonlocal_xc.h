#ifndef NONLOCAL_XC_H
#define NONLOCAL_XC_H

#include <vector>

class Charge;
struct XCFunctionalParameters;
namespace ModuleBase
{
class matrix;
template <class T>
class Vector3;
} // namespace ModuleBase

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

/// The same, reusing the plane-wave gradient of the total density (valence plus core) that the caller
/// has already computed, which saves its four FFTs. nullptr computes the gradient here; the choice is
/// collective, so a process without grid points passes its empty vector, not nullptr.
void add_rvv10(const XCFunctionalParameters& parameters,
               const Charge& chr,
               const std::vector<ModuleBase::Vector3<double>>* total_gradient,
               double& etxc,
               double& vtxc,
               ModuleBase::matrix& v);

/// Add the rVV10 stress (Ry/Bohr^3) that the XC diagonal and potential terms do not contain.
void add_rvv10_stress(const XCFunctionalParameters& parameters, const Charge& chr, ModuleBase::matrix& sigma);

} // namespace vdw

#endif
