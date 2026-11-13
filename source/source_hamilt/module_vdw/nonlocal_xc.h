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

// rVV10 nonlocal correlation of vdw_method rvv10, evaluated by libvdwxc. The three evaluation functions do
// nothing unless parameters.rvv10 is set. They use the total density of chr (valence plus core) and are
// collective over the pool communicator of chr.rhopw. The libvdwxc object, its FFT plans and its buffers
// are kept between calls (Rvv10Workspace) and released by finalize_rvv10.
// Without libvdwxc the four functions are empty, so that the callers need no preprocessor branch;
// vdw_method rvv10 then stops in make_vdw.

#ifdef __LIBVDWXC
/**
 * @brief Add the rVV10 energy and potential to the result of an XC evaluation.
 * @param parameters  rVV10 switch and the parameters b and C
 * @param chr         density; the valence density plus rho_core is used, distributed over chr.rhopw
 * @param[in,out] etxc  XC energy (Ry), already reduced over the pool
 * @param[in,out] vtxc  integral of v_xc n (Ry), already reduced over the pool
 * @param[in,out] v     XC potential (Ry) on the nrxx points of this process; for nspin = 4 only v(0) is changed
 */
void add_rvv10(const XCFunctionalParameters& parameters,
               const Charge& chr,
               double& etxc,
               double& vtxc,
               ModuleBase::matrix& v);

/**
 * @brief The same, reusing the plane-wave gradient of the total density that the caller already holds.
 * @param total_gradient  gradient of (valence plus core) density on the nrxx points of this process, which
 *                        saves the adapter's own four FFTs; nullptr computes it here. The choice is collective,
 *                        so a process without grid points passes its empty vector, not nullptr.
 */
void add_rvv10(const XCFunctionalParameters& parameters,
               const Charge& chr,
               const std::vector<ModuleBase::Vector3<double>>* total_gradient,
               double& etxc,
               double& vtxc,
               ModuleBase::matrix& v);

/**
 * @brief Add the rVV10 stress (Ry/Bohr^3) that the XC diagonal and potential terms do not contain.
 * @param parameters  rVV10 switch and the parameters b and C
 * @param chr         density, as in add_rvv10
 * @param[in,out] sigma  3 x 3 stress tensor
 */
void add_rvv10_stress(const XCFunctionalParameters& parameters, const Charge& chr, ModuleBase::matrix& sigma);

/// Release the libvdwxc object of the run. Call before MPI_Finalize; a later evaluation builds it again.
void finalize_rvv10();
#else
inline void add_rvv10(const XCFunctionalParameters&, const Charge&, double&, double&, ModuleBase::matrix&)
{
}
inline void add_rvv10(const XCFunctionalParameters&,
                      const Charge&,
                      const std::vector<ModuleBase::Vector3<double>>*,
                      double&,
                      double&,
                      ModuleBase::matrix&)
{
}
inline void add_rvv10_stress(const XCFunctionalParameters&, const Charge&, ModuleBase::matrix&)
{
}
inline void finalize_rvv10()
{
}
#endif

} // namespace vdw

#endif // NONLOCAL_XC_H
