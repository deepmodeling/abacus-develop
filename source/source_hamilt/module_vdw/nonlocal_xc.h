#ifndef ABACUS_NONLOCAL_XC_H
#define ABACUS_NONLOCAL_XC_H

#include <memory>
#include <string>
#include <vector>

namespace ModulePW
{
class PW_Basis;
}

namespace vdw
{

/// rVV10 parameters published for the GGA_X_RPW86+GGA_C_PBE pairing
/// (R. Sabatini, T. Gorni and S. de Gironcoli, Phys. Rev. B 87, 041108(R) (2013)).
/// The current interface accepts only this pairing, so the values are fixed
/// here instead of being INPUT parameters.
const double RVV10_B_RPW86_PBE = 6.3;
const double RVV10_C_RPW86_PBE = 0.0093;

struct NonlocalXCResult
{
    double energy;
    double vtxc;
    std::vector<double> potential;
};

/**
 * @brief Density-grid interface for a nonlocal exchange-correlation backend.
 *
 * The host potential code only knows this interface. Concrete methods such as
 * rVV10, and future libvdwxc-backed methods, remain in the vdW module.
 */
class NonlocalXC
{
  public:
    virtual ~NonlocalXC() = default;

    /**
     * @param pw density plane-wave basis the densities are defined on
     * @param total_density valence plus core density on the real-space grid
     * @param valence_density valence density used for the vtxc contraction
     * @return nonlocal energy, vtxc and potential, all in Ry
     */
    virtual NonlocalXCResult evaluate(const ModulePW::PW_Basis& pw,
                                      const std::vector<double>& total_density,
                                      const std::vector<double>& valence_density) const = 0;
};

/// Create the backend for `method` ("rvv10"). Quits for unknown methods and
/// for builds without the required library.
std::unique_ptr<NonlocalXC> make_nonlocal_xc(const ModulePW::PW_Basis* pw, const std::string& method);

} // namespace vdw

#endif // ABACUS_NONLOCAL_XC_H
