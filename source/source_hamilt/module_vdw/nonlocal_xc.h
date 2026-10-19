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

    virtual NonlocalXCResult evaluate(const ModulePW::PW_Basis& pw,
                                      const std::vector<double>& total_density,
                                      const std::vector<double>& valence_density) const = 0;
};

std::unique_ptr<NonlocalXC> make_nonlocal_xc(const ModulePW::PW_Basis* pw,
                                             const std::string& method,
                                             double b,
                                             double c);

} // namespace vdw

#endif // ABACUS_NONLOCAL_XC_H
