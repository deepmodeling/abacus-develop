#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_LR_GS_POT_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_LR_GS_POT_H

#include <memory>
#include <string>

class Charge;
class Parallel_Grid;
class Structure_Factor;
class UnitCell;
class pseudopot_cell_vl;
class surchem;
struct Input_para;
namespace ModulePW
{
class PW_Basis;
}
namespace elecstate
{
class Potential;
}

namespace LR
{
/// Owns the ground-state potential used by the LR force path.
class GsPot
{
  public:
    GsPot();
    ~GsPot();

    GsPot(const GsPot&) = delete;
    GsPot& operator=(const GsPot&) = delete;

    /// Rebuild the potential on `chg`. `refresh_sfac` is true on the file (`lr`) path, where the
    /// solver owns `sf`/`vloc`; on `ks-lr` the ground-state solver has already refreshed them.
    void refresh(const Input_para& inp, const std::string& xc_kernel, const UnitCell& ucell,
                 ModulePW::PW_Basis* pw_rhod, ModulePW::PW_Basis* pw_rho, const Parallel_Grid& pgrid, const Charge& chg,
                 Structure_Factor& sf, pseudopot_cell_vl& vloc, surchem& solvent, bool refresh_sfac,
                 bool has_float_data);

    elecstate::Potential* pot_hartree() const;

  private:
    std::unique_ptr<elecstate::Potential> pot_;
    std::unique_ptr<elecstate::Potential> pot_hartree_; ///< Hartree only, used by the test_force path
    double etxc_ = 0.0;
    double vtxc_ = 0.0;
};

} // namespace LR

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_LR_GS_POT_H
