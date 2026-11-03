#ifndef SCCS_FIXED_SRC_H
#define SCCS_FIXED_SRC_H

#include <memory>
#include <vector>

class UnitCell;
namespace ModuleBase { template <typename T> class Vector3; }
namespace ModulePW { class PW_Basis; }
namespace unitcell { struct AtomData; }
namespace ModuleSccs { struct SccsConfig; }
namespace elecstate
{
// Owns ion/core Gaussian densities and Cartesian grid coordinates. The cell,
// grid and SCCS configuration are fixed for one PotSccs, which is rebuilt every
// ionic step; update rebuilds only for different atoms and returns true when
// the previous sources were reused.
class SccsFixedSources
{
public:
    SccsFixedSources();
    ~SccsFixedSources();
    bool update(const UnitCell& cell,
                const ModulePW::PW_Basis& basis,
                const std::vector<unitcell::AtomData>& atoms,
                const ModuleSccs::SccsConfig& config);
    // Empty until the first update; views live until the next update.
    const std::vector<double>& ionic_density() const;
    const std::vector<double>& core_density() const;
    const std::vector<ModuleBase::Vector3<double>>& positions() const;
private:
    struct Data;
    std::unique_ptr<Data> data_;
};
} // namespace elecstate
#endif
