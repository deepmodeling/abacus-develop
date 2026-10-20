#include "sccs_fixed_src.h"
#include "source_base/timer.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_grid_geometry.h"
#include "source_cell/cell_tools.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/module_sccs/sccs_ion_charge.h"
#include "source_hamilt/module_sccs/sccs_parameters.h"

namespace
{
bool same_atoms(const std::vector<unitcell::AtomData>& left, const std::vector<unitcell::AtomData>& right)
{
    if (left.size() != right.size())
    {
        return false;
    }
    for (std::size_t i = 0; i < left.size(); ++i)
    {
        const ModuleBase::Vector3<double> displacement = left[i].position - right[i].position;
        if (displacement.norm2() != 0.0 || left[i].valence_charge != right[i].valence_charge
            || left[i].atomic_number != right[i].atomic_number)
        {
            return false;
        }
    }
    return true;
}
}
namespace elecstate
{
struct SccsFixedSources::Data
{
    bool valid = false;
    std::vector<unitcell::AtomData> atoms;
    std::vector<double> ionic_density;
    std::vector<double> core_density;
    std::vector<ModuleBase::Vector3<double>> positions;
};
SccsFixedSources::SccsFixedSources()
{
    data_.reset(new Data);
}
SccsFixedSources::~SccsFixedSources() = default;

bool SccsFixedSources::update(const UnitCell& cell,
                              const ModulePW::PW_Basis& basis,
                              const std::vector<unitcell::AtomData>& atoms,
                              const ModuleSccs::SccsConfig& config)
{
    ModuleBase::timer::start("SccsFixedSources", "update");
    // Atoms are replicated, so every pool rank takes the same decision.
    const bool reuse = data_->valid && same_atoms(data_->atoms, atoms);
    if (!reuse)
    {
        ModulePW::grid_positions(basis, cell.latvec, cell.lat0, data_->positions);
        ModuleSccs::gauss_ion_rho(atoms, basis, cell.tpiba, ModuleSccs::gaussian_ion_spread,
                                  data_->ionic_density);
        if (config.core_electrons)
        {
            ModuleSccs::gauss_core_rho(atoms, basis, cell.tpiba, config.core_spreads, data_->core_density);
        }
        data_->atoms = atoms;
        data_->valid = true;
    }
    ModuleBase::timer::end("SccsFixedSources", "update");
    return reuse;
}
const std::vector<double>& SccsFixedSources::ionic_density() const { return data_->ionic_density; }
const std::vector<double>& SccsFixedSources::core_density() const { return data_->core_density; }
const std::vector<ModuleBase::Vector3<double>>& SccsFixedSources::positions() const { return data_->positions; }
} // namespace elecstate
