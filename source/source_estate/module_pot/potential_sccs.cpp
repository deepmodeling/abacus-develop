#include "potential_new.h"
#include "pot_sccs.h"

#include <ostream>

namespace elecstate
{
bool Potential::update_model(int electronic_iteration, double density_residual, std::ostream* log)
{
    PotSccs* sccs = this->sccs_component();
    if (sccs == nullptr)
    {
        return false;
    }
    const bool activated = sccs->try_activate(electronic_iteration, density_residual);
    if (activated && log != nullptr)
    {
        *log << " SCCS activated at electronic iteration " << electronic_iteration
             << ", DRHO = " << density_residual << std::endl;
    }
    return activated;
}

bool Potential::model_pending() const
{
    const PotSccs* sccs = this->sccs_component();
    return sccs != nullptr && !sccs->is_active();
}

void Potential::write_corr_iter(std::ostream& output, int level, double residual) const
{
    const PotSccs* sccs = this->sccs_component();
    if (level == 0 || sccs == nullptr)
    {
        return;
    }
    // Active SCCS reports the coupled PCC energy too; a deferred one its start condition.
    const double pcc_energy = this->pcc_energy_rydberg();
    sccs->write_iter(output, level, residual, pcc_energy);
}

const PotSccs* Potential::sccs_component() const
{
    for (const PotBase* component : this->components)
    {
        const PotSccs* sccs = dynamic_cast<const PotSccs*>(component);
        if (sccs != nullptr)
        {
            return sccs;
        }
    }
    return nullptr;
}

PotSccs* Potential::sccs_component()
{
    const Potential& self = *this;
    return const_cast<PotSccs*>(self.sccs_component());
}

SccsResume Potential::get_resume() const
{
    const PotSccs* sccs = this->sccs_component();
    if (sccs == nullptr)
    {
        return SccsResume();
    }
    return sccs->resume_state();
}
} // namespace elecstate
