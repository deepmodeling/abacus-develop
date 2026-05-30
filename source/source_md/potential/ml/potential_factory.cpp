#include "potential_factory.h"

#include "ml/dpmd_potential.h"
#include "ml/nep_potential.h"

#include <stdexcept>

namespace ModuleMD
{

std::unique_ptr<PotentialBase> PotentialFactory::create(
    const std::string& potential_type,
    const std::string& model_file)
{
    std::unique_ptr<PotentialBase> potential;

    if (potential_type == "dpmd" || potential_type == "deepmd")
    {
        potential = std::make_unique<DPMDPotential>();
    }
    else if (potential_type == "nep")
    {
        potential = std::make_unique<NEPPotential>();
    }
    else
    {
        throw std::runtime_error(
            "Unknown MD potential type: " + potential_type);
    }

    potential->init(model_file);
    return potential;
}

} // namespace ModuleMD