#ifndef PCC_MOMENTS_H
#define PCC_MOMENTS_H

#include "pcc_0d.h"
#include "pcc_2d.h"

#include <vector>

namespace ModuleSurchem
{
class ChargeReduction;
}

namespace ModulePcc
{

// Sum local multipole moments over the processes that share the grid.
MultipoleMoments reduce_pcc_moments(MultipoleMoments moments,
                                    const ModuleSurchem::ChargeReduction& reduction);

Pcc2dMoments reduce_pcc_2d_moments(Pcc2dMoments moments,
                                   const ModuleSurchem::ChargeReduction& reduction);

// Moments of a distributed grid density about the PCC origin.
MultipoleMoments reduced_pcc_density_moments(
    const std::vector<double>& density,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    double volume_element,
    const PccGeometry& geometry,
    const ModuleSurchem::ChargeReduction& reduction);

Pcc2dMoments reduced_pcc_2d_density_moments(
    const std::vector<double>& density,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    double volume_element,
    const Pcc2dGeometry& geometry,
    const ModuleSurchem::ChargeReduction& reduction);

} // namespace ModulePcc

#endif
