#ifndef SCCS_CAV_DERIV_H
#define SCCS_CAV_DERIV_H

#include "source_base/vector3.h"

#include <vector>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{
struct CavityParameters;
struct SccsResponse;

struct CavityDerivatives
{
    std::vector<double> coefficient;
    std::vector<double> filter;
    std::vector<ModuleBase::Vector3<double>> gradient;
};

// Collective over the PW pool. Periodic cells use density chain derivatives;
// open boundaries differentiate solute on the FFT grid.
void make_cav_derivs(const std::vector<double>& density,
                     const CavityParameters& cavity,
                     const ModulePW::PW_Basis& basis,
                     double tpiba,
                     bool open_boundary,
                     SccsResponse& response,
                     CavityDerivatives& derivatives);
} // namespace ModuleSccs

#endif
