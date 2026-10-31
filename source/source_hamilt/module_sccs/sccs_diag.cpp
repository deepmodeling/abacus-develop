#include "sccs_diag.h"
#include "sccs_response.h"

#include "source_base/constants.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <complex>

namespace ModuleSccs
{
void cont_pol_charge(const std::vector<double>& charge,
                     const SccsResponse& response,
                     const ModulePW::PW_Basis& basis,
                     double tpiba,
                     std::vector<double>& density)
{
    std::vector<std::complex<double>> potential_g(basis.npw);
    basis.real2recip(response.polarization.potential.data(), potential_g.data());
    std::vector<ModuleBase::Vector3<double>> gradient(charge.size());
    XC_Functional::grad_rho(potential_g.data(), gradient.data(), &basis, tpiba);
    density.resize(charge.size());
    for (std::size_t i = 0; i < charge.size(); ++i)
    {
        const double contraction = response.grad_log_epsilon[i] * gradient[i];
        density[i] = contraction / ModuleBase::FOUR_PI + charge[i] * (1.0 / response.epsilon[i] - 1.0);
    }
}
} // namespace ModuleSccs
