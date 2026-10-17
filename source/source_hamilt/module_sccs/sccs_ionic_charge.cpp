#include "sccs_ionic_charge.h"

#include "source_base/constants.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/cell_tools.h"

#include <cmath>
#include <complex>

namespace ModuleSccs
{
void gaussian_ionic_density(const std::vector<unitcell::AtomData>& atoms,
                            const ModulePW::PW_Basis& basis,
                            double tpiba,
                            double spread,
                            std::vector<double>& density)
{
    std::vector<std::complex<double>> charge(basis.npw, 0.0);
    const double tpiba2 = tpiba * tpiba;
    for (const unitcell::AtomData& atom : atoms)
    {
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            const double exponent = -0.25 * spread * spread * tpiba2 * basis.gg[ig];
            const double phase = tpiba * (basis.gcar[ig] * atom.position);
            const std::complex<double> phase_factor = ModuleBase::NEG_IMAG_UNIT * phase;
            charge[ig] += atom.valence_charge / basis.omega * std::exp(exponent) * std::exp(phase_factor);
        }
    }
    density.resize(basis.nrxx);
    basis.recip2real(charge.data(), density.data());
}
} // namespace ModuleSccs
