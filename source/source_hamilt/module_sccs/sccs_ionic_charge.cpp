#include "sccs_ionic_charge.h"
#include "sccs_pw_coulomb.h"

#include "source_base/constants.h"
#include "source_base/parallel_reduce.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/cell_tools.h"

#include <cmath>
#include <complex>

namespace ModuleSccs
{
bool gaussian_ionic_density(const std::vector<unitcell::AtomData>& atoms,
                            const ModulePW::PW_Basis& basis,
                            double tpiba,
                            double spread,
                            std::vector<double>& density,
                            std::string& error)
{
    if (!validate_pw_grid(basis, tpiba, error))
    {
        return false;
    }
    double invalid = 0.0;
    if (!std::isfinite(spread) || spread <= 0.0)
    {
        invalid = 1.0;
    }
    for (const unitcell::AtomData& atom : atoms)
    {
        if (!std::isfinite(atom.valence_charge) || atom.valence_charge < 0.0
            || !std::isfinite(atom.position.x) || !std::isfinite(atom.position.y)
            || !std::isfinite(atom.position.z))
        {
            invalid = 1.0;
        }
    }
    Parallel_Reduce::reduce_max_pool(basis.poolnproc, invalid);
    if (invalid != 0.0)
    {
        error = "SCCS Gaussian ions require finite positions, nonnegative valence charges and positive finite spread";
        return false;
    }

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
    std::vector<double> candidate(basis.nrxx);
    basis.recip2real(charge.data(), candidate.data());
    if (!validate_grid_values(candidate, basis, error))
    {
        return false;
    }
    density.swap(candidate);
    return true;
}
} // namespace ModuleSccs
