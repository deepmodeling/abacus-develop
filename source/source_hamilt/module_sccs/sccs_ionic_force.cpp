#include "sccs_ionic_force.h"
#include "sccs_pw_coulomb.h"

#include "source_base/constants.h"
#include "source_base/parallel_reduce.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/cell_tools.h"

#include <cmath>
#include <complex>

namespace ModuleSccs
{
bool gaussian_ionic_force(const std::vector<unitcell::AtomData>& atoms,
                          const std::vector<double>& reaction_potential,
                          const ModulePW::PW_Basis& basis,
                          double tpiba,
                          double spread,
                          std::vector<ModuleBase::Vector3<double>>& forces,
                          std::string& error)
{
    if (!validate_pw_grid(basis, tpiba, error)
        || !validate_grid_values(reaction_potential, basis, error))
    {
        return false;
    }
    double invalid = 0.0;
    if (!std::isfinite(spread) || spread <= 0.0) { invalid = 1.0; }
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
        error = "SCCS Gaussian forces require finite atoms and a positive finite spread";
        return false;
    }
    std::vector<std::complex<double>> potential_g(basis.npw);
    basis.real2recip(reaction_potential.data(), potential_g.data());
    std::vector<ModuleBase::Vector3<double>> candidate(atoms.size());
    const double tpiba2 = tpiba * tpiba;
    // F = -integral(v_reaction d rho_ion/d R), using ABACUS's normalized FFT.
    // The cell volume cancels the Z/omega in the Gaussian charge coefficients.
    for (std::size_t ia = 0; ia < atoms.size(); ++ia)
    {
        const unitcell::AtomData& atom = atoms[ia];
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            const double exponent = -0.25 * spread * spread * tpiba2 * basis.gg[ig];
            const double gaussian = std::exp(exponent);
            const double phase = tpiba * (basis.gcar[ig] * atom.position);
            const std::complex<double> phase_argument = ModuleBase::NEG_IMAG_UNIT * phase;
            const std::complex<double> phase_factor = std::exp(phase_argument);
            const std::complex<double> weighted = std::conj(potential_g[ig]) * phase_factor;
            const double factor = -atom.valence_charge * tpiba * gaussian * weighted.imag();
            candidate[ia] += basis.gcar[ig] * factor;
        }
        Parallel_Reduce::reduce_pool(candidate[ia].x);
        Parallel_Reduce::reduce_pool(candidate[ia].y);
        Parallel_Reduce::reduce_pool(candidate[ia].z);
        if (!std::isfinite(candidate[ia].x) || !std::isfinite(candidate[ia].y)
            || !std::isfinite(candidate[ia].z))
        {
            error = "SCCS Gaussian ionic force is not finite";
            return false;
        }
    }
    forces.swap(candidate);
    return true;
}
} // namespace ModuleSccs
