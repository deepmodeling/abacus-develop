#include "sccs_ionic_force.h"

#include "source_base/constants.h"
#include "source_base/parallel_reduce.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/cell_tools.h"

#include <cmath>
#include <complex>

namespace ModuleSccs
{
void gaussian_ionic_force(const std::vector<unitcell::AtomData>& atoms,
                          const std::vector<double>& reaction_potential,
                          const ModulePW::PW_Basis& basis,
                          double tpiba,
                          double spread,
                          std::vector<ModuleBase::Vector3<double>>& forces)
{
    std::vector<std::complex<double>> potential_g(basis.npw);
    basis.real2recip(reaction_potential.data(), potential_g.data());
    forces.assign(atoms.size(), ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
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
            forces[ia] += basis.gcar[ig] * factor;
        }
        Parallel_Reduce::reduce_pool(forces[ia].x);
        Parallel_Reduce::reduce_pool(forces[ia].y);
        Parallel_Reduce::reduce_pool(forces[ia].z);
    }
}
} // namespace ModuleSccs
