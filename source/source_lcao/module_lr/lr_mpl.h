#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_LR_MPL_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_LR_MPL_H

#include "source_lcao/module_lr/utils/lr_util.h"

#include <functional>
#include <fstream>
#include <vector>

class UnitCell;

namespace LR
{
/// Linear vibronic coupling data of one degenerate multiplet.
struct MplLvc
{
    int ispin = 0;             ///< which spin channel (index into `spin_types`)
    std::vector<int> states;   ///< the multiplet's state indices, ascending
    double omega0 = 0.0;       ///< the common excitation energy (Ry)
    double omega_spread = 0.0; ///< max - min over the members (Ry); 0 if exact
    /// G[k][l], symmetric, each a (nat, 3) force matrix
    std::vector<std::vector<ModuleBase::matrix>> g;
};

/// Evaluate the forces of a batch of Z-window vectors: (vectors, omega, label) -> forces.
/// The caller owns the Z-vector solver, which is the only solver-side ingredient left.
using MplEval = std::function<std::vector<ModuleBase::matrix>(const ct::Tensor&, const std::vector<double>&, int)>;

/// @brief Assemble the gradient matrix of one multiplet from its diagonal forces.
///
/// Wraps the polarization identity of `grad_degen.h`: normalize the pairs, evaluate them in
/// one batch, and combine with the diagonal.
/// @param xz     the multiplet's members widened into the Z window, one block each
/// @param omega  their excitation energies, aligned with `states`
/// @param diag   their per-state forces, i.e. the diagonal of the matrix
/// @param eval   evaluates the forces of arbitrary vectors in the same window
template <typename T>
MplLvc build_mpl(const ct::Tensor& xz, const std::vector<double>& omega,
                 const std::vector<ModuleBase::matrix>& diag, int nloc, int ispin,
                 const std::vector<int>& states, const UnitCell& ucell, const MplEval& eval, std::ofstream& log);

/// Force of the steepest-descending branch of a multiplet, plus its mixing.
struct JtForce
{
    ModuleBase::matrix force;
    std::vector<double> mixing;
};

/// @brief Descend the Jahn-Teller branch: find the direction, log the split, return its force.
JtForce eval_jt_force(const std::vector<std::vector<ModuleBase::matrix>>& g, std::ostream& log);

} // namespace LR

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_LR_MPL_H
