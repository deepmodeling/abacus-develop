#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_LR_SFORCE_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_LR_SFORCE_H

#include "source_lcao/module_lr/lr_amp.h"
#include "source_lcao/module_lr/lr_mpl.h"
#include "source_lcao/module_lr/lr_zeq.h"

#include <vector>

namespace LR
{
/// Per-state side of the force path: the amplitudes widened into the Z window, plus the energies.
template <typename T>
struct StateCtx
{
    const std::vector<ct::Tensor>& X;
    const std::vector<Parallel_Blkcyclic>& para_x;
    const std::vector<Parallel_Blkcyclic>& para_x_z;
    const std::vector<int>& nocc;
    const std::vector<int>& nvirt;
    int nk;
    int nloc;
    int nloc_z;
    int nstates;
    bool openshell;
    const double* ekb;
};

/// Forces of the closed-shell states widened into the Z window; `istate_only < 0` means all.
template <typename T>
std::vector<ModuleBase::matrix> state_forces(const GradientInputs<T>& in, const GsForceCtx<T>& gs,
                                             const StateCtx<T>& st, int ispin, int istate_only);

/// Open-shell counterpart: both spin channels live in one vector, so the channel is always 0.
template <typename T>
std::vector<ModuleBase::matrix> state_forces_os(const GradientInputs<T>& in, const GsForceCtx<T>& gs,
                                                const StateCtx<T>& st, int istate_only);

/// One multiplet's gradient matrix, assembled from its diagonal forces.
template <typename T>
MplLvc build_mpl_state(const GradientInputs<T>& in, const GsForceCtx<T>& gs, const StateCtx<T>& st, int ispin,
                       const std::vector<int>& group, const std::vector<ModuleBase::matrix>& diag,
                       std::ofstream& ofs);

/// All multiplets of one channel among `forces`, with `degen_thr` as the grouping threshold.
template <typename T>
std::vector<MplLvc> multiplet_grad(const GradientInputs<T>& in, const GsForceCtx<T>& gs, const StateCtx<T>& st,
                                   int ispin, double degen_thr, const std::vector<ModuleBase::matrix>& forces,
                                   std::ofstream& ofs);

} // namespace LR

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_LR_SFORCE_H
