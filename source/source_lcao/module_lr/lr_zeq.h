#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_LR_ZEQ_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_LR_ZEQ_H

#include "source_lcao/module_lr/gradient_inputs.h"

#include <string>
#include <vector>

namespace LR
{
/// The ground-state side of a force evaluation: density-matrix ingredients and force terms.
/// Kept separate from `GradientInputs`, which describes the Z-window side.
template <typename T>
struct GsForceCtx
{
    const psi::Psi<T>& psi_all;
    const ModuleBase::matrix& wg_all;
    const Parallel_Orbitals& pv_all;
    const ModulePW::PW_Basis& rho;
    const ModulePW::PW_Basis& rhod;
    const pseudopot_cell_vl& vloc;
    Structure_Factor& sfac;
    const TwoCenterBundle& tcb;
    ForceSettings settings;
    bool openshell;
    const std::string& grad_solver;
    const std::string& in_dir;
};

/// Solve the Z-vector equation for `nst` independent blocks of `xz`.
template <typename T>
ct::Tensor solve_zeq(const GradientInputs<T>& in, const ct::Tensor& xz, int ispin, int nst, bool openshell,
                     const std::string& grad_solver, const std::string& in_dir,
                     const std::vector<std::vector<ct::Tensor>>& window_density);

/// Forces of `nst` closed-shell states widened into the Z window.
template <typename T>
std::vector<ModuleBase::matrix> force_xz(const GradientInputs<T>& in, const GsForceCtx<T>& gs, int ispin,
                                         const ct::Tensor& xz, const std::vector<double>& omega, int label_begin);

/// Open-shell counterpart: both spin channels live in one vector, so the channel is always 0.
template <typename T>
std::vector<ModuleBase::matrix> force_xz_os(const GradientInputs<T>& in, const GsForceCtx<T>& gs, const ct::Tensor& xz,
                                            const std::vector<double>& omega, int label_begin);

} // namespace LR

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_LR_ZEQ_H
