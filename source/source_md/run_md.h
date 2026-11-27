#ifndef RUN_MD_H
#define RUN_MD_H

class MDCell;
class DomainDecomposition;
struct Parameter;

namespace ModuleESolver
{
class ESolver;
}

/**
 * @brief the md loop line
 *
 */
namespace Run_MD
{
/// Returns true when a per-step action with cadence frequency `freq` should
/// fire at `current_step`. A frequency <= 0 disables the action (never
/// fires). Keeping output and restart cadences behind this predicate makes
/// their independence directly unit-testable.
inline bool is_cadence_step(const int current_step, const int freq) { return freq > 0 && current_step % freq == 0; }

/// Periodic structure output is independent of the MD restart cadence.
inline bool is_structure_output_step(const int current_step, const int out_freq_ion, const bool out_stru)
{
    return out_stru && is_cadence_step(current_step, out_freq_ion);
}

/// Returns true when `current_step` is the final ionic step of the run.
/// current_step is the global step (local step + restart offset); md_nstep is
/// the global step count, so the last processed step is md_nstep - 1. The
/// final step is always written even if it does not sit on the output cadence.
inline bool is_last_step(const int current_step, const int md_nstep)
{
    return md_nstep > 0 && current_step == md_nstep - 1;
}

/**
 * @brief the md loop line
 *
 * @param cell cell information
 * @param p_esolver energy solver
 * @param md_para input parameters used in md
 */
void md_line(
    MDCell &mdcell,
    ModuleESolver::ESolver *p_esolver,
    const Parameter &param_in,
    DomainDecomposition &decomp);
} // namespace Run_MD

#endif
