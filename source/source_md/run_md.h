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

/// Returns true when the MD loop should write a STRU frame at
/// `current_step`. A frame is produced on every md_restartfreq step (so
/// Restart_md.txt can always be resumed from a matching structure) and
/// additionally on output cadence steps when out_stru > 0. Effective output
/// frequency follows out_freq_ion with md_restartfreq fallback, matching
/// md_line(). This predicate exists so the "checkpoint frame tracks
/// Restart_md.txt" invariant stays unit-testable.
inline bool is_frame_write_step(
    const int current_step,
    const int out_freq_ion,
    const int md_restartfreq,
    const bool out_stru)
{
    const int eff_output_freq = (out_freq_ion > 0) ? out_freq_ion : md_restartfreq;
    const bool do_output_step = is_cadence_step(current_step, eff_output_freq);
    const bool do_restart_step = is_cadence_step(current_step, md_restartfreq);
    const bool want_output_frame = do_output_step && out_stru;
    const bool want_checkpoint_frame = do_restart_step;
    return want_output_frame || want_checkpoint_frame;
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
