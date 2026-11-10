#ifndef RUN_MD_H
#define RUN_MD_H

class MDCell;
class DomainDecomposition;
class UnitCell;
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
bool is_cadence_step(int current_step, int freq);

/// Returns true when the MD loop should write a STRU frame at
/// `current_step`. A frame is produced on every md_restartfreq step (so
/// Restart_md.txt can always be resumed from a matching structure) and
/// additionally on output cadence steps when out_stru > 0. Effective output
/// frequency follows out_freq_ion with md_restartfreq fallback, matching
/// md_line(). This predicate exists so the "checkpoint frame tracks
/// Restart_md.txt" invariant stays unit-testable.
bool is_frame_write_step(int current_step, int out_freq_ion, int md_restartfreq, bool out_stru);

void prepare_mdcell(MDCell& mdcell, const Parameter& param_in, DomainDecomposition& decomp);

void prepare_mdcell(MDCell& mdcell, UnitCell& ucell, DomainDecomposition& decomp);

/**
 * @brief the md loop line
 *
 * @param cell cell information
 * @param p_esolver energy solver
 * @param md_para input parameters used in md
 */
void md_line(MDCell& mdcell,
             ModuleESolver::ESolver* p_esolver,
             const Parameter& param_in,
             DomainDecomposition& decomp);
} // namespace Run_MD

#endif
